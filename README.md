# **HULK** — High-Volume Bulk RNA-seq Pipeline

**HULK** is a fully automated, containerized pipeline for bulk RNA-seq. It takes data either from SRA accessions or from local FASTQ files and produces gene-level quantification, QC reports, and R-based post-processing outputs (matrices and plots) with a single command.

> Fetch / Import → QC/Trim → Quantify → Post-process → Report

---

## Table of Contents

- [Features](#features)
- [Supported Platforms](#supported-platforms)
- [Installation](#installation)
  - [Option 1 — Official installer (recommended)](#option-1--official-installer-recommended)
  - [Option 2 — Run via Docker (no install)](#option-2--run-via-docker-no-install)
- [Quickstart](#quickstart)
- [Input Table (SRA RunInfo)](#input-table-sra-runinfo)
- [FASTQ Mode](#fastq-mode)
- [CLI Reference](#cli-reference)
- [Usage Examples](#usage-examples)
- [Output Structure](#output-structure)
- [Tools used inside HULK](#tools-used-inside-hulk)
- [License](#license)
- [Citation](#citation)
- [Links](#links)


## Features

- **End-to-end automation**
  - Fetches SRA runs (`prefetch`) **or** imports local FASTQ files (FASTQ mode).
  - Converts SRA accessions to FASTQ (`fasterq-dump`) in SRA mode.
  - Performs QC and trimming with `fastp`.
  - Quantifies transcript expression with `kallisto`.
  - Generates per-project and global `MultiQC` reports.

- **Two input modes**
  - **SRA mode**: runs are driven by an SRA RunInfo table with Run and BioProject information.
  - **FASTQ mode**: runs directly on user-provided FASTQ files grouped in a directory, sharing the same trimming, quantification and post-processing steps.

- **R-based post-processing**
  - Uses **tximport** to import quantification results into R.
  - Uses **DESeq2** (optional) to compute:
    - Normalized gene-level count matrices.
    - Variance-stabilizing transform (VST) expression matrices.  
    (DESeq2 is currently run without metadata/design formulas.)
  - Produces:
    - Variance heatmaps.
    - Expression heatmaps.
    - PCA plots at the sample level.
  - Saves both the figures and the corresponding matrix files for downstream reuse.

- **Gene-level expression focus**
  - Outputs are standardized on **gene-level counts**.
  - When starting from transcript-level quantifications, HULK aggregates transcripts to genes using a `tx2gene` mapping file provided by the user.
  - Gene-level matrices are designed to feed directly into downstream statistical and network analysis tools.

- **BioProject- and sample-aware outputs**
  - Per-BioProject plots and DESeq2 outputs in dedicated folders.
  - Per-BioProject log files summarizing what happened for each project.
  - Per-sample log files to simplify debugging of individual runs.
  - A structured directory layout that groups results by BioProject and by sample.

- **Concurrency and performance**
  - Separates **prefetch** and **processing** workers in SRA mode:
    - A dedicated prefetch loop continuously downloads SRR data into a shared cache.
    - Processing workers operate only on already-prefetched data, reducing idle time and improving throughput.
  - Uses CPU-aware concurrency to distribute work across BioProjects and runs.

- **Containerized & reproducible**
  - Distributed as a Docker image.
  - No local installation of tools is required beyond Docker.
  - Designed to be resumable and safe to re-run in partially processed directories.

---

## Supported Platforms

| Vendor / Family        | Models                                                                 |
|------------------------|-------------------------------------------------------------------------|
| **Illumina NovaSeq**   | 6000, X, X Plus                                                         |
| **Illumina HiSeq**     | X Ten, 4000, 3000, 2500, 2000, 1500, 1000                              |
| **Illumina NextSeq**   | 1000, 500, 550                                                          |
| **BGI / MGI**          | DNBSEQ-G400, DNBSEQ-T7, BGISEQ-500, MGISEQ-2000RS                      |
| **Legacy Illumina**    | Genome Analyzer II / IIX                                               |

---

## Requirements

- **OS:** Linux  
- **Container:** Docker installed and usable  
- **Network (for SRA mode):**
  - To pull the image from Docker Hub.
  - To access NCBI/ENA for `prefetch` / `fasterq-dump`.
- **Compute:**  
  - Designed for servers/HPC; large datasets can be computationally intensive and long-running.

---

## Installation

### Option 1 — Official installer (recommended)

```bash
curl -fsSL https://raw.githubusercontent.com/m13paiva/hulk/main/install_hulk.sh | bash
```

This script:

- Pulls the latest Docker image `m13paiva/hulk:latest`.
- Installs a wrapper at `/usr/local/bin/hulk`.
- Allows you to run `hulk` directly on the host.

**Uninstall wrapper:**
```bash
sudo rm /usr/local/bin/hulk
```

### Option 2 — Run via Docker (no install)

```bash
docker run --rm -v "$PWD":/data -w /data m13paiva/hulk:latest   -i SraRunInfo.csv -r transcripts.fa.gz -o results
```

---

## Quickstart

### SRA mode (RunInfo table)

```bash
# Minimal run using a transcriptome FASTA
hulk -i RunInfo.csv -r transcripts.fa.gz -o results

# Using a prebuilt kallisto index
hulk -i RunInfo.csv -r transcripts.idx -o results
```

In SRA mode, `-i/--input` must point to a **RunInfo-style table** (CSV/TSV/TXT) with at least the columns:

- `Run`
- `BioProject`
- `Model`

The table is typically obtained from NCBI SRA (“Send to → File → RunInfo”).  
The `Run` values are SRR accessions, `BioProject` is used to group runs into projects, and `Model` describes the sequencer/instrument.  
The `-r/--reference` option takes either a transcriptome FASTA (`.fa/.fasta[.gz]`) or a prebuilt `kallisto` index (`.idx`).

---

### FASTQ mode (local FASTQ directory)

```bash
# Local FASTQ directory with one subdirectory per sample
hulk -i fastq_runs/ -r transcripts.fa.gz -o results_fastq   --seq-tech "ILLUMINA NovaSeq 6000"
```

In FASTQ mode, `-i/--input` must point to a **directory**:

- Each **direct subdirectory** under this directory is treated as **one sample**.
- The **sample ID** is the name of the subdirectory.
- Each sample directory must contain:
  - **One** FASTQ/FASTQ.GZ file → treated as **single-end**, or
  - **Two** FASTQ/FASTQ.GZ files → treated as **paired-end**.
- Layout (single vs paired) is inferred automatically from the number of FASTQ files present.
- Deeper nested subdirectories are ignored; only FASTQ files directly inside the sample folder are considered.

Because there is no RunInfo metadata in FASTQ mode, you must provide `--seq-tech` with the sequencing technology (e.g. `ILLUMINA NovaSeq 6000`).  
HULK uses this to select sensible fragment-length defaults for `kallisto` single-end quantification.

---

## Input Table (SRA RunInfo)

To generate an SRA RunInfo table from NCBI:

1. Go to the **SRA Search**: https://www.ncbi.nlm.nih.gov/sra
2. Search by **BioProject** or paste SRA accessions (e.g. `PRJNA1141930`).
3. Select your runs of interest.
4. Click **“Send to:” → “File” → “RunInfo”**.
5. Download `RunInfo.csv`.

By default the **RunInfo** table already includes the columns requested by HULK:

- `Run`
- `BioProject`
- `Model`

**Example (`RunInfo.csv`):**
```csv
Run,BioProject,Model, ...
SRR30141434,PRJNA1141930,ILLUMINA HiSeq 2500, ...
SRR30141435,PRJNA1141930,ILLUMINA HiSeq 2500, ...
```

---

## FASTQ Mode

In FASTQ mode, HULK skips SRA download and runs directly on local FASTQ files.  
Instead of a RunInfo table, you provide a **directory** via `-i / --input`.

### Directory layout

Let’s say you have:

```text
fastq_runs/
├── SAMPLE_01/
│   ├── r1.fastq.gz
│   └── r2.fastq.gz
├── SAMPLE_02/
│   └── r.fastq.gz
├── SAMPLE_03/
│   ├── lane1/
│   │   └── something_else.txt
│   └── r1.fastq.gz
└── misc/
    └── not_a_sample.txt
```

- `fastq_runs/` is passed to HULK with `-i fastq_runs/`.
- HULK will look only at **direct subdirectories** of `fastq_runs/`:
  - `SAMPLE_01/`
  - `SAMPLE_02/`
  - `SAMPLE_03/`
  - `misc/`
- For each direct subdirectory:
  - If it contains **exactly one** FASTQ/FASTQ.GZ file (e.g. `r.fastq.gz`):
    - The sample is treated as **single-end**.
  - If it contains **exactly two** FASTQ/FASTQ.GZ files (e.g. `r1.fastq.gz`, `r2.fastq.gz`):
    - The sample is treated as **paired-end**.
  - If it contains **zero** FASTQ/FASTQ.GZ files, or **more than two**, that directory is ignored.
- The **sample ID** is taken from the **name of the subdirectory** (e.g. `SAMPLE_01`, `SAMPLE_02`, `SAMPLE_03`).
- Any deeper nested subdirectories (like `SAMPLE_03/lane1/`) are not interpreted as separate samples. Only the FASTQ files directly inside the sample directory are considered.

Both uncompressed (`.fastq`) and compressed (`.fastq.gz`) files are supported.

### Sequencing technology (`--seq-tech`)

In SRA mode, HULK can infer the sequencing technology from the `Model` column in the RunInfo table.  
In FASTQ mode, there is no such metadata, so the user must provide it:

- `--seq-tech TEXT`  
  A string describing the sequencing technology / instrument model, for example:
  - `"ILLUMINA NovaSeq 6000"`
  - `"ILLUMINA HiSeq 2500"`
  - `"DNBSEQ-G400"`

HULK uses `--seq-tech` in FASTQ mode for Kallisto parametrization and for basic descriptive labeling.  
If your FASTQ directory mixes data from multiple sequencing technologies, it is recommended to run HULK separately for each technology, with the appropriate `--seq-tech` value each time.

---

## CLI Reference

HULK uses a Click-based CLI with a main command and persistent configuration subcommands.

> Subcommand options are stored in a small JSON file inside the container  
> at `/app/.hulk.json`.

### Main command — `hulk [OPTIONS]`

```text
Required I/O
  -i, --input PATH
      SRA mode: RunInfo table (.csv, .tsv, or .txt) with columns: Run, BioProject, Model.
      FASTQ mode: directory containing sample subdirectories with FASTQ/FASTQ.GZ files.  [required]

  -r, --reference PATH
      Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).  [required]

Outputs & performance
  -o, --output PATH
      Output directory.  [default: output]

  --min-threads INTEGER
      Minimum number of threads per SRR/sample.  [default: 4]

  -t, --max-threads INTEGER
      Maximum total threads.  [default: 10]

Behaviour flags
  --verbosity / --no-verbosity
      Enable/disable live progress bars and console messages.  [default: --verbosity]

  -y, --yes
      Assume 'yes' to prompts and run without asking.

  -f, --force, --overwrite
      Force re-run: overwrite totally/partially processed SRRs/samples.

  -n, --dry-run
      Validate inputs and configuration, print plan, and exit without running tools.

Quantification and post-processing
  -g, --gene-counts PATH
      Enable gene-level counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'.
      When set, transcript-level quantifications are aggregated to genes automatically.

  --deseq2 / --no-deseq2
      Enable DESeq2 normalization + VST expression matrices for downstream plots/network exports.
      When disabled, kallisto quantifications still run but no DESeq2/VST matrices or DESeq2 plots
      are produced.  [default: --deseq2]

  --seq-tech TEXT
      Sequencing technology/platform (e.g. 'ILLUMINA NovaSeq 6000').
      Required in FASTQ mode; used to choose kallisto fragment-length parameters.

Caching and FASTQ retention
  --no-cache
      Disable SRA caching (do not allocate a shared cache volume).

  -c, --cache INTEGER
      Maximum SRA cache size (GiB).  Default: auto (300 GiB or based on free space).

  --keep-fastq
      Keep trimmed FASTQ files instead of deleting them after processing.

Misc
  -V, --version
      Show version and exit.

  -h, --help
      Show this message and exit.
```

If no subcommand is specified, running `hulk` with the options above executes the full pipeline.

---

### Subcommand — `hulk trim [OPTIONS]` (fastp defaults)

Configure default `fastp` trimming parameters used by the main pipeline.

```text
Options
  -ws, --window-size INTEGER
      fastp sliding window size.

  -mq, --mean-quality INTEGER
      fastp mean quality threshold.

  --reset-defaults
      Reset trim options to built-in defaults.

  -h, --help
      Show help and exit.
```

If neither `--window-size` nor `--mean-quality` is provided and `--reset-defaults` is not set, the command prints a short usage hint.

---

### Subcommand — `hulk tximport [OPTIONS]` (aggregation/normalization)

Configure how `tximport` aggregates and normalizes quantifications in the R post-processing step.

```text
Options
  -m, --mode [raw_counts|length_scaled_tpm|scaled_tpm|dtu_scaled_tpm]
      tximport aggregation/normalization mode.

  --ignore-tx-version
      Strip transcript version suffixes before matching (toggle stored as a boolean).

  --reset-defaults
      Reset tximport options to built-in defaults.

  -h, --help
      Show help and exit.
```

If no option is provided (and `--reset-defaults` is not set), the command prints a short usage hint.

---

### Subcommand — `hulk align [OPTIONS]` (alignment/quantification backend)

Configure the alignment/quantification backend (currently `kallisto`) and its options.

```text
Options
  --method [kallisto]
      Alignment/quantification backend.  [default: kallisto]

  -b, --bootstrap INTEGER
      Number of bootstrap samples for kallisto quant.  [default: 100]

  --reset-defaults
      Reset alignment options to built-in defaults.

  -h, --help
      Show help and exit.
```

---

### Subcommand — `hulk plot [OPTIONS]` (plotting behaviour)

Configure which plots are requested in post-processing.  
Plotting only takes effect when DESeq2/VST is enabled in the main command (`--deseq2`).

```text
Options
  --global-pca BOOL
      Enable/disable global PCA plot across all samples (true/false).  [default: true]

  --global-heatmap BOOL
      Enable/disable global expression heatmap across all samples (true/false).  [default: true]

  --global-var-heatmap BOOL
      Enable/disable global variance heatmap (gene x BioProject) (true/false).  [default: true]

  --bp-pca BOOL
      Enable/disable per-BioProject PCA plots (one PCA per BioProject) (true/false).  [default: false]

  --bp-heatmap BOOL
      Enable/disable per-BioProject expression heatmaps (one heatmap per BioProject) (true/false).  [default: false]

  --reset-defaults
      Reset plot options to built-in defaults.

  -h, --help
      Show help and exit.
```

---

## Usage Examples

```bash
# 1) Basic SRA-based run with a transcriptome FASTA
hulk -i RunInfo.csv -r transcripts.fa.gz -o results

# 2) Using a prebuilt kallisto index
hulk -i RunInfo.csv -r transcripts.idx -o results

# 3) Force re-run of everything (including post-processing)
hulk -i RunInfo.csv -r transcripts.fa.gz -o results -f

# 4) Generate gene counts (requires tx2gene mapping)
hulk -i RunInfo.csv -r transcripts.fa.gz -g tx2gene.csv -o results

# 5) Dry run (validate configuration without running tools)
hulk -i RunInfo.csv -r transcripts.fa.gz -o results -n --verbosity

# 6) Configure fastp defaults and then run
hulk trim -ws 30 -mq 20
hulk -i RunInfo.csv -r transcripts.fa.gz -o results

# 7) Configure tximport defaults and then run
hulk tximport -m raw_counts --ignore-tx-version
hulk -i RunInfo.csv -r transcripts.fa.gz -g tx2gene.csv -o results

# 8) Configure kallisto bootstraps
hulk align -b 250
hulk -i RunInfo.csv -r transcripts.fa.gz -o results

# 9) Configure plotting behaviour
hulk plot --global-pca true --global-heatmap true --global-var-heatmap true           --bp-pca false --bp-heatmap false

# 10) Run via Docker directly (no wrapper install)
docker run --rm -v "$PWD":/data -w /data m13paiva/hulk:latest   -i RunInfo.csv -r transcripts.fa.gz -o results --verbosity

# 11) FASTQ mode (local FASTQ directory)
# fastq_runs/ has one subdirectory per sample, each with 1 (SE) or 2 (PE) FASTQ/FASTQ.GZ files.
hulk -i fastq_runs/ -r transcripts.fa.gz -o results_fastq   --seq-tech "ILLUMINA NovaSeq 6000" --verbosity
```

---

## Output Structure

The exact layout can vary slightly depending on the user's inputs, but a typical directory for multiple BioProjects looks like:

```text
outdir/
├── PRJNA1141930/
│   ├── SRR30141434/
│   │   ├── fastp.json             ← fastp QC and trimming report 
│   │   ├── run_info.json          ← kallisto run metadata for this sample
│   │   ├── abundance.tsv          ← kallisto transcript-level quantification for this sample
│   │   └── SRR30141434.log        ← Sample-level log file
│   ├── SRRxxxxxxxx/               ← Additional samples (same structure)
│   ├── multiqc_PRJNA1141930.html  ← Project-level MultiQC report
│   ├── plots/                     ← Per-BioProject plots and their underlying matrices
│   ├── deseq2/                    ← DESeq2-normalized counts and VST matrices
│   ├── tximport/                  ← tximport summaries and intermediate artifacts
│   └── PRJNA1141930.log           ← BioProject-level log file
├── PRJNAxxxxxxx/                  ← Additional BioProjects (same structure)
└── shared/                        ← Shared (across all Bioprojects and Samples) output
    ├── plots/                     ← Global/cross-project plots and their underlying matrices
    ├── deseq2/                    ← Global DESeq2 outputs
    ├── tximport/                  ← Global tables from tximport
    ├── log.txt                    ← Global log file
    └── multiqc_global.html        ← Global MultiQC report 
```

---

## Tools used inside HULK

HULK orchestrates several third-party tools. If you use HULK in your research, please also consider citing the underlying tools where appropriate:

- **NCBI SRA Toolkit**  
    https://github.com/ncbi/sra-tools <br>
  Used for fetching and converting SRA data (`prefetch`, `fasterq-dump`) in SRA mode.

- **fastp**  
  https://github.com/OpenGene/fastp<br>
  Used for read QC and trimming.

- **kallisto**  
    https://github.com/pachterlab/kallisto<br>
  Used for pseudoalignment and transcript quantification.

- **MultiQC**  
    https://github.com/MultiQC/MultiQC <br>
  Used for aggregating QC metrics into per-project and global reports.

- **R and Bioconductor**  
  - **tximport** — imports transcript-level quantifications and summarizes them at the gene 
  level.  <br> https://github.com/thelovelab/tximport
  - **DESeq2** — computes normalized gene counts and VST matrices, and supports generation of diagnostic plots. <br> https://github.com/thelovelab/DESeq2

---

## License

**MIT License** © 2025 Manuel Paiva de Almeida

---

## Citation

If you use HULK in your research, please cite the software as:

>Paiva de Almeida, M., & Barros, P. **HULK: High-volume bulk RNA-seq pipeline.** https://github.com/m13paiva/hulk

In addition to citing HULK, please also cite the underlying tools used in the pipeline (SRA Toolkit, fastp, kallisto, MultiQC, tximport, DESeq2) as appropriate.

---

## Links

- GitHub: https://github.com/m13paiva/hulk  
- Docker Hub: https://hub.docker.com/r/m13paiva/hulk
