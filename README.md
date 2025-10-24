# HULK — High-Volume Bulk RNA-seq Pipeline

HULK is a fully automated, containerized pipeline for processing **bulk RNA-seq data directly from SRA accessions**.  
It handles everything from **data fetching** to **quantification and report generation** — with just one command.

---

## Features

- **Automated workflow**:
  - Fetches runs from SRA (`prefetch`)
  - Converts to FASTQ (`fasterq-dump`)
  - Performs QC and trimming with **fastp**
  - Quantifies transcript abundances with **kallisto**
  - Generates MultiQC reports automatically (project-level and global)
- **Supports paired-end and single-end reads**
- **Supports major sequencing platforms**, including:
  - **Illumina NovaSeq series** (6000, X, X Plus)
  - **Illumina HiSeq family** (X Ten, 4000, 3000, 2500, 2000, 1500, 1000)
  - **Illumina NextSeq series** (1000, 500, 550)
  - **BGI/MGI platforms** (DNBSEQ-G400, DNBSEQ-T7, BGISEQ-500, MGISEQ-2000RS)
  - **Legacy Illumina Genome Analyzer II / IIX**
- **Aggregates TPM and gene counts** across BioProjects
- **Containerized** — no environment setup needed
- **CPU-aware concurrency** with automatic workload distribution
- **Clean, resumable, and safe overwrites**

---

## 🧰 Installation

### Option 1 — Install via the official installer (recommended)
```bash
curl -fsSL https://raw.githubusercontent.com/m13paiva/hulk/main/install_hulk.sh | bash
```

This will:
- Pull the latest HULK Docker image (`m13paiva/hulk:latest`)
- Create a global wrapper (`/usr/local/bin/hulk`)
- Allow you to run `hulk` directly like a normal command

✅ Once installed, just run:
```bash
hulk -i SraRunInfo.csv -r transcripts.fa.gz -o results
```

To uninstall:
```bash
sudo rm /usr/local/bin/hulk
```

---

### Option 2 — Run directly from Docker (no installation)
```bash
docker run --rm -v "$PWD":/data -w /data m13paiva/hulk:latest   -i SraRunInfo.csv -r transcripts.fa.gz -o results
```

---

## Preparing Your Input Table

You don’t need to manually create any metadata — NCBI can generate it for you automatically.

1. Go to the [NCBI SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/).  
2. Search for your **BioProject** or list of SRA accessions (e.g., `PRJNA1141930`).
3. Select the runs you want to include.
4. Click **“Send to:” → “File” → “RunInfo”**.
5. Download the file — it will be named something like `SraSraRunInfo.csv`.  

That file is **ready to use as input** for HULK — it already contains the required columns:
`Run`, `BioProject`, and `Model`.

Example (`SraRunInfo.csv`):
```csv
Run,BioProject,Model
SRR30141434,PRJNA1141930,ILLUMINA HiSeq 2500
SRR30141435,PRJNA1141930,ILLUMINA HiSeq 2500
```

---

## Running the Pipeline

HULK accepts:
- Transcriptome FASTA (`.fa`, `.fasta`, `.fa.gz`)
- or a pre-built kallisto index (`.idx`)

Example command:
```bash
hulk -i SraRunInfo.csv -r transcripts.fa.gz -o results
```

Optional flags:
```bash
  -f, --force          Re-run and overwrite previous outputs
  -a, --aggregate      Generate global TPM and gene count tables (and global gene counts if enabled)
  -g, --gene-counts    Provide a tx2gene.csv file for gene-level output
  -n, --dry-run        Validate input and show plan without running
  -v, --verbose        Show progress bars and tool logs
  -V, --version        Show current version
```

---

## 📂 Output Structure

After successful execution, your output folder will look like this:

```
results/
├── PRJNA1141930/
│   ├── SRR30141434/
│   │   ├── fastp.json
│   │   ├── run_info.json
│   │   ├── abundance.tsv
│   │   └── kallisto.log
│   ├── _mqc_inputs/
│   ├── multiqc_PRJNA1141930.html      ← Project-level MultiQC report
│   └── read_metrics.tsv               ← Read count summary
├── aggregated/
│   ├── tpm_matrix.tsv                 ← Global TPM table (if -a)
│   └── gene_counts.tsv                ← Global gene counts (if -a and -g)
└── multiqc_global.html                ← Global MultiQC report (all projects)
```

---

## License

MIT License © 2025 **Manuel Paiva de Almeida**

---

## Citation

If you use HULK in your research, please cite:
>Almeida, M., & Barros, P. (2025). *HULK: A reproducible, scalable pipeline for high-volume bulk RNA-seq data processing (included in the Rice2B project)*. Instituto de Tecnologia Química e Biológica António Xavier, Universidade Nova de Lisboa (ITQB NOVA), GPlantS Lab.

---

## 🔗 Links

- GitHub: [https://github.com/m13paiva/hulk](https://github.com/m13paiva/hulk)  
- Docker Hub: [https://hub.docker.com/r/m13paiva/hulk](https://hub.docker.com/r/m13paiva/hulk)

---


