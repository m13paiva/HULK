# **HULK** Changelog
---

## [1.3.0] - 2025-12-17
This release focuses on stability at scale and fixing critical logic errors that affected large datasets and specific SRA formats. It also introduces parallelized plotting and targeted gene analysis.

### Fixed (Critical)
- **SRA "False Paired-End" Detection**: Resolved a critical bug where `fasterq-dump` would occasionally produce a second, empty FASTQ file for single-end reads (SRR), tricking the pipeline into treating them as paired-end. The layout detection logic now explicitly checks file sizes to verify if the second read file is valid before classifying the sample.
- **Large Dataset Scalability (DESeq2)**: Switched DESeq2 size factor estimation to `type="poscounts"`. The previous default (geometric mean) caused pipeline failures on large, sparse datasets (>1000 samples) with the error "every gene contains at least one zero".
- **Visual Rendering Crashes**: Heatmaps for massive datasets are now rasterized (`use_raster=TRUE`) with dynamic font sizing (6pt). This prevents the R process from running out of memory or generating illegible PDFs when plotting thousands of samples.

### Fixed (Logic & CLI)
- **Stale Cache Persistence**: Fixed an issue where the R script would load an old `txi.rds` cache file and ignore newly processed samples. The pipeline now forces a re-import (`--force-txi`) during the compute phase.
- **Variance Threshold Ignored**: The `var_threshold` parameter was previously parsed but not applied. It is now correctly used to filter low-variance genes before plotting.
- **CLI Variable Scope**: Fixed a `NameError` in `cli.py` where plot configuration variables were referenced before definition.
- **Argument Passing**: Corrected logic to ensure multiple `--target-genes` files are correctly joined and passed to the R script.

### Added
- **Targeted Gene Analysis** (`--target-genes`): Users can now provide one or more text files containing gene IDs. HULK will generate dedicated expression and variance heatmaps for these specific gene subsets, alongside the global plots.
- **New Quality Control Visualizations**:
  - **Sample-to-Sample Correlation Heatmap** (`--sample-cor`): Pearson correlation heatmap for outlier detection.
  - **DESeq2 Dispersion Plot** (`--dispersion`): Diagnostic plot for Negative Binomial model fit.
- **Configurable "Top N" Genes** (`--top-n`): Customize the number of top variable genes used for PCA and global heatmaps (default: 500).
- **Parallelized Post-Processing**:
  - **Split-Phase Execution**: Separated matrix computation (sequential) from plotting (concurrent) to save memory while speeding up visualization.
  - **Concurrent BioProjects**: Multiple BioProjects are now post-processed simultaneously using a thread pool.
- **Enhanced Configuration**: `hulk plot` now saves defaults for all new plot types (`--sample-cor`, `--dispersion`) and settings (`--top-n`) to `.hulk.json`.

### Changed
- **`hulk report` Workflow**:
  - Removed redundant plotting flags from the `report` command to prevent "flag blindness".
  - The command now strictly follows the saved configuration from `hulk plot`, defaulting to `True` for standard plots if no config exists.
  - Added verbose logging to startup to explicitly list active plots.

---

## [1.2.0] 2025-12-04
This will be the last feature release focused on standard bulk RNA-seq workflows (apart from potential bug-fix releases). Future versions of **HULK** will be centered on gene co-expression network (GCN) analysis using **Seidr**.

### Added
- New **FASTQ mode**:
  - Allows running HULK directly on user-provided FASTQ files instead of SRR accessions.
  - Automatically detects FASTQ layout and integrates with the existing pipeline flow.

- Clear separation between **prefetch** and **processing** workers:
  - A dedicated prefetch loop continuously downloads SRR data into a shared cache, independent of the processing pool.
  - Processing workers only operate on already-prefetched data, reducing idle time caused by I/O waits.
  - Leads to better thread utilization, smoother load distribution across BioProjects, and noticeably improved runtimes on large or heavily parallel runs.

- New **R-based post-processing** pipeline (replacing the previous `pytximport`-based flow):
  - Uses **tximport** to import quantification results.
  - Runs **DESeq2** to generate normalized count matrices and variance-stabilizing transform (VST) expression matrices (no metadata usage on DESeq2 supported). Controlled via a `--deseq2 / --no-deseq2` flag.
  - Produces a set of plots (variance heatmaps, expression heatmaps, and PCA) and saves both figures and underlying matrices for downstream reuse.

- Clearer output layout for post-processing results directly under `shared/`:
  - `shared/plots/` for all post-processing plots and their associated matrix files (variance, expression, etc.).
  - `shared/deseq2/` for DESeq2-normalized counts and VST matrices.
  - `shared/tximport/` for `tximport` summaries and intermediate artifacts.

- More BioProject- and sample-aware outputs:
  - Per-BioProject plots and summaries in post-processing.
  - Per-BioProject log files.
  - Per-sample log files.
  - The reorganized output layout:
    - `bioproject_id/plots/`
    - `bioproject_id/deseq2/`
    - `bioproject_id/tximport/`

### Changed
- Aggregation semantics and expression modes:
  - The old `aggregate` option has been removed.
  - The pipeline is now standardized on **gene-level counts**:
    - When starting from transcript-level quantifications, aggregation to genes is always performed automatically in gene-count mode.
    - TPM aggregation functionality has been dropped; output for downstream analysis is now consistently based on gene counts.

### Removed
- The `aggregate` option and TPM aggregation mode:
  - No more separate TPM aggregation path.
  - The pipeline now consistently produces gene-level counts, automatically aggregating transcripts to genes when needed.

- The old `pytximport`-based post-processing:
  - All responsibilities formerly handled by `pytximport` have been migrated to the R + `tximport` + DESeq2 workflow described above.

### Fixed
- Improved robustness of repeated post-processing runs:
  - Re-running post-processing on an already processed output directory no longer corrupts or duplicates results.
  - Interrupted post-processing stages are handled more gracefully, reusing completed outputs whenever possible.

- Logging and progress bars:
  - The state revealed by progress bars and log messages is more tightly synchronized with what actually happens per BioProject and per sample, reducing confusion during long runs.


---

## [1.1.1] - 2025-10-29
### Fixed
- **Fixed bug in** `_detect_fastq_layout_` **from** `_utils.py_`:
Previously, since the option `--split-files` is always used with `fasterq-dump`, even for single-end SRRs, a `_1` suffix was appended to the output file name. This caused `_detect_fastq_layout_` to treat the SRR as incomplete and abort further processing.
In this patched version, `_detect_fastq_layout_` now recognizes a single FASTQ file with the `_1` suffix as a valid single-end SRR and automatically corrects the filename by removing the suffix.


---

## [1.1.0] - 2025-10-26
### Added
- Switched CLI from **argparse** to **Click** with a custom help renderer.
- **Subcommands** to configure defaults saved locally:
  - `hulk trim` with options:
    - `-ws, --window-size <int>`
    - `-mq, --mean-quality <int>`
  - `hulk tximport` with options:
    - `-m, --mode {raw_counts,length_scaled_tpm,scaled_tpm,dtu_scaled_tpm}`
    - `--ignore-tx-version`
  - Subcommand settings are persisted in a JSON file (default `./.hulk.json`, overridable via `HULK_CONFIG`) and automatically applied when running `hulk ...`.

---

## [1.0.0] - 2025-10-21
### Added
- Initial public release of **HULK**.

---

[1.3.0]: https://github.com/m13paiva/hulk/releases/tag/v1.3.0
[1.2.0]: https://github.com/m13paiva/hulk/releases/tag/v1.2.0
[1.1.1]: https://github.com/m13paiva/hulk/releases/tag/v1.1.1
[1.1.0]: https://github.com/m13paiva/hulk/releases/tag/v1.1.0
[1.0.0]: https://github.com/m13paiva/hulk/releases/tag/v1.0.0
