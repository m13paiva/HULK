#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)  # for heatmaps
})

`%||%` <- function(x, y) if (is.null(x)) y else x

## ===========================
## CLI OPTIONS
## ===========================

option_list <- list(
  make_option(c("-s", "--search-dir"), dest = "search_dir",
              type = "character", default = "output",
              help = "Root folder that contains sample folders with abundance.tsv [default: %default]"),
  make_option(c("-x", "--exclude-dir"), dest = "exclude_dir",
              type = "character", action = "append", default = c("_mqc_inputs"),
              help = "Folder name to exclude (can be given multiple times). [default: %default]"),
  make_option(c("-t", "--tx2gene"), dest = "tx2gene",
              type = "character", default = "tx2gene.clean.tsv",
              help = "Path to tx2gene file (transcript_id, gene_id) [default: %default]"),
  make_option(c("-o", "--out-dir"), dest = "out_dir",
              type = "character", default = ".",
              help = "Output root directory (HULK 'shared' or BioProject dir) [default: %default]"),
  make_option(c("-p", "--prefix"), dest = "prefix",
              type = "character", default = "",
              help = "Optional output prefix (currently unused for filenames) [default: '%default']"),
  make_option(c("--bioproject"), dest = "bioproject",
              type = "character", default = NULL,
              help = "Optional BioProject ID to restrict analysis to (per-BP mode). [default: all]"),
  make_option(c("--ignore-tx-version"), dest = "ignore_tx_version",
              action = "store_true", default = FALSE,
              help = "Ignore transcript version suffixes in tximport [default: %default]"),
  make_option(c("--counts-from-abundance"), dest = "counts_from_abundance",
              type = "character", default = "no",
              help = "tximport countsFromAbundance: no | scaledTPM | lengthScaledTPM | dtuScaledTPM [default: %default]"),
  make_option(c("--force-txi"), dest = "force_txi",
              action = "store_true", default = FALSE,
              help = "Force re-running tximport even if cached RDS exists [default: %default]"),
  make_option(c("--use-matrix"), dest = "use_matrix",
              type = "character", default = "vst",
              help = "Matrix for Seidr-style export: vst | normalized [default: %default]"),
  make_option(c("--no-drop-nonvarying"), dest = "drop_nonvarying",
              action = "store_false", default = TRUE,
              help = "Do NOT drop non-varying genes in Seidr-style exports [default: drop]"),
  make_option(c("--pca"), dest = "pca",
              action = "store_true", default = FALSE,
              help = "Generate PCA plot (VST) [default: %default]"),
  make_option(c("--heatmap"), dest = "heatmap",
              action = "store_true", default = FALSE,
              help = "Generate global / targeted expression heatmap (grouped by BioProject or sample) [default: %default]"),
  make_option(c("--var-heatmap"), dest = "var_heatmap",
              action = "store_true", default = FALSE,
              help = "Generate gene x BioProject variance heatmap (VST); disabled in --bioproject mode and in FASTQ mode [default: %default]"),
  make_option(c("--plots-only"), dest = "plots_only",
              action = "store_true", default = FALSE,
              help = "Only build PCA/heatmaps from existing VST file; skip tximport/DESeq2/Seidr [default: %default]"),
  make_option(c("--tximport-only"), dest = "tximport_only",
              action = "store_true", default = FALSE,
              help = "Only run tximport and write gene x sample counts table; skip DESeq2, VST, plotting, and Seidr [default: %default]"),
  make_option(c("--target-genes"), dest = "target_genes",
              type = "character", default = NULL,
              help = "Path to file with target gene IDs (one per line) for targeted plots [default: %default]"),
  make_option(c("--mode"), dest = "mode",
              type = "character", default = "SRR",
              help = "Input mode: SRR or FASTQ (affects plot labels and variance heatmap) [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

## Map options
SEARCH_DIR        <- opt$search_dir
EXCLUDE_DIRS      <- opt$exclude_dir %||% character(0)
TX2GENE_PATH      <- opt$tx2gene
OUT_DIR           <- opt$out_dir
PREFIX            <- opt$prefix  # kept for potential future use, but NOT used in filenames
BIOPROJECT_ID     <- opt$bioproject
IGNORE_TX_VERSION <- isTRUE(opt$ignore_tx_version)
COUNTS_FROM_AB    <- opt$counts_from_abundance
FORCE_TXI         <- isTRUE(opt$force_txi)
USE_MATRIX        <- tolower(opt$use_matrix)
DROP_NONVARYING   <- isTRUE(opt$drop_nonvarying)
DO_PCA            <- isTRUE(opt$pca)
DO_HEATMAP        <- isTRUE(opt$heatmap)
DO_VAR_HM         <- isTRUE(opt$var_heatmap)
PLOTS_ONLY        <- isTRUE(opt$plots_only)
TXIMPORT_ONLY     <- isTRUE(opt$tximport_only)
TARGET_GENES_PATH <- opt$target_genes
MODE              <- toupper(opt$mode %||% "SRR")

# Per-BioProject mode flag (used to treat plots like FASTQ: sample-centred)
PER_BP_MODE <- !is.null(BIOPROJECT_ID)

## Basic sanity check
if (PLOTS_ONLY && TXIMPORT_ONLY) {
  stop("Cannot use --plots-only together with --tximport-only.")
}

## Disable variance heatmap in FASTQ mode regardless of CLI flag
if (MODE == "FASTQ" && DO_VAR_HM) {
  message("[VarHeatmap] Disabled in FASTQ mode (--mode=FASTQ).")
  DO_VAR_HM <- FALSE
}

## ===========================
## OUTPUT PATHS & LAYOUT
## ===========================

# OUT_DIR is now the root (global: <OUTPUT>/shared; BP: <OUTPUT>/<BioProject>)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

plots_dir    <- file.path(OUT_DIR, "plots")
deseq2_dir   <- file.path(OUT_DIR, "deseq2")
tximport_dir <- file.path(OUT_DIR, "tximport")

dir.create(plots_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(deseq2_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(tximport_dir, showWarnings = FALSE, recursive = TRUE)

# tximport-related
txi_rds        <- file.path(tximport_dir, "txi.rds")
txi_counts_tsv <- file.path(tximport_dir, "tximport_counts.tsv")

# DESeq2 / VST / normalized
dds_rds  <- file.path(deseq2_dir, "dds.rds")
vst_tsv  <- file.path(deseq2_dir, "vst.tsv")
norm_tsv <- file.path(deseq2_dir, "normalized_counts.tsv")
sf_tsv   <- file.path(deseq2_dir, "size_factors.tsv")

# Plots (and their sidecar matrices)
pca_pdf         <- file.path(plots_dir, "PCA.pdf")
heatmap_pdf     <- file.path(plots_dir, "expression_heatmap.pdf")
var_heatmap_pdf <- file.path(plots_dir, "variance_heatmap.pdf")

## ===========================
## HELPERS
## ===========================

msg <- function(...) cat(sprintf(...), "\n", sep = "")

find_abundance_files <- function(search_dir, exclude_dirs) {
  all <- list.files(search_dir, pattern = "^abundance\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(exclude_dirs)) {
    keep <- !vapply(all, function(p) {
      any(basename(dirname(p)) %in% exclude_dirs |
            basename(dirname(dirname(p))) %in% exclude_dirs)
    }, logical(1))
    all <- all[keep]
  }
  sort(all)
}

read_tx2gene <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, col_types = cols(.default = "c")))
  if (ncol(df) < 2) stop("tx2gene must have at least 2 columns")
  colnames(df)[1:2] <- c("transcript_id", "gene_id")
  df[, 1:2, drop = FALSE]
}

is_nonvarying <- function(x) {
  x <- as.numeric(x); x[!is.finite(x)] <- NA_real_
  length(unique(na.omit(x))) < 2
}

write_seidr <- function(M, out_dir, drop_nonvarying = FALSE) {
  # out_dir is now typically dese2_dir (no "Seidr" folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  M <- as.data.frame(M)        # samples x genes
  M[] <- lapply(M, as.numeric)
  genes <- make.unique(colnames(M))

  expr_path  <- file.path(out_dir, "expression.tsv")
  genes_path <- file.path(out_dir, "genes.txt")
  readr::write_tsv(M, expr_path, col_names = FALSE)
  readr::write_lines(genes, genes_path)
  msg("[ExprExport] Wrote %s (shape %d x %d) and %s (n=%d)",
      expr_path, nrow(M), ncol(M), genes_path, length(genes))

  if (drop_nonvarying) {
    keep <- vapply(M, function(col) !is_nonvarying(col), logical(1))
    Mf <- M[, keep, drop = FALSE]
    expr_f  <- file.path(out_dir, "expression.nonzvar.tsv")
    genes_f <- file.path(out_dir, "genes.nonzvar.txt")
    readr::write_tsv(Mf, expr_f, col_names = FALSE)
    readr::write_lines(make.unique(colnames(Mf)), genes_f)
    msg("[ExprExport] Dropped %d non-varying genes; kept %d", ncol(M) - ncol(Mf), ncol(Mf))
    msg("[ExprExport] Wrote %s and %s", expr_f, genes_f)
  }
}

# Read target genes from file (one gene ID per line)
read_target_genes <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    stop(sprintf("Target gene list not found: %s", path))
  }
  genes <- readLines(path, warn = FALSE)
  genes <- unique(trimws(genes))
  genes <- genes[genes != ""]
  genes
}

## ===========================
## PLOTTING
## ===========================

make_plots <- function(vst_mat,
                       sample_info,
                       DO_PCA,
                       DO_EXPR_HM,
                       DO_VAR_HM,
                       pca_pdf,
                       expr_heatmap_pdf,
                       var_heatmap_pdf,
                       target_genes = NULL) {

  if (!DO_PCA && !DO_EXPR_HM && !DO_VAR_HM) {
    return(invisible(NULL))
  }

  # Align sample_info to cols of vst_mat
  sample_info <- sample_info %>%
    dplyr::filter(sample %in% colnames(vst_mat)) %>%
    dplyr::arrange(match(sample, colnames(vst_mat)))

  if (!nrow(sample_info)) {
    stop("No overlap between vst matrix samples and sample_info samples.")
  }

  # -------- Decide when we need top-variable genes --------
  need_top <- DO_PCA ||
    ((DO_EXPR_HM || DO_VAR_HM) &&
       (is.null(target_genes) || length(target_genes) == 0L))

  vst_top <- NULL
  if (need_top) {
    gene_var <- apply(vst_mat, 1, var, na.rm = TRUE)
    ord_genes <- order(gene_var, decreasing = TRUE)
    if (!length(ord_genes)) stop("No genes with finite variance for PCA/heatmap.")
    top_n <- min(500L, length(ord_genes))
    keep_genes <- ord_genes[seq_len(top_n)]
    vst_top <- vst_mat[keep_genes, , drop = FALSE]  # genes x samples
  }

  # Helper: choose gene set for variance heatmap
  get_gene_set <- function() {
    if (!is.null(target_genes) && length(target_genes) > 0L) {
      g <- intersect(target_genes, rownames(vst_mat))
      if (!length(g)) {
        msg("[Plots] WARNING: none of the target genes were found in the VST matrix.")
        return(character(0))
      }
      return(g)
    } else {
      return(rownames(vst_top))
    }
  }

  ## ---------------- PCA ----------------
  if (DO_PCA) {
    msg("[PCA] Computing PCA (VST, top variable genes)…")
    mat_pca <- t(vst_top)  # samples x genes
    pca <- prcomp(mat_pca, center = TRUE, scale. = FALSE)

    pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    pca_df$sample <- rownames(pca_df)
    pca_df <- pca_df %>%
      dplyr::left_join(sample_info, by = "sample")

    if (MODE == "FASTQ" || PER_BP_MODE) {
      # FASTQ or per-BioProject SRR: colour by sample, legend "Samples"
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = sample)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_bw() +
        labs(
          title  = "PCA (VST-normalized expression)",
          colour = "Samples"
        )
    } else {
      # Global SRR mode: colour by BioProject
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = bioproject)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_bw() +
        labs(
          title  = "PCA (VST-normalized expression)",
          colour = "BioProject"
        )
    }

    ggsave(pca_pdf, plot = p, width = 6, height = 5)
    msg("[PCA] Saved PCA plot to %s", pca_pdf)
  }

  ## ---------------- Expression heatmap ----------------
  if (DO_EXPR_HM) {
    # Decide grouping for columns: BioProject (global SRR) vs Samples (FASTQ or per-BP)
    if (MODE == "FASTQ" || PER_BP_MODE) {
      # FASTQ and per-BioProject SRR: order/group by sample
      ord_samples     <- order(sample_info$sample)
      sample_info_ord <- sample_info[ord_samples, , drop = FALSE]

      group_vec      <- sample_info_ord$sample
      group_title    <- "Samples"
      legend_suffix  <- ".sample_legend.tsv"
    } else {
      # Global SRR: order by BioProject then sample, group by BioProject
      ord_samples     <- order(sample_info$bioproject, sample_info$sample)
      sample_info_ord <- sample_info[ord_samples, , drop = FALSE]

      group_vec      <- sample_info_ord$bioproject
      group_title    <- "BioProject group"
      legend_suffix  <- ".bp_legend.tsv"
    }

    grp_levels <- unique(group_vec)
    group_ids  <- seq_along(grp_levels)
    group_factor <- factor(
      match(group_vec, grp_levels),
      levels = group_ids,
      labels = as.character(group_ids)
    )

    # Choose matrix: targeted vs global
    heat_mat <- NULL
    if (!is.null(target_genes) && length(target_genes) > 0L) {
      msg("[Heatmap] Building TARGETED expression heatmap for supplied gene list…")
      idx <- match(target_genes, rownames(vst_mat))
      idx <- idx[!is.na(idx)]
      if (!length(idx)) {
        msg("[Heatmap] WARNING: None of the target genes were found in the VST matrix; skipping expression heatmap.")
      } else {
        heat_mat <- vst_mat[idx, ord_samples, drop = FALSE]
        msg("[Heatmap] Using %d of %d requested target genes.", nrow(heat_mat), length(target_genes))
      }
    } else {
      msg("[Heatmap] Building GLOBAL expression heatmap (top variable genes)…")
      heat_mat <- vst_top[, ord_samples, drop = FALSE]
    }

    if (!is.null(heat_mat)) {
      # Legend table written to disk (same folder as heatmap_pdf, i.e. plots/)
      if (MODE == "FASTQ" || PER_BP_MODE) {
        legend_df <- data.frame(
          group_id = group_ids,
          sample   = grp_levels,
          stringsAsFactors = FALSE
        )
      } else {
        legend_df <- data.frame(
          group_id   = group_ids,
          bioproject = grp_levels,
          stringsAsFactors = FALSE
        )
      }

      legend_path <- sub("\\.pdf$", legend_suffix, heatmap_pdf)
      readr::write_tsv(legend_df, legend_path)
      msg("[Heatmap] Saved %s numeric legend to %s", group_title, legend_path)

      # Draw heatmap: column_split gives numeric group labels 1..N
      msg("[Heatmap] Saving expression heatmap to %s", heatmap_pdf)
      grDevices::pdf(heatmap_pdf, width = 10, height = 6)
      ht <- ComplexHeatmap::Heatmap(
        heat_mat,
        name = "vst",
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_split = group_factor,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          group = group_factor,
          show_annotation_name = FALSE,  # hide the extra "group" label near the bar
          annotation_legend_param = list(
            group = list(
              title  = group_title,
              at     = as.character(group_ids),
              labels = paste(group_ids, grp_levels, sep = " = ")
            )
          )
        ),
        column_title_gp = grid::gpar(fontsize = 6)
      )
      ComplexHeatmap::draw(ht, merge_legend = TRUE)
      grDevices::dev.off()
      msg("[Heatmap] Saved expression heatmap to %s", heatmap_pdf)
    }
  }

  ## ---------------- Variance heatmap (gene x BioProject) ----------------
  if (DO_VAR_HM) {
    msg("[VarHeatmap] Building variance heatmap (gene x BioProject)…")
    genes_use <- get_gene_set()
    if (!length(genes_use)) {
      msg("[VarHeatmap] No genes available for variance heatmap; skipping.")
    } else {
      vst_sub <- vst_mat[genes_use, , drop = FALSE]

      df <- t(vst_sub) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("sample") %>%
        dplyr::left_join(sample_info, by = "sample") %>%
        tidyr::pivot_longer(
          cols = all_of(genes_use),
          names_to = "gene",
          values_to = "vst"
        )

      var_tbl <- df %>%
        dplyr::group_by(gene, bioproject) %>%
        dplyr::summarise(var_vst = stats::var(vst, na.rm = TRUE), .groups = "drop")

      var_mat <- var_tbl %>%
        tidyr::pivot_wider(names_from = bioproject, values_from = var_vst) %>%
        tibble::column_to_rownames("gene") %>%
        as.matrix()

      # Save matrix next to variance_heatmap.pdf (plots/)
      var_tsv <- sub("\\.pdf$", ".variance_matrix.tsv", var_heatmap_pdf)
      readr::write_tsv(
        as.data.frame(var_mat) %>% tibble::rownames_to_column("gene"),
        var_tsv
      )
      msg("[VarHeatmap] Saved variance matrix to %s", var_tsv)

      # Plot
      grDevices::pdf(var_heatmap_pdf, width = 8, height = 6)
      ht_var <- ComplexHeatmap::Heatmap(
        var_mat,
        name = "var(VST)",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE
      )
      ComplexHeatmap::draw(ht_var)
      grDevices::dev.off()
      msg("[VarHeatmap] Saved variance heatmap to %s", var_heatmap_pdf)
    }
  }

  invisible(NULL)
}

## ===========================
## COLLECT abundance.tsv FILES
## ===========================

files <- find_abundance_files(SEARCH_DIR, EXCLUDE_DIRS)
if (!length(files)) stop("No 'abundance.tsv' found under SEARCH_DIR (after exclusions).")
msg("Found %d abundance.tsv files under %s", length(files), SEARCH_DIR)

# Sample name = parent dir; BioProject = grandparent dir
sample_names_all <- basename(dirname(files))
bioprojects_all  <- basename(dirname(dirname(files)))

sample_info <- tibble(
  sample = sample_names_all,
  bioproject = bioprojects_all
)

names(files) <- sample_info$sample

# Optional per-BioProject filter (BP mode)
if (!is.null(BIOPROJECT_ID)) {
  keep <- sample_info$bioproject == BIOPROJECT_ID
  if (!any(keep)) {
    stop(sprintf("No samples for BioProject '%s' found under SEARCH_DIR.", BIOPROJECT_ID))
  }
  sample_info <- sample_info[keep, , drop = FALSE]
  files <- files[sample_info$sample]
  msg("[Mode] Restricting analysis to BioProject '%s' (%d samples).",
      BIOPROJECT_ID, nrow(sample_info))

  if (DO_VAR_HM) {
    msg("[VarHeatmap] Disabled in per-BioProject mode (--bioproject=%s).", BIOPROJECT_ID)
    DO_VAR_HM <- FALSE
  }
}

# Update sample_names to current set (possibly filtered)
sample_names <- sample_info$sample

## ===========================
## TARGET GENE LIST (optional)
## ===========================

TARGET_GENES <- NULL
if (!is.null(TARGET_GENES_PATH)) {
  TARGET_GENES <- read_target_genes(TARGET_GENES_PATH)
  msg("Loaded %d target genes from %s", length(TARGET_GENES), TARGET_GENES_PATH)
}

## ===========================
## PLOTS-ONLY MODE
## ===========================

if (PLOTS_ONLY) {
  msg("[Mode] plots-only: using existing VST file to build PCA/heatmaps.")

  if (!DO_PCA && !DO_HEATMAP && !DO_VAR_HM) {
    msg("All plot options are disabled; nothing to do in plots-only mode.")
    quit(status = 0)
  }

  if (!file.exists(vst_tsv)) {
    stop(sprintf("plots-only mode requires existing VST file: %s (run once without --plots-only)", vst_tsv))
  }

  vst_df <- readr::read_tsv(vst_tsv, col_types = cols())

  # New format: first column 'sample'
  if ("sample" %in% colnames(vst_df)) {
    sample_col <- vst_df$sample
    mat <- as.matrix(vst_df[, -1, drop = FALSE])
    rownames(mat) <- sample_col
    vst_mat <- t(mat)  # genes x samples
  } else {
    mat <- as.matrix(vst_df)
    rownames(mat) <- sample_names
    vst_mat <- t(mat)  # genes x samples
  }

  # If per-BP mode in plots-only, restrict to that BP's samples
  if (!is.null(BIOPROJECT_ID)) {
    common <- intersect(colnames(vst_mat), sample_info$sample)
    if (!length(common)) {
      stop(sprintf("No overlap between VST matrix samples and BioProject '%s' samples.", BIOPROJECT_ID))
    }
    vst_mat <- vst_mat[, common, drop = FALSE]
    sample_info <- sample_info[sample_info$sample %in% common, , drop = FALSE]
  }

  make_plots(
    vst_mat,
    sample_info,
    DO_PCA,
    DO_HEATMAP,
    DO_VAR_HM,
    pca_pdf,
    heatmap_pdf,
    var_heatmap_pdf,
    target_genes = TARGET_GENES
  )
  msg("Done (plots-only).")
  quit(status = 0)
}

## ===========================
## tximport + (optional) DESeq2 + plots + Seidr
## ===========================

tx2gene <- read_tx2gene(TX2GENE_PATH)

if (file.exists(txi_rds) && !FORCE_TXI) {
  msg("[tximport] Loading cached %s", txi_rds)
  txi <- readRDS(txi_rds)

  # If cached txi was built on more samples, subset to current sample_names
  if (!all(colnames(txi$counts) == sample_names)) {
    common <- intersect(colnames(txi$counts), sample_names)
    if (!length(common)) {
      stop("Cached txi object shares no samples with current dataset.")
    }
    txi$counts    <- txi$counts[, common, drop = FALSE]
    txi$abundance <- txi$abundance[, common, drop = FALSE]
    txi$length    <- txi$length[, common, drop = FALSE]
  }
} else {
  msg("[tximport] Running tximport (kallisto)…")
  txi <- tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = COUNTS_FROM_AB,
    ignoreTxVersion = IGNORE_TX_VERSION
  )
  saveRDS(txi, txi_rds)
  msg("[tximport] Saved %s", txi_rds)
}

## ===========================
## tximport-only mode
## ===========================

if (TXIMPORT_ONLY) {
  msg("[Mode] tximport-only: writing gene x sample counts and exiting.")
  counts_mat <- txi$counts
  counts_df <- as.data.frame(counts_mat)
  counts_df <- tibble::rownames_to_column(counts_df, var = "gene_id")
  readr::write_tsv(counts_df, txi_counts_tsv)
  msg("[tximport-only] Wrote tximport counts table to %s", txi_counts_tsv)
  msg("Done (tximport-only).")
  quit(status = 0)
}

## ===========================
## DESeq2 + VST + plots + Seidr
## ===========================

# Align sample_info to txi columns
sample_info <- sample_info %>%
  dplyr::filter(sample %in% colnames(txi$counts)) %>%
  dplyr::arrange(match(sample, colnames(txi$counts)))

coldata <- data.frame(
  sample = colnames(txi$counts),
  row.names = colnames(txi$counts),
  stringsAsFactors = FALSE
)
coldata$bioproject <- sample_info$bioproject[match(rownames(coldata), sample_info$sample)]

dds <- suppressMessages(DESeqDataSetFromTximport(
  txi = txi,
  colData = coldata,
  design = ~ 1
))

# Filter genes with no counts
keep <- rowSums(counts(dds)) > 0
dds  <- dds[keep, ]

# Size factors
dds <- estimateSizeFactors(dds, type = "poscounts")
sf <- tryCatch(sizeFactors(dds), error = function(e) NULL)
if (is.null(sf)) {
  sf <- rep(1, ncol(dds))
  names(sf) <- colnames(dds)
}

# Normalized counts (genes x samples)
norm_counts <- counts(dds, normalized = TRUE)

# VST (genes x samples)
vst_obj <- suppressMessages(vst(dds, blind = TRUE))
vst_mat <- assay(vst_obj)

# Save matrices with samples in rows; include sample column for reloadability
vst_out <- as.data.frame(t(vst_mat))
vst_out <- tibble::rownames_to_column(vst_out, var = "sample")
readr::write_tsv(vst_out, file = vst_tsv, col_names = TRUE)

norm_out <- as.data.frame(t(norm_counts))
norm_out <- tibble::rownames_to_column(norm_out, var = "sample")
readr::write_tsv(norm_out, file = norm_tsv, col_names = TRUE)

readr::write_tsv(
  data.frame(sample = names(sf), size_factor = as.numeric(sf)),
  file = sf_tsv
)

saveRDS(dds, dds_rds)
msg("Saved:")
msg("  VST matrix:        %s", vst_tsv)
msg("  Normalized counts: %s", norm_tsv)
msg("  Size factors:      %s", sf_tsv)
msg("  DDS RDS:           %s", dds_rds)

# Plots (PCA + expression / variance heatmaps)
make_plots(
  vst_mat,
  sample_info,
  DO_PCA,
  DO_HEATMAP,
  DO_VAR_HM,
  pca_pdf,
  heatmap_pdf,
  var_heatmap_pdf,
  target_genes = TARGET_GENES
)

# Seidr-ready exports (now placed directly under dese2_dir, no "Seidr" folder)
seidr_dir <- deseq2_dir
if (USE_MATRIX == "normalized") {
  M <- t(norm_counts)   # samples x genes
} else {
  M <- t(vst_mat)       # samples x genes
}
write_seidr(M, seidr_dir, drop_nonvarying = DROP_NONVARYING)

msg("Done.")
