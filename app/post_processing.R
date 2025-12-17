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
  library(ComplexHeatmap)
  library(grid) # Required for gpar()
})

`%||%` <- function(x, y) if (is.null(x)) y else x

## ===========================
## CLI OPTIONS
## ===========================

option_list <- list(
  make_option(c("-s", "--search-dir"), dest = "search_dir", type = "character", default = "output"),
  make_option(c("-x", "--exclude-dir"), dest = "exclude_dir", type = "character", action = "append", default = c("_mqc_inputs", "multiqc_data", "multiqc", "shared", "logs")),
  make_option(c("-t", "--tx2gene"), dest = "tx2gene", type = "character", default = "tx2gene.clean.tsv"),
  make_option(c("-o", "--out-dir"), dest = "out_dir", type = "character", default = "."),
  make_option(c("-p", "--prefix"), dest = "prefix", type = "character", default = ""),
  make_option(c("--bioproject"), dest = "bioproject", type = "character", default = NULL),
  make_option(c("--ignore-tx-version"), dest = "ignore_tx_version", action = "store_true", default = FALSE),
  make_option(c("--counts-from-abundance"), dest = "counts_from_abundance", type = "character", default = "no"),
  make_option(c("--force-txi"), dest = "force_txi", action = "store_true", default = FALSE),
  make_option(c("--use-matrix"), dest = "use_matrix", type = "character", default = "vst"),
  make_option(c("--no-drop-nonvarying"), dest = "drop_nonvarying", action = "store_false", default = TRUE),
  make_option(c("--var-threshold"), dest = "var_threshold", type = "numeric", default = 0.1),
  make_option(c("--top-n"), dest = "top_n", type = "integer", default = 500),
  make_option(c("--pca"), dest = "pca", action = "store_true", default = FALSE),
  make_option(c("--heatmap"), dest = "heatmap", action = "store_true", default = FALSE),
  make_option(c("--var-heatmap"), dest = "var_heatmap", action = "store_true", default = FALSE),
  make_option(c("--plots-only"), dest = "plots_only", action = "store_true", default = FALSE),
  make_option(c("--tximport-only"), dest = "tximport_only", action = "store_true", default = FALSE),
  make_option(c("--target-genes"), dest = "target_genes", type = "character", default = NULL),
  make_option(c("--mode"), dest = "mode", type = "character", default = "SRR")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Parse Comma-Separated List ---
TARGET_GENES_LIST <- NULL
if (!is.null(opt$target_genes)) {
  raw_list <- strsplit(opt$target_genes, ",")[[1]]
  TARGET_GENES_LIST <- trimws(raw_list)
  TARGET_GENES_LIST <- TARGET_GENES_LIST[TARGET_GENES_LIST != ""]

  message(sprintf("DEBUG: Parsed %d target gene file(s):", length(TARGET_GENES_LIST)))
  for (f in TARGET_GENES_LIST) message(sprintf("  - %s", f))
}

# Map options
SEARCH_DIR    <- opt$search_dir
EXCLUDE_DIRS  <- opt$exclude_dir
TX2GENE_PATH  <- opt$tx2gene
OUT_DIR       <- opt$out_dir
BIOPROJECT_ID <- opt$bioproject
TOP_N         <- opt$top_n
DO_PCA        <- isTRUE(opt$pca)
DO_HEATMAP    <- isTRUE(opt$heatmap)
DO_VAR_HM     <- isTRUE(opt$var_heatmap)
PLOTS_ONLY    <- isTRUE(opt$plots_only)
MODE          <- toupper(opt$mode %||% "SRR")
PER_BP_MODE   <- !is.null(BIOPROJECT_ID)

if (MODE == "FASTQ" && DO_VAR_HM) DO_VAR_HM <- FALSE

# Paths
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
plots_dir  <- file.path(OUT_DIR, "plots")
deseq2_dir <- file.path(OUT_DIR, "deseq2")
tximport_dir <- file.path(OUT_DIR, "tximport")
invisible(lapply(list(plots_dir, deseq2_dir, tximport_dir), dir.create, showWarnings=FALSE, recursive=TRUE))

txi_rds  <- file.path(tximport_dir, "txi.rds")
txi_tsv  <- file.path(tximport_dir, "tximport_counts.tsv")
dds_rds  <- file.path(deseq2_dir, "dds.rds")
vst_tsv  <- file.path(deseq2_dir, "vst.tsv")
norm_tsv <- file.path(deseq2_dir, "normalized_counts.tsv")
sf_tsv   <- file.path(deseq2_dir, "size_factors.tsv")

pca_pdf     <- file.path(plots_dir, "PCA.pdf")
heatmap_pdf <- file.path(plots_dir, "expression_heatmap_global.pdf")
var_hm_pdf  <- file.path(plots_dir, "variance_heatmap.pdf")

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# --- Improved File Finder with Exclusion Logic ---
find_files <- function(dir, excludes) {
  all_files <- list.files(dir, pattern="^abundance\\.tsv$", recursive=TRUE, full.names=TRUE)

  if (length(excludes) > 0) {
      keep <- !vapply(all_files, function(f) {
          any(sapply(excludes, function(e) grepl(e, f)))
      }, logical(1))
      all_files <- all_files[keep]
  }
  return(sort(all_files))
}

read_tx2gene <- function(p) {
  d <- suppressMessages(readr::read_tsv(p, col_types=readr::cols(.default="c")))
  if(ncol(d)<2) stop("tx2gene needs 2 cols")
  colnames(d)[1:2] <- c("transcript_id","gene_id"); d[,1:2]
}
read_targets <- function(p) {
  if(is.null(p)||!file.exists(p)) return(NULL)
  unique(trimws(readLines(p, warn=FALSE)))
}

# --- Plotting Function ---
make_plots <- function(vst_mat, info) {

  # --- DIAGNOSTIC BLOCK ---
  msg("\n--- DIAGNOSTICS: BioProject Alignment ---")

  vst_samples <- colnames(vst_mat)
  info_samples <- info$sample
  common <- intersect(vst_samples, info_samples)

  msg("Samples in VST Matrix (Loaded): %d", length(vst_samples))
  msg("Samples in Directory Scan (On Disk): %d", length(info_samples))
  msg("Overlapping Samples: %d", length(common))

  info <- info[match(common, info$sample), ]
  vst_mat <- vst_mat[, common, drop=FALSE]

  counts <- info %>% group_by(bioproject) %>% summarize(n_samples = n())
  msg("\nBioProjects in Matrix: %d", nrow(counts))

  singletons <- counts %>% filter(n_samples < 2)
  if(nrow(singletons) > 0) {
      msg("\nWARNING: %d BioProjects have < 2 matched samples (undefined variance):", nrow(singletons))
      print(as.data.frame(singletons))
  }
  msg("------------------------------------------\n")
  # --- END DIAGNOSTICS ---

  if(!DO_PCA && !DO_HEATMAP && !DO_VAR_HM) return()

  # Top N
  vars <- apply(vst_mat, 1, var, na.rm=TRUE)
  top <- order(vars, decreasing=TRUE)[seq_len(min(TOP_N, length(vars)))]
  vst_top <- vst_mat[top, , drop=FALSE]

  # --- PCA ---
  if(DO_PCA && ncol(vst_top)>1) {
    pca <- prcomp(t(vst_top[apply(vst_top,1,var)>0,]), center=TRUE)
    df <- as.data.frame(pca$x[,1:2]) %>% tibble::rownames_to_column("sample") %>%
      left_join(info, by="sample")

    col_by <- if(PER_BP_MODE || MODE=="FASTQ") "sample" else "bioproject"
    p <- ggplot(df, aes(x=PC1, y=PC2, color=.data[[col_by]])) +
         geom_point(size=2) + theme_bw() + labs(title="PCA (VST)")

    ggsave(pca_pdf, p, width=6, height=5)
    pca_tsv <- sub("\\.pdf$", ".tsv", pca_pdf)
    readr::write_tsv(df, pca_tsv)
    msg("[PCA] Saved plot to %s and data to %s", pca_pdf, pca_tsv)
  }

  # Heatmap Setup
  if(PER_BP_MODE || MODE=="FASTQ") {
    grp_vec <- info$sample; title <- "Samples"
  } else {
    grp_vec <- info$bioproject; title <- "BioProject"
  }

  draw_hm <- function(m, t, f) {
    grDevices::pdf(f, width=12, height=8) # Increased size slightly

    # VISUAL UPDATE: Font size 6pt
    ht <- ComplexHeatmap::Heatmap(m, name="vst",
      show_row_names=(nrow(m)<150),  # Show more rows if possible
      show_column_names=FALSE,
      column_split=grp_vec, column_title=t,
      row_names_gp = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 8, fontface="bold")
    )
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()

    mat_file <- sub("\\.pdf$", ".tsv", f)
    out_df <- as.data.frame(m) %>% tibble::rownames_to_column("gene_id")
    readr::write_tsv(out_df, mat_file)
    msg("[Heatmap] Saved plot to %s", f)
  }

  if(DO_HEATMAP) {
      msg("[Heatmap] Generating Global Heatmap")
      draw_hm(vst_top, "Global Expression", heatmap_pdf)
  }

  # --- TARGETED HEATMAPS ---
  targets_vec <- NULL
  bases <- NULL
  if (!is.null(TARGET_GENES_LIST)) {
    targets_vec <- as.character(unlist(TARGET_GENES_LIST))
    bases <- make.unique(tools::file_path_sans_ext(basename(targets_vec)), sep="_")
  }

  if(DO_HEATMAP && !is.null(targets_vec)) {
    msg("[TargetPlot] Processing %d target lists for Expression...", length(targets_vec))
    for(i in seq_along(targets_vec)) {
      g <- read_targets(targets_vec[i]); idx <- which(rownames(vst_mat) %in% g)
      if(length(idx) > 0) {
        f <- file.path(plots_dir, sprintf("expression_heatmap_%s.pdf", bases[i]))
        draw_hm(vst_mat[idx, , drop=FALSE], paste("Targeted:", bases[i]), f)
      }
    }
  }

  # --- VARIANCE HEATMAPS ---
  draw_var_hm <- function(genes, t, f) {
    df <- t(vst_mat[genes, , drop=FALSE]) %>% as.data.frame() %>%
      tibble::rownames_to_column("sample") %>% left_join(info, by="sample") %>%
      tidyr::pivot_longer(cols = all_of(genes), names_to="gene", values_to="vst")

    var_mat <- df %>% group_by(gene, bioproject) %>%
      summarise(var_vst = stats::var(vst, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = bioproject, values_from = var_vst) %>%
      column_to_rownames("gene") %>% as.matrix()

    # Detect and remove all-NA columns
    na_cols <- colSums(is.na(var_mat)) == nrow(var_mat)
    if (any(na_cols)) {
        var_mat <- var_mat[, !na_cols, drop=FALSE]
    }

    if (ncol(var_mat) == 0) return()

    grDevices::pdf(f, width = 10, height = 8)

    # VISUAL UPDATE: Font size 6pt for Rows AND Columns
    ht <- ComplexHeatmap::Heatmap(var_mat, name="var(VST)",
        show_row_names=(nrow(var_mat)<150),
        show_column_names=TRUE, # Always show BP names
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        column_title = t
    )
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()

    mat_file <- sub("\\.pdf$", ".tsv", f)
    out_df <- as.data.frame(var_mat) %>% tibble::rownames_to_column("gene_id")
    readr::write_tsv(out_df, mat_file)
    msg("[VarHeatmap] Saved plot to %s", f)
  }

  if(DO_VAR_HM) {
    msg("[VarHeatmap] Generating Global Variance Heatmap")
    draw_var_hm(rownames(vst_top), "Global Variance", var_hm_pdf)

    if (!is.null(targets_vec)) {
      msg("[VarHeatmap] Processing %d target lists for Variance...", length(targets_vec))
      for(i in seq_along(targets_vec)) {
        g <- read_targets(targets_vec[i]); idx <- rownames(vst_mat)[rownames(vst_mat) %in% g]
        if(length(idx) > 0) {
          f <- file.path(plots_dir, sprintf("variance_heatmap_%s.pdf", bases[i]))
          draw_var_hm(idx, paste("Targeted Var:", bases[i]), f)
        }
      }
    }
  }
}

# Main
files <- find_files(SEARCH_DIR, EXCLUDE_DIRS)
info <- tibble(sample=basename(dirname(files)), bioproject=basename(dirname(dirname(files))))
names(files) <- info$sample

if(PER_BP_MODE) {
  info <- info[info$bioproject == BIOPROJECT_ID,]
  files <- files[info$sample]
}

if(PLOTS_ONLY) {
  if (!file.exists(vst_tsv)) stop("VST file missing.")
  vst_df <- readr::read_tsv(vst_tsv, col_types=cols())
  if ("sample" %in% colnames(vst_df)) {
      mat <- as.matrix(vst_df[,-1]); rownames(mat) <- vst_df$sample; vst_mat <- t(mat)
  } else {
      mat <- as.matrix(vst_df); rownames(mat) <- info$sample; vst_mat <- t(mat)
  }

  make_plots(vst_mat, info)
  quit(status=0)
}

# Pipeline
if(file.exists(txi_rds) && !opt$force_txi) {
    msg("[tximport] LOADING CACHED DATA (txi.rds).")
    txi <- readRDS(txi_rds)
} else {
    msg("[tximport] RECALCULATING from %d files found on disk...", length(files))
    txi <- tximport(files, type="kallisto", tx2gene=read_tx2gene(TX2GENE_PATH),
           countsFromAbundance=opt$counts_from_abundance, ignoreTxVersion=opt$ignore_tx_version)
}
if(!file.exists(txi_rds)) saveRDS(txi, txi_rds)

if(opt$tximport_only) {
  readr::write_tsv(as.data.frame(txi$counts) %>% tibble::rownames_to_column("gene"), txi_tsv)
  quit(status=0)
}

coldata <- info[match(colnames(txi$counts), info$sample),]
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~1)
dds <- dds[rowSums(counts(dds))>0,]
msg("[DESeq2] Estimating size factors using type='poscounts'...")
dds <- estimateSizeFactors(dds, type="poscounts")

vst_mat <- assay(vst(dds, blind=TRUE))
vst_out <- tibble::rownames_to_column(as.data.frame(t(vst_mat)), "sample")
readr::write_tsv(vst_out, vst_tsv)

make_plots(vst_mat, coldata)
msg("Done.")