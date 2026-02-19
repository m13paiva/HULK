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
  library(grid)
  library(circlize)
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
  make_option(c("--bioproject"), dest = "bioproject", type = "character", default = NULL, help="Single ID or comma-separated list"),
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
  make_option(c("--sample-cor"), dest = "sample_cor", action = "store_true", default = FALSE),
  make_option(c("--dispersion"), dest = "dispersion", action = "store_true", default = FALSE),
  make_option(c("--plots-only"), dest = "plots_only", action = "store_true", default = FALSE),
  make_option(c("--tximport-only"), dest = "tximport_only", action = "store_true", default = FALSE),
  make_option(c("--target-genes"), dest = "target_genes", type = "character", action = "append", default = NULL),
  make_option(c("--mode"), dest = "mode", type = "character", default = "SRR")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- CRITICAL LOCALE FIX ---
Sys.setlocale("LC_NUMERIC", "C")
options(scipen = 999)
options(OutDec = ".")

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# --- MANUAL WRITER FUNCTION ---
write_seidr_manual <- function(mat, outfile) {
  msg("[ManualWrite] Starting manual write to %s", outfile)
  mat[is.na(mat)] <- 0
  con <- file(outfile, "w")
  on.exit(close(con))
  rows <- nrow(mat); cols <- ncol(mat)
  msg("[ManualWrite] Processing %d rows x %d cols...", rows, cols)
  for (i in seq_len(rows)) {
    vals <- mat[i, ]
    vals_fmt <- formatC(vals, format = "f", digits = 6, decimal.mark = ".")
    line_str <- paste(vals_fmt, collapse = "\t")
    cat(line_str, "\n", file = con, sep = "")
    if (i %% 50 == 0) cat(".")
  }
  cat("\n")
  msg("[ManualWrite] Finished writing.")

  msg("[DEBUG] Verifying first 100 characters of output file:")
  raw_head <- readLines(outfile, n = 1, warn = FALSE)
  msg("'%s'", substr(raw_head, 1, 100))
}

# ... [Setup functions] ...
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

# --- GLOBAL VARIABLES ---
# FIX: Handle double nested plots folder
OUT_DIR <- opt$out_dir
if (basename(OUT_DIR) == "plots") {
  OUT_DIR <- dirname(OUT_DIR)
}

SEARCH_DIR    <- opt$search_dir
EXCLUDE_DIRS  <- opt$exclude_dir
TX2GENE_PATH  <- opt$tx2gene

# <--- UPDATED: Handle comma-separated list of BioProjects
BIOPROJECT_IDS <- if(!is.null(opt$bioproject)) unlist(strsplit(opt$bioproject, ",")) else NULL
PER_BP_MODE   <- !is.null(BIOPROJECT_IDS)

TOP_N         <- opt$top_n
DO_PCA        <- isTRUE(opt$pca)
DO_HEATMAP    <- isTRUE(opt$heatmap)
DO_VAR_HM     <- isTRUE(opt$var_heatmap)
DO_SAMPLE_COR <- isTRUE(opt$sample_cor)
DO_DISP       <- isTRUE(opt$dispersion)
PLOTS_ONLY    <- isTRUE(opt$plots_only)
MODE          <- toupper(opt$mode %||% "SRR")

if (MODE == "FASTQ" && DO_VAR_HM) DO_VAR_HM <- FALSE

TARGET_GENES_LIST <- NULL
if (!is.null(opt$target_genes)) {
  TARGET_GENES_LIST <- unlist(strsplit(opt$target_genes, ","))
  msg("[Setup] Target genes detected: %d files/lists provided.", length(TARGET_GENES_LIST))
} else {
  msg("[Setup] No target genes provided (Global mode only).")
}

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
vst_transposed_tsv <- file.path(deseq2_dir, "vst_transposed.tsv")
sf_tsv   <- file.path(deseq2_dir, "size_factors.tsv")
disp_pdf <- file.path(plots_dir, "deseq2_dispersion.pdf")

seidr_genes <- file.path(deseq2_dir, "genes.txt")
seidr_expr  <- file.path(deseq2_dir, "expression.tsv")

pca_pdf    <- file.path(plots_dir, "PCA.pdf")
heatmap_pdf <- file.path(plots_dir, "expression_heatmap_global.pdf")
var_hm_pdf <- file.path(plots_dir, "variance_heatmap.pdf")
scor_pdf   <- file.path(plots_dir, "sample_correlation.pdf")


# --- Plotting Function ---
# UPDATED: Now accepts `vst_top` directly to avoid re-calculation
make_plots <- function(vst_mat, vst_top, info, targets = NULL) {
  if(!DO_PCA && !DO_HEATMAP && !DO_VAR_HM && !DO_SAMPLE_COR) return()

  n_plot <- nrow(vst_top)

  # 1. PCA
  if(DO_PCA && ncol(vst_top)>1) {
    # Use vst_top directly
    pca <- prcomp(t(vst_top[apply(vst_top,1,var)>0,]), center=TRUE)
    df <- as.data.frame(pca$x[,1:2]) %>% tibble::rownames_to_column("sample") %>%
      left_join(info, by="sample")

    df$numeric_id <- seq_len(nrow(df))
    df$legend_label <- paste0(df$numeric_id, ": ", df$sample)
    df$legend_label <- factor(df$legend_label, levels = df$legend_label)

    p <- ggplot(df, aes(x=PC1, y=PC2)) +
         geom_text(aes(label=numeric_id), hjust=-0.3, vjust=-0.3, size=3, show.legend=FALSE) +
         theme_bw() +
         labs(title=sprintf("PCA (Top %d Var Genes)", n_plot))

    if(nrow(df) <= 60) {
       p <- p + geom_point(aes(color=legend_label), size=2) +
            guides(color = guide_legend(ncol = 2, title="Sample Index"))
    } else {
       col_by <- if(PER_BP_MODE || MODE=="FASTQ") "sample" else "bioproject"
       p <- p + geom_point(aes(color=.data[[col_by]]), size=2) +
            labs(subtitle="Legend shows groupings (too many samples for individual index list)")
    }

    ggsave(pca_pdf, p, width=10, height=7)
    readr::write_tsv(df, sub("\\.pdf$", ".tsv", pca_pdf))
    msg("[PCA] Saved %s", pca_pdf)
  }

  # Heatmap Helper
  grp_vec <- if(PER_BP_MODE || MODE=="FASTQ") info$sample else info$bioproject
  title   <- if(PER_BP_MODE || MODE=="FASTQ") "Samples" else "BioProject"

  draw_hm <- function(m, t, f, raster=FALSE) {
    q01 <- quantile(m, 0.01, na.rm = TRUE)
    q50 <- quantile(m, 0.50, na.rm = TRUE)
    q99 <- quantile(m, 0.99, na.rm = TRUE)
    col_fun <- colorRamp2(c(q01, q50, q99), c("blue", "yellow", "red"))

    grDevices::pdf(f, width=12, height=8)
    ht <- ComplexHeatmap::Heatmap(m, name="vst",
      col = col_fun,
      show_row_names=(nrow(m)<150),
      show_column_names=FALSE,
      column_split=grp_vec, column_title=t,
      row_names_gp = gpar(fontsize = 6),
      column_title_gp = gpar(fontsize = 8, fontface="bold"),
      use_raster = raster
    )
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()

    out_df <- as.data.frame(m) %>% tibble::rownames_to_column("gene_id")
    readr::write_tsv(out_df, sub("\\.pdf$", ".tsv", f))
    msg("[Heatmap] Saved %s", f)
  }

  # 2. GLOBAL HEATMAP (Uses vst_top)
  if(DO_HEATMAP) draw_hm(vst_top, "Global Expression (Top N)", heatmap_pdf)

  # 3. TARGETED HEATMAPS (Still searches full vst_mat)
  if(DO_HEATMAP && !is.null(targets)) {
     targets_vec <- as.character(unlist(targets))
     bases <- make.unique(tools::file_path_sans_ext(basename(targets_vec)), sep="_")
     for(i in seq_along(targets_vec)) {
        g <- read_targets(targets_vec[i]); idx <- which(rownames(vst_mat) %in% g)
        if(length(idx) > 0) {
           f <- file.path(plots_dir, sprintf("expression_heatmap_%s.pdf", bases[i]))
           draw_hm(vst_mat[idx, , drop=FALSE], paste("Targeted:", bases[i]), f)
        } else {
           msg("[Heatmap] Warning: No matching genes found for target file %s", targets_vec[i])
        }
     }
  }

  # 4. SAMPLE CORRELATION (Still uses full vst_mat for robustness)
  if(DO_SAMPLE_COR) {
      msg("[SampleCor] Computing Pearson correlation...")
      cor_mat <- cor(vst_mat)

      col_fun_cor <- colorRamp2(c(-1, 0, 1), c("blue", "yellow", "red"))

      grDevices::pdf(scor_pdf, width=10, height=8)
      ht <- ComplexHeatmap::Heatmap(cor_mat,
         name="Pearson",
         col = col_fun_cor,
         show_row_names = (ncol(cor_mat) < 80),
         show_column_names = (ncol(cor_mat) < 80),
         column_title = "Sample-to-Sample Correlation",
         row_names_gp = gpar(fontsize = 6),
         column_names_gp = gpar(fontsize = 6),
         use_raster = TRUE
      )
      ComplexHeatmap::draw(ht)
      grDevices::dev.off()
      readr::write_tsv(as.data.frame(cor_mat) %>% tibble::rownames_to_column("sample"), sub("\\.pdf$", ".tsv", scor_pdf))
      msg("[SampleCor] Saved %s", scor_pdf)
  }

  # 5. VARIANCE HEATMAP (Global uses vst_top)
  if(DO_VAR_HM) {
      draw_var_hm <- function(genes, t, f) {
        # Note: We still access the full matrix for plotting, but `genes` limits it
        df <- t(vst_mat[genes, , drop=FALSE]) %>% as.data.frame() %>%
          tibble::rownames_to_column("sample") %>% left_join(info, by="sample") %>%
          tidyr::pivot_longer(cols = all_of(genes), names_to="gene", values_to="vst")

        var_mat <- df %>% group_by(gene, bioproject) %>%
          summarise(var_vst = stats::var(vst, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = bioproject, values_from = var_vst) %>%
          column_to_rownames("gene") %>% as.matrix()

        na_cols <- colSums(is.na(var_mat)) == nrow(var_mat)
        if (any(na_cols)) var_mat <- var_mat[, !na_cols, drop=FALSE]
        if (ncol(var_mat) == 0) return()

        v_max <- quantile(var_mat, 0.99, na.rm=TRUE)
        col_fun_var <- colorRamp2(c(0, v_max/2, v_max), c("blue", "yellow", "red"))

        grDevices::pdf(f, width = 10, height = 8)
        ht <- ComplexHeatmap::Heatmap(var_mat, name="var(VST)",
            col = col_fun_var,
            show_row_names=(nrow(var_mat)<150),
            show_column_names=TRUE,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6),
            column_title = t
        )
        ComplexHeatmap::draw(ht)
        grDevices::dev.off()
        readr::write_tsv(as.data.frame(var_mat) %>% tibble::rownames_to_column("gene_id"), sub("\\.pdf$", ".tsv", f))
        msg("[VarHeatmap] Saved %s", f)
      }

      msg("[VarHeatmap] Generating Global Variance Heatmap")
      draw_var_hm(rownames(vst_top), "Global Variance", var_hm_pdf)

      # --- TARGETED VARIANCE HEATMAPS ---
      if (!is.null(targets)) {
         targets_vec <- as.character(unlist(targets))
         bases <- make.unique(tools::file_path_sans_ext(basename(targets_vec)), sep="_")
         for(i in seq_along(targets_vec)) {
            g <- read_targets(targets_vec[i]); idx <- which(rownames(vst_mat) %in% g)
            if(length(idx) > 0) {
               f <- file.path(plots_dir, sprintf("variance_heatmap_%s.pdf", bases[i]))
               draw_var_hm(rownames(vst_mat)[idx], paste("Targeted Variance:", bases[i]), f)
            } else {
               msg("[VarHeatmap] Warning: No matching genes found for target file %s", targets_vec[i])
            }
         }
      }
  }
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

files <- find_files(SEARCH_DIR, EXCLUDE_DIRS)
info <- tibble(sample=basename(dirname(files)), bioproject=basename(dirname(dirname(files))))
names(files) <- info$sample

# <--- UPDATED: Check for membership in the provided LIST of IDs, not just equality to a single string
if(PER_BP_MODE) {
  info <- info[info$bioproject %in% BIOPROJECT_IDS,]
  files <- files[info$sample]
}

# --- LOAD / CALCULATE ---

# Shared logic for selecting Top N genes
select_top_genes <- function(mat, n_genes) {
    all_vars <- apply(mat, 1, var, na.rm=TRUE)
    n_select <- min(n_genes, length(all_vars))
    top_idx <- order(all_vars, decreasing=TRUE)[seq_len(n_select)]

    # Return list with both the sub-matrix and the variance dataframe for saving
    list(
        mat = mat[top_idx, , drop=FALSE],
        df  = tibble(gene_id = rownames(mat)[top_idx], variance = all_vars[top_idx])
    )
}

if(PLOTS_ONLY) {
  if (!file.exists(vst_tsv)) stop("VST file missing for fast mode.")
  vst_df <- readr::read_tsv(vst_tsv, col_types=cols())
  if ("sample" %in% colnames(vst_df)) {
      mat <- as.matrix(vst_df[,-1]); rownames(mat) <- vst_df$sample; vst_mat <- t(mat)
  } else {
      mat <- as.matrix(vst_df); rownames(mat) <- info$sample; vst_mat <- t(mat)
  }

  vst_transposed_out <- tibble::rownames_to_column(as.data.frame(vst_mat), "gene_id")
  readr::write_tsv(vst_transposed_out, vst_transposed_tsv)
  msg("[Output] Saved transposed VST matrix to %s", vst_transposed_tsv)

  # --- CALCULATE TOP N ONCE ---
  res <- select_top_genes(vst_mat, TOP_N)
  vst_top <- res$mat

  # Save the Top N file (even in plots-only mode)
  top_var_outfile <- file.path(deseq2_dir, sprintf("top_%d_variable_genes.tsv", TOP_N))
  readr::write_tsv(res$df, top_var_outfile)
  msg("[Analysis] Saved Top %d Variable Genes to %s", nrow(res$df), top_var_outfile)

  make_plots(vst_mat, vst_top, info, TARGET_GENES_LIST)
  quit(status=0)
}

if(file.exists(txi_rds) && !opt$force_txi) {
    msg("[tximport] LOADING CACHED DATA.")
    txi <- readRDS(txi_rds)
} else {
    msg("[tximport] RECALCULATING from %d files...", length(files))
    txi <- tximport(files, type="kallisto", tx2gene=read_tx2gene(TX2GENE_PATH),
           countsFromAbundance=opt$counts_from_abundance, ignoreTxVersion=opt$ignore_tx_version)
}
if(!file.exists(txi_rds)) saveRDS(txi, txi_rds)

msg("[tximport] Saving raw counts to %s", txi_tsv)
readr::write_tsv(as.data.frame(txi$counts) %>% tibble::rownames_to_column("gene"), txi_tsv)

if(opt$tximport_only) quit(status=0)

coldata <- info[match(colnames(txi$counts), info$sample),]
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~1)

# --- FILTERING ---
# Calculate 10% of sample size
min_samples <- round(ncol(dds) * 0.1)
dds <- dds[rowSums(counts(dds) >= 10) >= min_samples, ]
msg("[DESeq2] Estimating size factors (poscounts)...")
dds <- estimateSizeFactors(dds, type="poscounts")

if(DO_DISP) {
    msg("[DESeq2] Estimating Dispersions & Plotting...")
    dds <- estimateDispersions(dds)
    grDevices::pdf(disp_pdf, width=6, height=6)
    plotDispEsts(dds, main="DESeq2 Dispersion Estimates")
    grDevices::dev.off()
    msg("[DispPlot] Saved %s", disp_pdf)
}

vst_mat <- assay(vst(dds, blind=TRUE))

if (opt$drop_nonvarying) {
    msg("[Filter] Removing low-variance genes (threshold < %.2f)...", opt$var_threshold)
    all_vars <- apply(vst_mat, 1, var, na.rm=TRUE)
    keep_genes <- !is.na(all_vars) & (all_vars >= opt$var_threshold)
    vst_mat <- vst_mat[keep_genes, , drop=FALSE]
    msg("[Filter] Genes Remaining: %d", nrow(vst_mat))

    if(nrow(vst_mat) == 0) {
        msg("[Error] No genes left after filtering. Try lowering --var-threshold.")
        quit(status=0)
    }
}

vst_out <- tibble::rownames_to_column(as.data.frame(t(vst_mat)), "sample")
readr::write_tsv(vst_out, vst_tsv)

vst_transposed_out <- tibble::rownames_to_column(as.data.frame(vst_mat), "gene_id")
readr::write_tsv(vst_transposed_out, vst_transposed_tsv)
msg("[Output] Saved transposed VST matrix to %s", vst_transposed_tsv)

# ==============================================================================
# CALCULATE TOP N (MAIN MODE)
# ==============================================================================
msg("[Analysis] Calculating Top %d Most Variable Genes...", TOP_N)
res <- select_top_genes(vst_mat, TOP_N)
vst_top <- res$mat

top_var_outfile <- file.path(deseq2_dir, sprintf("top_%d_variable_genes.tsv", TOP_N))
readr::write_tsv(res$df, top_var_outfile)
msg("[Analysis] Saved Top %d Variable Genes to %s", nrow(res$df), top_var_outfile)


# ==============================================================================
# SEIDR OUTPUT
# ==============================================================================
msg("[Seidr] Saving %d genes for Network Inference...", nrow(vst_mat))
clean_genes <- gsub("\\s+", "_", rownames(vst_mat))
writeLines(clean_genes, seidr_genes)
seidr_mat <- t(vst_mat)
write_seidr_manual(seidr_mat, seidr_expr)

# ==============================================================================
# PLOTS (Passing pre-calculated vst_top)
# ==============================================================================
make_plots(vst_mat, vst_top, coldata, TARGET_GENES_LIST)
msg("Done.")