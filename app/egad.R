#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(EGAD)
  library(Matrix)
})

# --- DEBUG LOGGER ---
log_msg <- function(...) {
  print(sprintf(...))
}

option_list <- list(
  make_option(c("-n", "--network"), type="character", default=NULL, help="Network file"),
  make_option(c("-m", "--mapman"), type="character", default=NULL, help="MapMan file"),
  make_option(c("-g", "--go_file"), type="character", default=NULL, help="BioMart GO export file (TSV)"),
  make_option(c("-o", "--output"), type="character", default="egad_results.tsv", help="Output TSV file"),
  make_option(c("-c", "--min_genes"), type="integer", default=10, help="Min genes per term"),
  make_option(c("--curves"), type="character", default=NULL, help="Prefix path to save curve data"),
  make_option(c("--auroc"), action="store_true", default=FALSE, help="Calculate AUROC"),
  make_option(c("--aupr"), action="store_true", default=FALSE, help="Calculate AUPR")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$network)) stop("Network file (-n) is required.")

use_mapman <- !is.null(opt$mapman)
use_go <- !is.null(opt$go_file)

if (!use_mapman && !use_go) stop("You must provide either a MapMan file (-m) or a GO BioMart file (-g).")

log_msg("---------------------------------------------------")
log_msg("Starting EGAD Analysis (Micro & Macro Mode)")

# ------------------------------------------------------------------------------
# 1. Loading Network
# ------------------------------------------------------------------------------
log_msg(">> Loading Network...")
edges <- fread(opt$network)

if (ncol(edges) >= 3) {
  colnames(edges)[1:3] <- c("source", "target", "weight")
} else {
  stop("Network file must have at least 3 columns.")
}

edges[, source := toupper(source)]
edges[, target := toupper(target)]
edges[, weight := as.numeric(weight)]
edges <- edges[!is.na(weight)]

all_genes <- unique(c(edges$source, edges$target))
n_genes <- length(all_genes)
gene_map <- setNames(1:n_genes, all_genes)

edges[, idx1 := gene_map[source]]
edges[, idx2 := gene_map[target]]
edges[, i := pmin(idx1, idx2)]
edges[, j := pmax(idx1, idx2)]
edges_unique <- unique(edges, by = c("i", "j"))
rm(edges); gc()

net_sparse_init <- Matrix::sparseMatrix(
  i = edges_unique$i,
  j = edges_unique$j,
  x = edges_unique$weight,
  dims = c(n_genes, n_genes),
  dimnames = list(all_genes, all_genes),
  symmetric = TRUE
)
rm(edges_unique); gc()

# ------------------------------------------------------------------------------
# 2. Annotation Parsing
# ------------------------------------------------------------------------------
annotations_list <- list()

if (use_mapman) {
  log_msg(">> Parsing MapMan file...")
  mm <- fread(opt$mapman, header = TRUE, quote = "", fill = TRUE, sep = "\t")
  colnames(mm) <- gsub("'", "", colnames(mm))
  if (!"IDENTIFIER" %in% colnames(mm)) {
    if(ncol(mm) >= 3) {
      mm <- mm[, 1:3]
      colnames(mm) <- c("BINCODE", "NAME", "IDENTIFIER")
    }
  }
  mm_clean <- mm[IDENTIFIER != "" & IDENTIFIER != "''"]
  mm_clean[, IDENTIFIER := toupper(gsub("'", "", IDENTIFIER))]
  mm_clean[, BINCODE := gsub("'", "", BINCODE)]
  mm_clean[, NAME := gsub("'", "", NAME)]

  unique_pairs <- unique(mm_clean[, .(IDENTIFIER, BINCODE)])
  get_parents <- function(bincode) {
    parts <- unlist(strsplit(bincode, "\\."))
    sapply(1:length(parts), function(i) paste(parts[1:i], collapse = "."))
  }

  expanded_dt <- unique_pairs[, .(expanded_bins = unlist(lapply(BINCODE, get_parents))), by = IDENTIFIER]
  bin_names <- unique(mm[, .(BINCODE = gsub("'", "", BINCODE), NAME = gsub("'", "", NAME))])
  final_annot_dt <- merge(expanded_dt, bin_names, by.x = "expanded_bins", by.y = "BINCODE", all.x = TRUE)
  final_annot_dt[is.na(NAME), NAME := expanded_bins]

  annot_matrix_long <- final_annot_dt[, .(IDENTIFIER, NAME)]
  setnames(annot_matrix_long, c("IDENTIFIER", "TERM"))
  annot_matrix_long[, val := 1]

  annotations_list[["MapMan"]] <- Matrix::sparseMatrix(
    i = as.numeric(as.factor(annot_matrix_long$IDENTIFIER)),
    j = as.numeric(as.factor(annot_matrix_long$TERM)),
    x = annot_matrix_long$val,
    dimnames = list(levels(as.factor(annot_matrix_long$IDENTIFIER)),
                    levels(as.factor(annot_matrix_long$TERM)))
  )
  rm(mm, mm_clean, expanded_dt, bin_names, final_annot_dt, annot_matrix_long); gc()
}

if (use_go) {
  log_msg(">> Parsing BioMart GO file...")
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("Missing GO.db")
  suppressPackageStartupMessages(library(GO.db))

  go_dt <- fread(opt$go_file, header=TRUE, sep="\t", quote="", fill=TRUE)
  if (ncol(go_dt) < 2) stop("BioMart file must have at least 2 columns.")
  setnames(go_dt, 1:2, c("IDENTIFIER", "GO_ID"))
  go_dt <- go_dt[GO_ID != "" & !is.na(GO_ID)]
  go_dt[, IDENTIFIER := toupper(IDENTIFIER)]

  unique_go <- unique(go_dt$GO_ID)
  bp_anc <- as.list(GOBPANCESTOR); mf_anc <- as.list(GOMFANCESTOR); cc_anc <- as.list(GOCCANCESTOR)

  go_anc_list <- lapply(unique_go, function(id) {
    ancs <- c(id)
    if (id %in% names(bp_anc)) ancs <- c(ancs, bp_anc[[id]])
    if (id %in% names(mf_anc)) ancs <- c(ancs, mf_anc[[id]])
    if (id %in% names(cc_anc)) ancs <- c(ancs, cc_anc[[id]])
    return(unique(ancs[!is.na(ancs) & ancs != "all"]))
  })
  names(go_anc_list) <- unique_go

  annot_matrix_long <- go_dt[, .(TERM = unlist(go_anc_list[GO_ID])), by = IDENTIFIER]
  annot_matrix_long <- unique(annot_matrix_long)
  annot_matrix_long[, val := 1]

  annotations_list[["GO"]] <- Matrix::sparseMatrix(
    i = as.numeric(as.factor(annot_matrix_long$IDENTIFIER)),
    j = as.numeric(as.factor(annot_matrix_long$TERM)),
    x = annot_matrix_long$val,
    dimnames = list(levels(as.factor(annot_matrix_long$IDENTIFIER)),
                    levels(as.factor(annot_matrix_long$TERM)))
  )
  rm(go_dt, annot_matrix_long, go_anc_list, bp_anc, mf_anc, cc_anc); gc()
}

# ------------------------------------------------------------------------------
# 3. Execution Loop
# ------------------------------------------------------------------------------
final_results <- list()
start_time_global <- Sys.time()

for (anno_source in names(annotations_list)) {
  log_msg("===================================================")
  log_msg(">> Processing Annotation Source: %s", anno_source)

  annot_sparse_init <- annotations_list[[anno_source]]
  common_genes <- intersect(rownames(net_sparse_init), rownames(annot_sparse_init))

  if (length(common_genes) < 10) next

  net_final_sparse <- net_sparse_init[common_genes, common_genes]
  annot_final_sparse <- annot_sparse_init[common_genes, ]

  valid_terms <- colSums(annot_final_sparse) >= opt$min_genes
  if (sum(valid_terms) < 2) next

  annot_final_sparse <- annot_final_sparse[, valid_terms, drop = FALSE]

  net_final <- as.matrix(net_final_sparse)
  annot_final <- as.matrix(annot_final_sparse)
  rm(net_final_sparse, annot_final_sparse); gc()

  mode(annot_final) <- "numeric"
  diag(net_final) <- 0
  net_final[is.na(net_final)] <- 0

  term_names <- colnames(annot_final)
  res_auc <- data.table(Term = term_names)
  res_aupr <- data.table(Term = term_names)

  tryCatch({
    if (opt$aupr) {
      log_msg(">> Running EGAD (AUPR mode)...")
      preds_aupr <- EGAD::predictions(genes.labels = annot_final, network = net_final)
      aupr_vals <- sapply(1:ncol(annot_final), function(i) {
        pr_curve <- get_prc(preds_aupr[, i], annot_final[, i])
        pr_curve <- pr_curve[order(pr_curve[, 1]), ]
        sum(diff(pr_curve[, 1]) * (pr_curve[-1, 2] + pr_curve[-nrow(pr_curve), 2]) / 2, na.rm = TRUE)
      })
      res_aupr <- data.table(Term = term_names, AUPR = aupr_vals)
    }

    if (opt$auroc) {
      log_msg(">> Running EGAD (AUROC mode)...")
      nv_res <- neighbor_voting(genes.labels = annot_final, network = net_final, nFold = 3, output = "AUROC")
      res_auc <- data.table(Term = rownames(nv_res), AUC = as.numeric(nv_res[, 1]), NodeDegreeAUC = as.numeric(nv_res[, 2]))
    }

    if (!is.null(opt$curves)) {
      log_msg(">> Generating MICRO-AVERAGED curve data...")
      if (!exists("preds_aupr")) preds <- EGAD::predictions(genes.labels = annot_final, network = net_final) else preds <- preds_aupr

      flat_preds <- as.vector(preds)
      flat_labels <- as.vector(annot_final)
      valid_idx <- !is.na(flat_preds) & !is.na(flat_labels)

      roc_data <- get_roc(flat_preds[valid_idx], flat_labels[valid_idx])
      prc_data <- get_prc(flat_preds[valid_idx], flat_labels[valid_idx])
      baseline_val <- sum(flat_labels[valid_idx]) / length(flat_labels[valid_idx])

      # Calculate Micro Integrations
      micro_auc <- sum(diff(roc_data[,1]) * (roc_data[-1,2] + roc_data[-nrow(roc_data),2]) / 2, na.rm = TRUE)
      prc_ordered <- prc_data[order(prc_data[,1]), ]
      micro_aupr <- sum(diff(prc_ordered[,1]) * (prc_ordered[-1,2] + prc_ordered[-nrow(prc_ordered),2]) / 2, na.rm = TRUE)

      if (opt$auroc) log_msg(">> %s [MICRO] Averaged AUROC: %.4f", anno_source, micro_auc)
      if (opt$aupr)  log_msg(">> %s [MICRO] Averaged AUPR: %.4f", anno_source, micro_aupr)

      fwrite(data.table(FPR=roc_data[,1], TPR=roc_data[,2]), paste0(opt$curves, "_", anno_source, "_roc.tsv"), sep="\t")
      fwrite(data.table(Recall=prc_data[,1], Precision=prc_data[,2]), paste0(opt$curves, "_", anno_source, "_prc.tsv"), sep="\t")
      fwrite(data.table(Baseline=baseline_val), paste0(opt$curves, "_", anno_source, "_baseline.tsv"), sep="\t")
    }
  }, error = function(e) {
    log_msg("!!! EGAD CRASHED on %s: %s", anno_source, e$message)
    return(NULL)
  })

  results_dt <- merge(res_auc, res_aupr, by="Term", all=TRUE)
  if (nrow(results_dt) > 0) {
    results_dt[, Annotation_Source := anno_source]
    final_results[[anno_source]] <- results_dt
    if (opt$auroc) log_msg(">> %s [MACRO] Averaged AUROC: %.4f", anno_source, mean(results_dt$AUC, na.rm=TRUE))
    if (opt$aupr)  log_msg(">> %s [MACRO] Averaged AUPR: %.4f", anno_source, mean(results_dt$AUPR, na.rm=TRUE))
  }

  rm(net_final, annot_final)
  if (exists("preds")) rm(preds)
  if (exists("preds_aupr")) rm(preds_aupr)
  gc()
}

if (length(final_results) > 0) {
  combined_dt <- rbindlist(final_results, use.names=TRUE, fill=TRUE)
  col_order <- c("Term", "Annotation_Source", setdiff(names(combined_dt), c("Term", "Annotation_Source")))
  setcolorder(combined_dt, col_order)
  fwrite(combined_dt, opt$output, sep="\t")
} else {
  fwrite(data.table(Term="None", Annotation_Source="None", AUC=NA, AUPR=NA, NodeDegreeAUC=NA), opt$output, sep="\t")
}
log_msg("===================================================")
log_msg(">> All Done. (Total Time: %.2fs)", as.numeric(difftime(Sys.time(), start_time_global, units="secs")))