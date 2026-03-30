#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(EGAD)
  library(Matrix)
})

# --- DEBUG LOGGER ---
log_msg <- function(...) {
  message(sprintf(...))
}

option_list <- list(
  make_option(c("-n", "--network"), type="character", default=NULL, help="Network file"),
  make_option(c("-m", "--mapman"), type="character", default=NULL, help="MapMan file"),
  make_option(c("-o", "--output"), type="character", default="egad_results.tsv", help="Output TSV file"),
  make_option(c("-c", "--min_genes"), type="integer", default=10, help="Min genes per term"),
  make_option(c("--curves"), type="character", default=NULL, help="Prefix path to save curve data"),
  make_option(c("--auroc"), action="store_true", default=FALSE, help="Calculate AUROC"),
  make_option(c("--aupr"), action="store_true", default=FALSE, help="Calculate AUPR")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$network) || is.null(opt$mapman)) stop("Files required.")

log_msg("---------------------------------------------------")
log_msg("Starting EGAD Analysis")
log_msg("AUROC requested: %s", opt$auroc)
log_msg("AUPR requested: %s", opt$aupr)

# ------------------------------------------------------------------------------
# 1. MapMan
# ------------------------------------------------------------------------------
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

final_annot <- merge(expanded_dt, bin_names, by.x = "expanded_bins", by.y = "BINCODE", all.x = TRUE)
final_annot[is.na(NAME), NAME := expanded_bins]

annot_matrix_long <- final_annot[, .(IDENTIFIER, NAME)]
annot_matrix_long[, val := 1]

annot_sparse <- Matrix::sparseMatrix(
  i = as.numeric(as.factor(annot_matrix_long$IDENTIFIER)),
  j = as.numeric(as.factor(annot_matrix_long$NAME)),
  x = annot_matrix_long$val,
  dimnames = list(levels(as.factor(annot_matrix_long$IDENTIFIER)),
                  levels(as.factor(annot_matrix_long$NAME)))
)

# ------------------------------------------------------------------------------
# 2. Network
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

rm(edges)
gc()

net_mat_general <- Matrix::sparseMatrix(
  i = edges_unique$i,
  j = edges_unique$j,
  x = edges_unique$weight,
  dims = c(n_genes, n_genes),
  dimnames = list(all_genes, all_genes),
  symmetric = FALSE
)
net_mat_sparse <- Matrix::forceSymmetric(net_mat_general, uplo = "U")

rm(edges_unique, net_mat_general)
gc()

# ------------------------------------------------------------------------------
# 3. Intersect & Prep
# ------------------------------------------------------------------------------
log_msg(">> Intersecting Genes...")
common_genes <- intersect(rownames(net_mat_sparse), rownames(annot_sparse))

if (length(common_genes) < 10) {
    log_msg("!!! CRITICAL: Overlap too small (<10). Returning empty result.")
    fwrite(data.table(Term="None", AUC=NA, AUPR=NA, NodeDegreeAUC=NA), opt$output, sep="\t")
    quit(save="no", status=0)
}

net_final_sparse <- net_mat_sparse[common_genes, common_genes]
annot_final_sparse <- annot_sparse[common_genes, ]

term_sizes <- colSums(annot_final_sparse)
valid_terms <- term_sizes >= opt$min_genes

if (sum(valid_terms) < 2) {
    log_msg("!!! CRITICAL: No valid terms remaining.")
    fwrite(data.table(Term="None", AUC=NA, AUPR=NA, NodeDegreeAUC=NA), opt$output, sep="\t")
    quit(save="no", status=0)
}

annot_final_sparse <- annot_final_sparse[, valid_terms, drop = FALSE]

net_final <- as.matrix(net_final_sparse)
annot_final <- as.matrix(annot_final_sparse)
mode(annot_final) <- "numeric"

net_final[is.na(net_final)] <- 0
diag(net_final) <- 0

# ------------------------------------------------------------------------------
# 4. Run EGAD
# ------------------------------------------------------------------------------
start_time <- Sys.time()
term_names <- colnames(annot_final)
res_auc <- data.table(Term = term_names)
res_aupr <- data.table(Term = term_names)

tryCatch({
    if (opt$aupr) {
        log_msg(">> Running EGAD (AUPR mode)...")
        nv_prc <- neighbor_voting(genes.labels = annot_final, network = net_final, nFold = 3, output = "PR")
        res_aupr <- data.table(Term = rownames(nv_prc), AUPR = as.numeric(nv_prc[, 1]))
    }

    if (opt$auroc) {
        log_msg(">> Running EGAD (AUROC mode)...")
        nv_res <- neighbor_voting(genes.labels = annot_final, network = net_final, nFold = 3, output = "AUROC")
        res_auc <- data.table(Term = rownames(nv_res), AUC = as.numeric(nv_res[, 1]), NodeDegreeAUC = as.numeric(nv_res[, 2]))
    }

    if (!is.null(opt$curves)) {
        log_msg(">> Generating raw network predictions for curves...")
        preds <- EGAD::predictions(genes.labels = annot_final, network = net_final)
    }
}, error = function(e) {
    log_msg("!!! EGAD CRASHED: %s", e$message)
    quit(save="no", status=1)
})

# Merge results based on whatever was requested
results_dt <- merge(res_auc, res_aupr, by="Term", all=TRUE)

# Only print means for what was actually calculated
if (opt$auroc) {
    valid_auc <- results_dt[!is.na(AUC)]
    log_msg(sprintf(">> Mean AUROC: %.4f", ifelse(nrow(valid_auc) > 0, mean(valid_auc$AUC), NA)))
}
if (opt$aupr) {
    valid_aupr <- results_dt[!is.na(AUPR)]
    log_msg(sprintf(">> Mean AUPR: %.4f", ifelse(nrow(valid_aupr) > 0, mean(valid_aupr$AUPR), NA)))
}

fwrite(results_dt, opt$output, sep="\t")

# Output Curve Data if requested
if (!is.null(opt$curves) && exists("preds")) {
    log_msg(">> Exporting micro-average curve data & baseline...")
    flat_preds <- as.vector(preds)
    flat_labels <- as.vector(annot_final)

    valid_idx <- !is.na(flat_preds) & !is.na(flat_labels)
    flat_preds <- flat_preds[valid_idx]
    flat_labels <- flat_labels[valid_idx]

    roc_data <- get_roc(flat_preds, flat_labels)
    prc_data <- get_prc(flat_preds, flat_labels)

    fwrite(data.table(FPR=roc_data[,1], TPR=roc_data[,2]), paste0(opt$curves, "_roc.tsv"), sep="\t")
    fwrite(data.table(Recall=prc_data[,1], Precision=prc_data[,2]), paste0(opt$curves, "_prc.tsv"), sep="\t")

    # Calculate global random baseline AUPR
    baseline_pr <- sum(flat_labels) / length(flat_labels)
    fwrite(data.table(Baseline=baseline_pr), paste0(opt$curves, "_baseline.tsv"), sep="\t")
}

end_time <- Sys.time()
log_msg(">> Done. (Time: %.2fs)", as.numeric(difftime(end_time, start_time, units="secs")))