#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(EGAD)
  library(Matrix)
})

# --- DEBUG LOGGER ---
# We use message() so it goes to stderr, which Python captures and prints
log_msg <- function(...) {
  message(sprintf(...))
}

option_list <- list(
  make_option(c("-n", "--network"), type="character", default=NULL, help="Network file"),
  make_option(c("-m", "--mapman"), type="character", default=NULL, help="MapMan file"),
  # HULK passes a specific FILE path for output, not a directory
  make_option(c("-o", "--output"), type="character", default="egad_results.tsv", help="Output TSV file"),
  make_option(c("-c", "--min_genes"), type="integer", default=5, help="Min genes per term")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$network) || is.null(opt$mapman)) stop("Files required.")

log_msg("---------------------------------------------------")
log_msg("Starting EGAD Analysis (Robust Version)")
log_msg("Network: %s", opt$network)
log_msg("MapMan:  %s", opt$mapman)
log_msg("Output:  %s", opt$output)

# ------------------------------------------------------------------------------
# 1. MapMan
# ------------------------------------------------------------------------------
log_msg(">> Parsing MapMan file...")
# Fill=TRUE handles irregular rows
mm <- fread(opt$mapman, header = TRUE, quote = "", fill = TRUE, sep = "\t")

# Clean column names
colnames(mm) <- gsub("'", "", colnames(mm))

# Fallback column detection if headers aren't standard
if (!"IDENTIFIER" %in% colnames(mm)) {
    # Try guessing based on position (Mercator standard: Bin=1, Name=2, ID=3)
    if(ncol(mm) >= 3) {
        log_msg("   [Info] Non-standard headers. Assuming Cols 1=Bin, 2=Name, 3=ID")
        mm <- mm[, 1:3]
        colnames(mm) <- c("BINCODE", "NAME", "IDENTIFIER")
    }
}

mm_clean <- mm[IDENTIFIER != "" & IDENTIFIER != "''"]

# Cleanup values
mm_clean[, IDENTIFIER := toupper(gsub("'", "", IDENTIFIER))]
mm_clean[, BINCODE := gsub("'", "", BINCODE)]
mm_clean[, NAME := gsub("'", "", NAME)]

# Ontology Expansion (Parents)
log_msg("   [Info] Expanding ontology hierarchy...")
unique_pairs <- unique(mm_clean[, .(IDENTIFIER, BINCODE)])
get_parents <- function(bincode) {
  parts <- unlist(strsplit(bincode, "\\."))
  sapply(1:length(parts), function(i) paste(parts[1:i], collapse = "."))
}

# This step can be slow, but essential for correct GO-like behavior
expanded_dt <- unique_pairs[, .(expanded_bins = unlist(lapply(BINCODE, get_parents))), by = IDENTIFIER]
bin_names <- unique(mm[, .(BINCODE = gsub("'", "", BINCODE), NAME = gsub("'", "", NAME))])

final_annot <- merge(expanded_dt, bin_names, by.x = "expanded_bins", by.y = "BINCODE", all.x = TRUE)
# Fill missing names with BinCode if name missing
final_annot[is.na(NAME), NAME := expanded_bins]

# Initial Matrix
annot_matrix_long <- final_annot[, .(IDENTIFIER, NAME)]
annot_matrix_long[, val := 1]

# Create Sparse Annotation Matrix
annot_sparse <- Matrix::sparseMatrix(
  i = as.numeric(as.factor(annot_matrix_long$IDENTIFIER)),
  j = as.numeric(as.factor(annot_matrix_long$NAME)),
  x = annot_matrix_long$val,
  dimnames = list(levels(as.factor(annot_matrix_long$IDENTIFIER)),
                  levels(as.factor(annot_matrix_long$NAME)))
)
log_msg("   [Info] Annotation Matrix: %d genes x %d terms", nrow(annot_sparse), ncol(annot_sparse))

# ------------------------------------------------------------------------------
# 2. Network
# ------------------------------------------------------------------------------
log_msg(">> Loading Network...")
edges <- fread(opt$network)

# Handle Seidr's variable headers
if (ncol(edges) >= 3) {
    # Assume 1=Source, 2=Target, 3=Weight if headers missing/weird
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

# Map to integer indices
edges[, idx1 := gene_map[source]]
edges[, idx2 := gene_map[target]]
edges[, i := pmin(idx1, idx2)]
edges[, j := pmax(idx1, idx2)]
edges_unique <- unique(edges, by = c("i", "j"))

rm(edges)
gc()

log_msg(">> Constructing Sparse Network Matrix...")
net_mat_general <- Matrix::sparseMatrix(
  i = edges_unique$i,
  j = edges_unique$j,
  x = edges_unique$weight,
  dims = c(n_genes, n_genes),
  dimnames = list(all_genes, all_genes),
  symmetric = FALSE
)
# Safe symmetric forcing
net_mat_sparse <- Matrix::forceSymmetric(net_mat_general, uplo = "U")

rm(edges_unique, net_mat_general)
gc()

log_msg("   [Info] Network Matrix: %d x %d", nrow(net_mat_sparse), ncol(net_mat_sparse))

# ------------------------------------------------------------------------------
# 3. Intersect & Prep
# ------------------------------------------------------------------------------
log_msg(">> Intersecting Genes...")
common_genes <- intersect(rownames(net_mat_sparse), rownames(annot_sparse))
log_msg("   [Info] Intersection size: %d genes", length(common_genes))

if (length(common_genes) < 10) {
    log_msg("!!! CRITICAL: Overlap too small (<10). Returning empty result.")
    fwrite(data.table(Term="None", AUC=0.5, NodeDegreeAUC=0.5), opt$output, sep="\t")
    quit(save="no", status=0)
}

# Subset
net_final_sparse <- net_mat_sparse[common_genes, common_genes]
annot_final_sparse <- annot_sparse[common_genes, ]

# Filter Terms
log_msg(">> Filtering terms based on intersection overlap (min_genes=%d)...", opt$min_genes)
term_sizes <- colSums(annot_final_sparse)
valid_terms <- term_sizes >= opt$min_genes
log_msg("   [Info] Terms retained: %d (out of %d)", sum(valid_terms), length(valid_terms))

if (sum(valid_terms) < 2) {
    log_msg("!!! CRITICAL: No valid terms remaining.")
    fwrite(data.table(Term="None", AUC=0.5, NodeDegreeAUC=0.5), opt$output, sep="\t")
    quit(save="no", status=0)
}

annot_final_sparse <- annot_final_sparse[, valid_terms, drop = FALSE]

# Convert to Dense (EGAD neighbor_voting usually prefers dense or specific sparse classes)
log_msg(">> Converting to Dense Matrix for EGAD...")
net_final <- as.matrix(net_final_sparse)
annot_final <- as.matrix(annot_final_sparse)
mode(annot_final) <- "numeric"

# Zero out NAs and Diagonal
net_final[is.na(net_final)] <- 0
diag(net_final) <- 0

log_msg("   [Info] Final evaluation set: %d genes, %d terms.", nrow(net_final), ncol(annot_final))

# ------------------------------------------------------------------------------
# 4. Run EGAD
# ------------------------------------------------------------------------------
log_msg(">> Running EGAD (neighbor_voting)...")
start_time <- Sys.time()

# Catch potential crashes in EGAD itself
tryCatch({
    nv_res <- neighbor_voting(genes.labels = annot_final, network = net_final, nFold = 3, output = "AUROC")
}, error = function(e) {
    log_msg("!!! EGAD CRASHED: %s", e$message)
    quit(save="no", status=1)
})

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs")

# [DATA FORENSICS]
auroc_scores <- numeric(0)
null_scores <- numeric(0)
term_names <- character(0)

if (is.list(nv_res) && !is.data.frame(nv_res) && !is.matrix(nv_res)) {
  # Scenario A: List
  auroc_scores <- as.numeric(nv_res$AUC)
  null_scores  <- as.numeric(nv_res$deg_AUC)
  term_names   <- names(nv_res$AUC)

} else if (is.matrix(nv_res)) {
  # Scenario B: Matrix
  cols <- colnames(nv_res)

  # Grab AUROC
  auc_col <- grep("auc", cols, ignore.case=TRUE, value=TRUE)
  auc_col <- auc_col[!grepl("null|degree", auc_col, ignore.case=TRUE)]

  if (length(auc_col) > 0) {
    auroc_scores <- as.numeric(nv_res[, auc_col[1]])
  } else {
    log_msg("   [Warn] Could not identify 'auc' column. Using 1st column.")
    auroc_scores <- as.numeric(nv_res[, 1])
  }

  # Grab Null
  null_col <- grep("null", cols, ignore.case=TRUE, value=TRUE)
  if (length(null_col) > 0) {
    null_scores <- as.numeric(nv_res[, null_col[1]])
  } else {
    null_scores <- rep(NA, length(auroc_scores))
  }

  term_names <- rownames(nv_res)

} else if (is.vector(nv_res)) {
  # Scenario C: Vector
  auroc_scores <- as.numeric(nv_res)
  null_scores  <- rep(NA, length(nv_res))
  term_names   <- names(nv_res)

} else {
  stop("EGAD returned an unrecognizable object type.")
}

if (is.null(term_names)) term_names <- paste0("Term_", 1:length(auroc_scores))

# Create Result Table
# We use 'AUC' as the column name to match what the Python wrapper expects
results_dt <- data.table(
  Term = term_names,
  AUC = auroc_scores,
  NodeDegreeAUC = null_scores
)

# Remove NA results
valid_results <- results_dt[!is.na(AUC)]

if (nrow(valid_results) == 0) {
  log_msg(">> [FAILURE] All terms resulted in NA/NaN.")
  mean_auroc <- 0.5
} else {
  mean_auroc <- mean(valid_results$AUC, na.rm = TRUE)
}

log_msg(sprintf("\n>> RESULTS: Mean AUROC: %.4f (Time: %.2fs)", mean_auroc, as.numeric(duration)))
log_msg(sprintf("   Terms Evaluated: %d", nrow(valid_results)))

# Write to the specific output file requested by Python wrapper
fwrite(results_dt, opt$output, sep="\t")
log_msg(">> Done.")