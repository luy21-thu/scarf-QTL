# ================================================================
# This script fits null gene-expression models (no SNP testing).
# For a given cell type (CT_idx) and chromosome (chr), this script:
#   (i) loads pre-split SCT expression (genes × cells) and covariates,
#   (ii) runs scarfQTL with Test_type = NULL to fit the null model only,
#   (iii) saves per-gene eta_hat0 to scarfQTL/MLE0/<ct>_chr<chr>.rds.
# Usage: Rscript 2_fit_null_hypothesis.R <CT_idx> <chr>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")
library(data.table)
library(dplyr)
library(stringr)

# Parse command-line arguments and cell-type labels
{
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 2)
  CT_idx <- as.integer(args[1])
  chr    <- as.integer(args[2])
  stopifnot(!is.na(CT_idx), !is.na(chr), CT_idx >= 1, CT_idx <= 14, chr >= 1, chr <= 22)
  
  # Example
  # CT_idx <- 14; chr <- 1
  CT_14 <- c("B IN", "B Mem", "CD4 ET", "CD4 NC", "CD4 SOX4", 
             "CD8 ET", "CD8 NC", "CD8 S100B", "DC", 
             "Mono C", "Mono NC", "NK", "NK R", "Plasma")
  CT <- CT_14[CT_idx]
  ct <- str_replace(tolower(CT), " ", "")
  
  DIR_MLE0 <- "results/realdata/intermediate/MLE0"
  dir.create(DIR_MLE0, recursive = TRUE, showWarnings = FALSE)
  out_file <- sprintf("%s/%s_chr%s.rds", DIR_MLE0, ct, chr)
}

if (!file.exists(out_file)) {
  # Read input data
  {
    # Cell-level and individual-level covariates
    {
      cell_meta <- readRDS(sprintf(
        "data/covariate/cell_pt_%s.rds", ct))
      INDI_C <- cell_meta$individual; INDI <- sort(unique(INDI_C))
      TIME <- cell_meta$pseudotime
      individual_meta <- readRDS(sprintf(
        "data/covariate/individual_%s.rds", ct))
      
      W_INDI_ <- t(scale(individual_meta[INDI, c(
        "sex", "age", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pf1", "pf2")]))
      W_CELL <- t(scale(cell_meta[,c("log_nUMI", "percent.mt")]))
      #Z_INDI_ <- matrix(1, ncol=length(INDI), dimnames=list("baseline", INDI))
      #Z_CELL <- NULL
    }
    # Genetic relatedness matrix
    {
      ids <- fread("data/genotype/GRM.rel.id", header = T)
      colnames(ids) <- c("FID", "IID")
      GRM <- as.matrix(fread("data/genotype/GRM.rel"))
      colnames(GRM) <- rownames(GRM) <- ids$IID
      GRM_INDI <- GRM[INDI, INDI]
      std_devs <- sqrt(diag(GRM_INDI))
      GRM_INDI <- GRM_INDI / (std_devs %*% t(std_devs))
    }
    # Gene expression
    {
      Y_mat <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
      Gene_CT <- readRDS(sprintf(
        "data/Gene_cis_SNP/Gene_%s.rds", ct))[[chr]]
      # Gene_SNP_list is only used to define the set of genes here
      Gene_SNP_list <- setNames(vector("list", length(Gene_CT)), Gene_CT)
      stopifnot(all(Gene_CT %in% rownames(Y_mat)))
    }
  }
  
  # Fit null model
  {
    cat(sprintf("[%s | chr%d] fitting null models for %d genes\n",
                ct, chr, length(Gene_SNP_list)))
    res_0 <- scarfQTL(
      Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
      W_INDI_ = W_INDI_, W_CELL = W_CELL, 
      # Z_INDI_ = Z_INDI_, Z_CELL = Z_CELL,
      TIME = TIME, 
      Test_type = NULL, Gene_SNP_list = Gene_SNP_list, 
      max_iter = 5000)
    eta_hat0 <- res_0$eta_hat0
  }
  # Save results
  {
    tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
    saveRDS(eta_hat0, file = tmp_file)
    file.rename(tmp_file, out_file)
  }
} else {
  cat(sprintf("[%s chr%d] null model exists, skipped.\n", ct, chr))
}

