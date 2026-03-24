# ================================================================
# This script estimates dynamic eQTL effect-size curves for significant signals.
# For each chromosome (chr = 1..22) and each of 14 immune cell types, this script:
#   (i) loads pre-split SCT expression (genes × cells), covariates,
#       significant SNP lists, per-gene null-fit eta_hat, and per-chr genotypes;
#   (ii) runs scarfQTL with Estimation = TRUE on significant eQTLs only;
#   (iii) saves estimated effect-size curves under scarf_QTL/Estimation/.
# Usage: Rscript 4_estimation.R <chr>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")
library(data.table)
library(dplyr)
library(stringr)

# Parse command-line arguments and cell-type labels
{
  args <- commandArgs(trailingOnly = TRUE)
  chr <- as.numeric(args[1])
  stopifnot(!is.na(chr), chr >= 1, chr <= 22)
  
  # Example
  # chr <- 1
  CT_14 <- c("B IN", "B Mem", "CD4 ET", "CD4 NC", "CD4 SOX4", 
             "CD8 ET", "CD8 NC", "CD8 S100B", "DC", 
             "Mono C", "Mono NC", "NK", "NK R", "Plasma")
  
  DIR_MLE0 <- "results/realdata/intermediate/MLE0"
  DIR_Test <- "results/realdata/intermediate/Test"
  DIR_Est <- "results/realdata/intermediate/Estimation"
  dir.create(DIR_Est, recursive = TRUE, showWarnings = FALSE)
}

# Read genotype data
G_SNP_mat <- readRDS(sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))

for (CT_idx in seq(14)) {
  CT <- CT_14[CT_idx]
  ct <- str_replace(tolower(CT), " ", "")
  
  out_file <- sprintf("%s/%s_chr%s.rds", DIR_Est, ct, chr)
  
  if (!file.exists(out_file)) {
    P_Significant <- readRDS(sprintf("%s/P_Significant/%s_chr%s.rds", DIR_Test, ct, chr))
    Gene_SNP_list <- lapply(P_Significant, names)
    
    if (length(Gene_SNP_list) > 0) {
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
        # Gene expression, cis-SNP lists, null-model fits, and genotypes
        {
          Y_mat <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
          eta_hat0 <- readRDS(sprintf("%s/%s_chr%s.rds", DIR_MLE0, ct, chr))
          if (!identical(names(Gene_SNP_list), names(eta_hat0))) {
            eta_hat0 <- eta_hat0[names(Gene_SNP_list)]
            names(eta_hat0) <- names(Gene_SNP_list)
          }
          Y_mat <- Y_mat[names(Gene_SNP_list), , drop = F]
          G_SNP_mat_ <- G_SNP_mat[, INDI, drop = F]
        }
      }
      
      # Run effect-size estimation
      {
        cat(sprintf("[%s | chr%d] estimating effect curves for %d genes\n",
                    ct, chr, length(Gene_SNP_list)))
        res_Estimation <- scarfQTL(
          Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
          W_INDI_ = W_INDI_, W_CELL = W_CELL, TIME = TIME, 
          Test_type = "RT", Gene_SNP_list = Gene_SNP_list, 
          G_SNP_mat = G_SNP_mat_, GRM_INDI = GRM_INDI, 
          eta_hat0 = eta_hat0, max_iter = 0,
          Estimation = T, ridge_lamb = c(0,0.1,1,10)
        )
      }
      
      # Save estimation results
      {
        tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
        saveRDS(res_Estimation$Estimation, file = tmp_file)
        file.rename(tmp_file, out_file)
      }
    }
  }else{
    cat(sprintf("[%s | chr%d] no significant eQTLs for effect estimation\n", ct, chr))
    saveRDS(NULL, file = out_file)
  }
}
