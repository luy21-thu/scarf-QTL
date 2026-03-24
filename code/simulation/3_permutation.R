# ================================================================
# This script runs permutation-based null analyses for scarf-QTL.
# For a selected simulation cell type (CT_idx), this script:
#   (i) loads chromosome-1 real-data inputs and fitted null-model parameters,
#   (ii) permutes individuals in genotype and GRM matrices,
#   (iii) runs retrospective, prospective, and pseudobulk-based tests,
#   (iv) saves per-permutation p-value tables for downstream summaries.
# Usage: Rscript 3_permutation.R <CT_idx>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")

library(data.table)
library(parallel)
library(stringr)

# ===== 1. Parse arguments and define simulation (permutation) settings =====
{
  DIR <- "results/simulation/intermediate/simulated_data"
  DIR_permut <- "results/simulation/intermediate/permutation"
  dir.create(DIR_permut, recursive = TRUE, showWarnings = FALSE)
  
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1)
  CT_idx <- as.integer(args[1])
  stopifnot(!is.na(CT_idx), CT_idx >= 1, CT_idx <= 3)
  
  #Example
  #CT_idx <- 3
  CT_simu <- c("CD4 NC", "B IN", "Plasma")
  names(CT_simu) <- str_replace(tolower(CT_simu), " ", "")
  CT <- CT_simu[CT_idx]
  ct <- str_replace(tolower(CT), " ", "")
}

# ===== 2. Load data and null-model inputs =====
{
  chr <- 1
  load(sprintf("%s/sd_noY_%s.RData", DIR, ct))
  Gene_SNP_list <- readRDS(sprintf("data/Gene_cis_SNP/%s_chr%s.rds", ct, chr))
  Y_mat <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
  G_SNP_mat <- readRDS(sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))
}

# ===== 3. Run permutation-based null analyses =====
parallel::mclapply(seq(100), function(seed){
  # Generate permuted genotype / GRM inputs
  {
    set.seed(seed)
    INDI_permute <- sample(INDI)
    GRM_INDI_permute <- GRM_INDI[INDI_permute, INDI_permute]
    colnames(GRM_INDI_permute) <- rownames(GRM_INDI_permute) <- INDI
    G_SNP_mat_permute <- G_SNP_mat[, INDI_permute, drop = F]
    colnames(G_SNP_mat_permute) <- INDI
  }
  
  # Align null-model parameters
  eta_hat0 <- readRDS(sprintf(
    "results/realdata/intermediate/MLE0/%s_chr%s.rds", ct, chr))
  if(!identical(names(Gene_SNP_list), names(eta_hat0))){
    eta_hat0 <- eta_hat0[names(Gene_SNP_list)]
    names(eta_hat0) <- names(Gene_SNP_list)
  }
  
  # Run association tests
  {
    cat(sprintf("[%s] permutation seed %d\n", ct, seed))
    res_Score <- scarfQTL(
      Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
      W_INDI_ = W_INDI_, W_CELL = W_CELL, TIME = TIME, 
      Test_type = c("RT", "Standard"), Gene_SNP_list = Gene_SNP_list, 
      G_SNP_mat = G_SNP_mat_permute, GRM_INDI = GRM_INDI_permute, 
      eta_hat0 = eta_hat0, max_iter = 0)
    
    res_Spearman <- PB_runtime(
      Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
      W_INDI_ = W_INDI_, 
      Gene_SNP_list = Gene_SNP_list, G_SNP_mat = G_SNP_mat_permute,
      Log2p = T, profile_steps = F)
    res_Score$P_Spearman <- res_Spearman
  }
  
  # Reformat and save permutation results
  {
    P <- res_Score$P_dynamic
    P_df <- cbind(data.frame(
      CELL_ID = ct,
      Chromosome = chr,
      GENE_ID = rep(names(P), sapply(P, length)),
      SNPID = (unlist(lapply(P, names)))),
      as.data.frame(lapply(res_Score[-1], unlist)))
    saveRDS(P_df, file = sprintf("%s/permutation_%s_seed%s.rds", DIR_permut, ct, seed))
  }
}, mc.cores = 50)

