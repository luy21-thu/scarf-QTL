# ================================================================
# This script benchmarks runtime of scarf-QTL and the pseudobulk baseline.
# Given input sample sizes (N_indi, N_cell, N_gene, N_snp), this script:
#   (i) resamples individuals, cells, genes, and SNPs from OneK1K-derived inputs,
#   (ii) runs scarf-QTL on the resampled dataset with profiling enabled,
#   (iii) runs pseudobulk-based testing with profiling enabled,
#   (iv) saves runtime profiles for downstream summaries.
# Usage: Rscript 4_runtime.R <N_indi> <N_cell> <N_gene> <N_snp>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")

library(data.table)
library(stringr)

# ===== 1. Parse arguments and define runtime settings =====
{
  DIR <- "results/simulation/intermediate/simulated_data"
  DIR_runtime <- "results/simulation/intermediate/runtime"
  dir.create(DIR_runtime, recursive = TRUE, showWarnings = FALSE)

  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 4)
  N_indi <- as.numeric(args[1])
  N_cell <- as.numeric(args[2])
  N_gene <- as.numeric(args[3])
  N_snp <- as.numeric(args[4])
  stopifnot(!is.na(N_indi), !is.na(N_cell), !is.na(N_gene), !is.na(N_snp))
  #N_indi <- 100; N_cell <- 500; N_gene <- 100; N_snp <- 3000
  
  CT <- "CD4 NC"; ct <- "cd4nc"; chr <- 1
}

# ===== 2. Load and resample inputs =====
{
  # Load base covariates, expression, and genotype data
  {
    cell_meta <- readRDS(sprintf("data/covariate/cell_pt_%s.rds", ct))
    individual_meta <- readRDS(sprintf("data/covariate/individual_%s.rds", ct))
    ids <- fread("data/genotype/GRM.rel.id", header = T)
    colnames(ids) <- c("FID", "IID")
    GRM <- as.matrix(fread("data/genotype/GRM.rel"))
    colnames(GRM) <- rownames(GRM) <- ids$IID
    
    Y_mat <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
    
    if (FALSE) {
      G_SNP_mat <- readRDS(
        sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))
      MAF <- readRDS("data/genotype/MAF_rename.rds")
      idx <- which(abs(MAF[[chr]]-0.5)<0.45,)
      G_SNP_mat_short <- G_SNP_mat[idx[seq(1e5)], , drop = F]
      saveRDS(G_SNP_mat_short, file = sprintf("%s/G_SNP_mat_short.rds", DIR_runtime))
    }else{
      G_SNP_mat <- readRDS(sprintf("%s/G_SNP_mat_short.rds", DIR_runtime))
    }
  }
  
  #resample individuals and cells
  {
    INDI_C <- cell_meta$individual; INDI <- sort(unique(INDI_C))
    
    INDI_sample <- sort(sample(INDI, N_indi, replace = (N_indi > length(INDI))))
    CELL_id <- rownames(cell_meta)
    CELL_sample_l <- lapply(INDI_sample, function(indi){
      CELL_i <- CELL_id[INDI_C == indi]
      CELL_sample_i <- sample(CELL_i, N_cell, replace = (N_cell > length(CELL_i)))
    })
    CELL_id_sample <- unlist(CELL_sample_l)
    INDI_C_sample <- rep(INDI_sample, each = N_cell)
  }
  
  # Build covariates and relatedness matrices for resampled data
  {
    W_INDI_ <- t(scale(individual_meta[INDI_sample, c(
      "sex", "age", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pf1", "pf2")]))
    W_CELL <- t(scale(cell_meta[CELL_id_sample, c("log_nUMI", "percent.mt")]))
    TIME <- cell_meta[CELL_id_sample, ]$pseudotime
  
    GRM_INDI <- GRM[INDI_sample, INDI_sample]
    std_devs <- sqrt(diag(GRM_INDI))
    GRM_INDI <- GRM_INDI / (std_devs %*% t(std_devs))
    G_SNP_mat_ <- G_SNP_mat[, INDI_sample, drop = F]
  }
  # Resample genes and SNPs
  {
    G_sample <- sort(sample(rownames(Y_mat), N_gene, replace = (N_gene>nrow(Y_mat))))
    Gene_SNP_list_sample <- lapply(G_sample, function(G){
      SNPs <- rownames(G_SNP_mat)
      SNPs_sample <- sort(sample(SNPs, N_snp, replace = (N_snp>length(SNPs))))
    })
    names(Gene_SNP_list_sample) <- G_sample
    Y_mat_sample <- Y_mat[G_sample, CELL_id_sample, drop = F]
  }
  # Finalize unique IDs and align names
  {
    INDI_sample_ <- make.unique(INDI_sample, ".")
    colnames(W_INDI_) <- colnames(GRM_INDI) <- rownames(GRM_INDI) <- 
      colnames(G_SNP_mat_) <- INDI_sample_
    INDI_C_sample <- rep(INDI_sample_, each = N_cell)
    
    CELL_id_sample_ <- CELL_id_sample
    CELL_id_sample_ <- make.unique(CELL_id_sample, ".")
    colnames(W_CELL) <- colnames(Y_mat_sample) <- CELL_id_sample_
  }
}

# ===== 3. Run runtime benchmarks =====
{
  # Run retrospective and prospective tests on resampled datasets
  cat(sprintf("[runtime] N_indi=%d, N_cell=%d, N_gene=%d, N_snp=%d\n",
              N_indi, N_cell, N_gene, N_snp))
  res_Score <- scarfQTL(
    Y_mat = Y_mat_sample, INDI = INDI_sample_, INDI_C = INDI_C_sample, 
    W_INDI_ = W_INDI_, W_CELL = W_CELL, TIME = TIME, 
    Test_type = "RT", 
    Gene_SNP_list = Gene_SNP_list_sample, G_SNP_mat = G_SNP_mat_, 
    GRM_INDI = GRM_INDI, 
    Y_RINT = T, max_iter = 5000, eta_hat0 = NULL, LM_cholesky = T, 
    print_step = F, profile_steps = T)
  prof_Score <- attr(res_Score, "profile")
  fn_Score <- sprintf("%s/time_RT_%s.rds", DIR_runtime, paste0(args, collapse = "_"))
  saveRDS(prof_Score, file = fn_Score)
  
  res_Spearman <- PB_runtime(
    Y_mat = Y_mat_sample, INDI = INDI_sample_, INDI_C = INDI_C_sample, 
    W_INDI_ = W_INDI_, 
    Gene_SNP_list = Gene_SNP_list_sample, G_SNP_mat = G_SNP_mat_,
    Log2p = T, profile_steps = T)
  prof_Spearman <- attr(res_Spearman, "profile")
  fn_Spearman <- sprintf("%s/time_PB_%s.rds", DIR_runtime, paste0(args, collapse = "_"))
  saveRDS(prof_Spearman, file = fn_Spearman)
}
