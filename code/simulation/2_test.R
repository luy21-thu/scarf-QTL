# ================================================================
# This script runs simulation-based association testing for scarf-QTL.
# For a selected simulation cell type (CT_idx), this script:
#   (i) evaluates null datasets under normal / Poisson / NB models,
#   (ii) evaluates causal datasets under multiple dynamic effect patterns,
#   (iii) saves p-value outputs and effect-size estimation results.
# Usage: Rscript 2_test.R <CT_idx>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")

library(data.table)
library(parallel)
library(stringr)

# ===== 1. Parse arguments and define simulation settings =====
{
  DIR <- "results/simulation/intermediate/simulated_data"
  DIR_res <- "results/simulation/intermediate/t1e_power"
  dir.create(DIR_res, recursive = TRUE, showWarnings = FALSE)
  
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
  
  models <- c("normal", "poisson", "nbinom")
  Gamma_patterns <- c("Static", "Linear", "TwoStage_onORoff", "TwoStage_switch", "Peak", 
                      "Cyclic_positive", "Cyclic_switch", "Cyclic_2period")
}

# ===== 2. Run null simulations =====
{
  load(sprintf("%s/sd_noY_%s.RData", DIR, ct))
  
  # Load or generate background SNP panels for null testing
  {
    if (FALSE) {
      G_SNP_mat_chr1 <- readRDS("data/genotype/chr1_genotype_matrix_rename.rds")
      SNP_df <- readRDS("data/genotype/rsid_AF.rds")
      #summary(SNP_df$A2_FREQ_ONEK1K)
      SNP_common <- intersect(SNP_df$SNPID, rownames(G_SNP_mat_chr1))
      set.seed(1)
      SNPs <- sample(SNP_common, 1e3)
      Gene_SNP_list_seed1 <- setNames(rep(list(SNPs), 1000), paste0("G", seq(1000)))
      G_SNP_mat_seed1 <- G_SNP_mat_chr1[SNPs, , drop = F]
      saveRDS(G_SNP_mat_seed1, file = "data/genotype/G_SNP_mat_seed1.rds")
    }else{
      G_SNP_mat_seed1 <- readRDS("data/genotype/G_SNP_mat_seed1.rds")
    }
    Gene_SNP_list_seed1 <- setNames(rep(list(rownames(G_SNP_mat_seed1)), 1000), 
                                    paste0("G", seq(1000)))
    G_SNP_mat_ <- G_SNP_mat_seed1[, INDI, drop = F]
  }
  
  # Run retrospective and prospective tests on null datasets
  Res_model <- mclapply(models, function(model){
    fn <- sprintf("%s-%s.rds", ct, model)
    if (!(sprintf("null_%s", fn) %in% list.files(DIR_res))) {
      Y_mat <- NA
      try(Y_mat <- readRDS(sprintf("%s/Ymat0-%s", DIR, fn)))
      if(!is.na(Y_mat)[1]){
        cat(sprintf("[%s] null simulation: %s\n", ct, fn))
        res_Score <- scarfQTL(
          Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
          W_INDI_ = W_INDI_, W_CELL = W_CELL, TIME = TIME, 
          Test_type = c("RT", "Standard"), 
          Gene_SNP_list = Gene_SNP_list_seed1, G_SNP_mat = G_SNP_mat_, 
          GRM_INDI = GRM_INDI, 
          Y_RINT = ifelse(model == "normal", FALSE, TRUE), eta_hat0 = NULL, 
          max_iter = 500, LM_cholesky = T, print_step = F, Estimation = F)
        saveRDS(res_Score, file = sprintf("%s/null_%s", DIR_res, fn))
      }
    }
  }, mc.cores = 3)
}

# ===== 3. Run causal simulations across dynamic patterns =====
{
  load(sprintf("%s/sd_noY_%s.RData", DIR, ct))
  Res_pattern <- mclapply(Gamma_patterns, function(gamma_pattern){
    # Define effect-size grid
    if (ct != "cd4nc") SNP_coef_list <- c(seq(5),seq(5)*10) else SNP_coef_list <- seq(5)
    Res_effectsize <- mclapply(SNP_coef_list, function(SNP_coef){
      Res_model <- lapply(models, function(model){
        # Run simulations across effect sizes and observation models
        fn <- sprintf("%s-%s-%s-%s", ct, SNP_coef, gamma_pattern, model)
        if (!(sprintf("eSNP_%s.rds", fn) %in% list.files(DIR_res))) {
          Y_mat <- NA
          try(Y_mat <- readRDS(sprintf("%s/Ymat-%s.rds", DIR, fn)))
          if (!is.na(Y_mat)[1]) {
            cat(sprintf("[%s] testing causal simulation: %s\n", ct, fn))
            
            ridge_lamb_seq <- c(0,0.1,1,10,100)
            res_Score <- scarfQTL(
              Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
              W_INDI_ = W_INDI_, W_CELL = W_CELL, 
              TIME = TIME, Test_type = c("RT"), 
              Gene_SNP_list = Gene_eSNP_list, G_SNP_mat = G_eSNP_mat, 
              GRM_INDI = GRM_INDI, 
              Y_RINT = ifelse(model == "normal", FALSE, TRUE), eta_hat0 = NULL, 
              max_iter = 500, LM_cholesky  =  T, print_step = F, 
              Estimation = T, ridge_lamb = ridge_lamb_seq)
            
            res_Spearman <- PB_runtime(
              Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
              W_INDI_ = W_INDI_, 
              Gene_SNP_list = Gene_eSNP_list, G_SNP_mat = G_eSNP_mat,
              Log2p = (model %in% c("poisson", "nbinom")), profile_steps = F)
            res_Score$P_Spearman <- res_Spearman
            
            saveRDS(res_Score[str_detect(names(res_Score), "P_")], 
                    file = sprintf("%s/eSNP_%s.rds", DIR_res, fn))
            saveRDS(res_Score["Estimation"], 
                    file = sprintf("%s/Estimation_%s.rds", DIR_res, fn))
          }
        }
        return(NULL)
      })
    }, mc.cores = 5)
  }, mc.cores = length(Gamma_patterns))
}


