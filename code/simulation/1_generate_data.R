# ================================================================
# This script generates simulated datasets for scarf-QTL evaluations.
# For a selected simulation cell type (CT_idx), this script:
#   (i) prepares shared covariates, relatedness matrices, and causal genotypes,
#   (ii) generates null expression matrices under normal / Poisson / NB models,
#   (iii) generates causal expression matrices under multiple dynamic effect patterns.
# Usage: Rscript 1_generate_data.R <CT_idx>
# ================================================================

setwd("~/sceQTL")
library(data.table)
library(parallel)
library(splines)
library(stringr)

# ===== 1. Parse arguments and define simulation settings =====
{
  DIR <- "results/simulation/intermediate/simulated_data"
  dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
  
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
  
  set.seed(42)
  N_sim <- 1000
  models <- c("normal", "poisson", "nbinom")
}

# ===== 2. Generate shared simulation inputs (excluding Y_mat) =====
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
  # Genotype of causal SNPs
  {
    G_sim <- paste0("G", seq(N_sim))
    Gene_eSNP_list <- setNames(lapply(G_sim, function(G) paste0("eSNP_", G)), G_sim)
    G_eSNP_mat <- matrix(rbinom(length(INDI) * N_sim, 2, 0.2), nrow = N_sim)
    colnames(G_eSNP_mat) <- INDI
    rownames(G_eSNP_mat) <- paste0("eSNP_", names(Gene_eSNP_list))
  }
  # save these data
  save(Gene_eSNP_list, G_eSNP_mat, GRM_INDI, INDI, INDI_C, W_INDI_, W_CELL, TIME,
       file = sprintf("%s/sd_noY_%s.RData", DIR, ct))
}

# ===== 3. Generate simulated expression matrices =====
{
  load(sprintf("%s/sd_noY_%s.RData", DIR, ct))
  
  # Construct spline basis and shared components
  {
    TIME_ <- rank(TIME)/(length(TIME)+1)
    K_Phi <- 10
    Phi <- bs(TIME_, knots = seq(0,1,length.out = K_Phi-3+1), Boundary.knots = c(0,1))
    Phi <- Phi[, -ncol(Phi), drop = F]; colnames(Phi) <- sprintf("bs%s", seq(K_Phi))
    
    covariate <- rbind(W_CELL, W_INDI_[,INDI_C])
    tau <- 1; eta <- 1*tau
  }
  
  # Generate null expression means
  {
    Y_mat_mu_null_ <- mclapply(names(Gene_eSNP_list), function(G) {
      beta_covariate <- rnorm(nrow(covariate))
      Y_covariate <- as.vector(beta_covariate %*% covariate)
      beta_baseline <- seq(K_Phi)/K_Phi
      Y_baseline <- as.vector(Phi %*% beta_baseline)
      U_indi <- matrix(rnorm(length(INDI)*K_Phi), ncol = K_Phi)
      rownames(U_indi) <- INDI
      Y_indi <- sqrt(eta/tau) * sapply(seq(length(INDI_C)), function(c_idx) {
        sum(Phi[c_idx,] * U_indi[INDI_C[c_idx],])
      })
      Y <- Y_covariate + Y_baseline + Y_indi
      return(Y)
    }, mc.cores = 50)
    Y_mat_mu_null_t <- sapply(Y_mat_mu_null_, c)
    colnames(Y_mat_mu_null_t) <- names(Gene_eSNP_list)
    rownames(Y_mat_mu_null_t) <- INDI_C
    saveRDS(Y_mat_mu_null_t, file = sprintf("%s/Y_mat_mu_null-%s.rds", DIR, ct))
  }
  
  # Generate null datasets
  {
    Y_mat_mu <- t(Y_mat_mu_null_t)
    fn_model <- sprintf("Ymat0-%s-%s.rds", ct, models)
    listfiles <- list.files(DIR, pattern = "Ymat0")
    {
      if (!(fn_model[1] %in% listfiles)) {
        Y_mat_norm <- Y_mat_mu + sqrt(1/tau) * rnorm(length(Y_mat_mu))
        cat(sprintf("[%s] generating %s\n", ct, fn_model[1]))
        saveRDS(Y_mat_norm, file = sprintf("%s/%s", DIR, fn_model[1]))
      }
      
      Y_mat_mu_sd <- sd(Y_mat_mu)
      if (Y_mat_mu_sd == 0 || is.na(Y_mat_mu_sd)) Y_mat_mu_sd <- 1
      
      if (!(fn_model[2] %in% listfiles)) {
        Y_mat_poisson <- matrix(rpois(length(Y_mat_mu), exp(Y_mat_mu/Y_mat_mu_sd)), nrow = nrow(Y_mat_mu))
        dimnames(Y_mat_poisson) <- dimnames(Y_mat_mu)
        cat(sprintf("[%s] generating %s\n", ct, fn_model[2]))
        saveRDS(Y_mat_poisson, file = sprintf("%s/%s", DIR, fn_model[2]))
      }
      
      if (!(fn_model[3] %in% listfiles)) {
        Y_mat_nbinom <- matrix(rnbinom(length(Y_mat_mu), mu = exp(Y_mat_mu/Y_mat_mu_sd), size=1), 
                               nrow = nrow(Y_mat_mu))
        dimnames(Y_mat_nbinom) <- dimnames(Y_mat_mu)
        cat(sprintf("[%s] generating %s\n", ct, fn_model[3]))
        saveRDS(Y_mat_nbinom, file = sprintf("%s/%s", DIR, fn_model[3]))
      }
    }
  }
  
  # Generate causal datasets
  {
    GammaSNP_raw <- list(
      Static = function(t) 1,
      Linear = function(t) t,
      TwoStage_onORoff = function(t) 1/(exp(-(20*(t-0.4)))+1),
      TwoStage_switch = function(t) 1/(exp(-(20*(t-0.4)))+1)-0.5,
      Peak = function(t) exp(-(20*(t-0.9)^2)),
      Cyclic_positive = function(t) sin(2*pi*t)+1,
      Cyclic_switch = function(t) sin(2*pi*t),
      Cyclic_2period = function(t) sin(4*pi*t)+1
    )
    integrate_gamma <- function(f, n_grid = 1000) {
      x <- (seq(n_grid)-0.5)/(n_grid)
      y <- sapply(x, f)
      c(raw = mean(y), abs = mean(abs(y)))
    }
    
    Y_pattern <- mclapply(names(GammaSNP_raw), function(gamma_pattern) {
      gammaSNP_raw <- GammaSNP_raw[[gamma_pattern]]
      gammaSNP_scale <- integrate_gamma(gammaSNP_raw)[2]
      gammaSNP_c <- sapply(TIME_, gammaSNP_raw) / gammaSNP_scale
      Y_mat_SNP1_t <- t(G_eSNP_mat[,INDI_C]) * gammaSNP_c
      cat(sprintf("[%s] pattern %s: generating causal datasets\n", ct, gamma_pattern))
      
      if (ct != "cd4nc") SNP_coef_list <- c(seq(5),seq(5)*10) else SNP_coef_list <- seq(5)
      
      res_effectsize <- mclapply(SNP_coef_list, function(SNP_coef) {
        Y_mat_mu <- t(Y_mat_mu_null_t + Y_mat_SNP1_t * (SNP_coef*0.02))
        rownames(Y_mat_mu) <- names(Gene_eSNP_list)
        colnames(Y_mat_mu) <- colnames(W_CELL)
        
        fn_model <- sprintf("Ymat-%s-%s-%s-%s.rds", ct, SNP_coef, gamma_pattern, models)
        listfiles <- list.files(DIR, pattern = "Ymat-")
        {
          if (!(fn_model[1] %in% listfiles)) {
            Y_mat_norm <- Y_mat_mu + sqrt(1/tau) * rnorm(length(Y_mat_mu))
            cat(sprintf("simulated_pattern/%s\n", fn_model[1]))
            saveRDS(Y_mat_norm, file = sprintf("%s/%s", DIR, fn_model[1]))
          }
          
          Y_mat_mu_sd <- sd(Y_mat_mu)
          if (Y_mat_mu_sd == 0 || is.na(Y_mat_mu_sd)) Y_mat_mu_sd <- 1
          
          if (!(fn_model[2] %in% listfiles)) {
            Y_mat_poisson <- matrix(rpois(length(Y_mat_mu), exp(Y_mat_mu/Y_mat_mu_sd)), 
                                    nrow = nrow(Y_mat_mu))
            dimnames(Y_mat_poisson) <- dimnames(Y_mat_mu)
            cat(sprintf("simulated_pattern/%s\n", fn_model[2]))
            saveRDS(Y_mat_poisson, file = sprintf("%s/%s", DIR, fn_model[2]))
          }
          
          if (!(fn_model[3] %in% listfiles)) {
            Y_mat_nbinom <- matrix(rnbinom(length(Y_mat_mu), mu = exp(Y_mat_mu/Y_mat_mu_sd), size=1), 
                                   nrow = nrow(Y_mat_mu))
            dimnames(Y_mat_nbinom) <- dimnames(Y_mat_mu)
            cat(sprintf("simulated_pattern/%s\n", fn_model[3]))
            saveRDS(Y_mat_nbinom, file = sprintf("%s/%s", DIR, fn_model[3]))
          }
        }
      }, mc.cores = 5)
    }, mc.cores = length(GammaSNP_raw))
  }
}
