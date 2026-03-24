# ================================================================
# This script summarizes simulation outputs from scarf-QTL analyses.
# Main tasks:
#   (i) summarize type I error across null simulations,
#   (ii) save the simulated dynamic effect patterns,
#   (iii) summarize power across causal simulations,
#   (iv) summarize effect-size estimation results across ridge penalties,
#   (v) summarize permutation-based null analyses,
#   (vi) summarize runtime benchmarks.
# ================================================================

setwd("~/sceQTL")
library(dplyr)
library(qvalue)
library(reshape2)
library(stringr)

# ===== 1. Define directories and simulation settings =====
{
  DIR_res <- "results/simulation/intermediate/t1e_power"
  DIR_permut <- "results/simulation/intermediate/permutation/"
  DIR_runtime <- "results/simulation/intermediate/runtime"
  DIR_summary <- "results/simulation/intermediate/summary"
  dir.create(DIR_summary, recursive = TRUE, showWarnings = FALSE)
  
  CT_simu <- c("CD4 NC", "B IN", "Plasma")
  names(CT_simu) <- str_replace(tolower(CT_simu), " ", "")
  
  models <- c("normal", "poisson", "nbinom")
  Gamma_patterns <- c("Static", "Linear", "TwoStage_onORoff", "TwoStage_switch", "Peak", 
                      "Cyclic_positive", "Cyclic_switch", "Cyclic_2period")
  {
    Time <- seq(0.05,0.95,0.01); names(Time) <- paste0("pt_", Time)
    K_Phi <- 10
    phi <- splines::bs(Time, knots = seq(0,1,length.out = K_Phi-3+1), 
                       Boundary.knots = c(0,1))
    phi <- phi[,-ncol(phi)]
  }
}

# ===== 2. Summarize type I error =====
{
  fns <- list.files(DIR_res, pattern = "null_")
  T1E_df <- do.call(rbind, lapply(fns, function(fn){
    res_P <- readRDS(sprintf("%s/%s", DIR_res, fn))
    res_P <- res_P[str_detect(names(res_P), "P_")]
    alphas <- c(0.05, 0.1^(2:4))
    
    Df <- do.call(rbind, lapply(names(res_P), function(m){
      pmat <- sapply(res_P[[m]], c)
      pmat[is.na(pmat)] <- 1
      T1E <- sapply(alphas, function(alph) mean(pmat < alph))
      T1E_df <- data.frame(value = T1E, alpha = alphas)
      T1E_df$method <- m
      fn_split <- str_split_fixed(str_split_fixed(fn, "_", 2)[1,2], "[-\\.]", 3)[1,]
      T1E_df$Cell.Type <- CT_simu[intersect(fn_split, names(CT_simu))]
      T1E_df$model <- intersect(models, fn_split)[1]
      T1E_df$N <- length(pmat)
      T1E_df$CR <- T1E_df$alpha + qnorm(0.99) * 
        sqrt(T1E_df$alpha*(1-T1E_df$alpha) / T1E_df$N)
      
      return(T1E_df)
    }))
    return(Df)
  }))
  saveRDS(T1E_df, file = sprintf("%s/T1E_df.rds", DIR_summary))
}

# ===== 3. Save simulated effect patterns =====
{
  Tim <- seq(0,1,0.01)
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
    integrate_gamma <- function(f, n_grid = 1000){
      x <- (seq(n_grid)-0.5)/(n_grid)
      y <- sapply(x, f)
      c(raw = mean(y), abs = mean(abs(y)), square = mean(y^2), cubic = mean(abs(y)^3))
    }
  }
  gammaSNP <- sapply(GammaSNP_raw, function(gammaSNP_raw){
    gammaSNP_scale <- integrate_gamma(gammaSNP_raw)[2]
    gammaSNP_c <- sapply(Tim, gammaSNP_raw) / gammaSNP_scale
  })
  rownames(gammaSNP) <- Tim
  gamma_df <- reshape2::melt(gammaSNP)
  colnames(gamma_df)[1:2] <- c("pseudotime", "gamma_pattern")
  saveRDS(gamma_df, file = sprintf("%s/Gamma_df.rds", DIR_summary))
}

# ===== 4. Summarize power =====
{
  fns <- list.files(DIR_res, pattern = "eSNP_")
  Power_df <- do.call(rbind, lapply(fns, function(fn){
    res_P <- readRDS(sprintf("%s/%s", DIR_res, fn))
    alphas <- c(0.05, 0.1^(2:4))
    Power <- sapply(res_P, function(P_list){
      ps <- unlist(P_list)
      ps[is.na(ps)] <- 1
      ret <- sapply(alphas, function(alph) mean(ps < alph))
      names(ret) <- alphas
      return(ret)
    })
    Power_df <- reshape2::melt(Power)
    colnames(Power_df)[1:2] <- c("alpha", "method")
    fn_split <- str_split_fixed(str_split_fixed(fn, "_", 2)[1,2], "[\\-\\.]", 5)[1,]
    Power_df$Cell.Type <- CT_simu[intersect(fn_split, names(CT_simu))]
    Power_df$gamma_pattern <- intersect(Gamma_patterns, fn_split)[1]
    Power_df$effect_size <- as.numeric(intersect(as.character(seq(100)), fn_split)[1])
    Power_df$model <- intersect(models, fn_split)[1]
    return(Power_df)
  }))
  saveRDS(Power_df, file = sprintf("%s/Power_df.rds", DIR_summary))
}

# ===== 5. Summarize effect-size estimation =====
{
  fns <- list.files(DIR_res, pattern = "Estimation_")
  Est_df <- do.call(rbind, lapply(fns, function(fn){
    res_Est <- readRDS(sprintf("%s/%s", DIR_res, fn))
    ridge_lamb_seq <- c(0, 0.1, 1, 10, 100)
    Est_summarise <- sapply(ridge_lamb_seq, function(ridge_lamb){
      Est_mat <- sapply(res_Est$Estimation, function(Estimation_G){
        beta <- Estimation_G[[as.character(ridge_lamb)]][,1]
        if (is.null(beta)) beta <- rep(NA, 10)
        return(beta)
      })
      Gamma_mat <- t(Est_mat) %*% t(phi)
      Gamma_mat_scale <- rowSums(abs(Gamma_mat))
      Gamma_mat_mean <- colMeans(Gamma_mat, na.rm = T)
      Gamma_mat_sd <- apply(Gamma_mat, 2, sd)
      c(Gamma_mat_mean, Gamma_mat_sd, 
        Gamma_mat_mean/mean(abs(Gamma_mat_mean)), Gamma_mat_sd/mean(abs(Gamma_mat_mean)))
    })
    
    Mean_effect <- Est_summarise[seq(length(Time)),]
    SD_effect <- Est_summarise[length(Time)+seq(length(Time)),]
    Mean_effect_scale <- Est_summarise[2*length(Time)+seq(length(Time)),]
    SD_effect_scale <- Est_summarise[3*length(Time)+seq(length(Time)),]
    rownames(Mean_effect) <- rownames(SD_effect) <- 
      rownames(Mean_effect_scale) <- rownames(SD_effect_scale) <- names(Time)
    colnames(Mean_effect) <- colnames(SD_effect) <- 
      colnames(Mean_effect_scale) <- colnames(SD_effect_scale) <- ridge_lamb_seq
    
    est_df <- reshape2::melt(Mean_effect)
    colnames(est_df) <- c("coef", "ridge_lamb", "mean_est")
    sd_df <- reshape2::melt(SD_effect)
    est_df$sd_est <- sd_df$value
    est_df2 <- reshape2::melt(Mean_effect_scale)
    est_df$mean_est_scale <- est_df2$value
    sd_df2 <- reshape2::melt(SD_effect_scale)
    est_df$sd_est_scale <- sd_df2$value
    est_df$ridge_lamb <- as.numeric(est_df$ridge_lamb)
    
    fn_split <- str_split_fixed(str_split_fixed(fn, "_", 2)[1,2], "[\\-\\.]", 5)[1,]
    est_df$Cell.Type <- CT_simu[intersect(fn_split, names(CT_simu))]
    est_df$gamma_pattern <- intersect(Gamma_patterns, fn_split)[1]
    est_df$effect_size <- as.numeric(intersect(as.character(seq(100)), fn_split)[1])
    est_df$model <- intersect(models, fn_split)[1]
    return(est_df)
  }))
  saveRDS(Est_df, file = sprintf("%s/Estimation_df.rds", DIR_summary))
}

# ===== 6. Summarize permutation-based null analyses =====
{
  Permute_res <- lapply(CT_simu, function(CT){
    ct <- str_replace(tolower(CT), " ", "")
    parallel::mclapply(seq(100), function(seed){
      fn <- sprintf("%s/permutation_%s_seed%s.rds", DIR_permut, ct, seed)
      P_df <- readRDS(fn)
      P_df_mat <- P_df[str_detect(colnames(P_df), "P_")]
      P_df_mat[is.na(P_df_mat)] <- 1
      
      FDR <- t(apply(P_df_mat, 2, function(ps){
        qobj <- qvalue(ps, lambda = 0.5)
        N_eSNP_g <- tapply(qobj$lfdr<0.05, P_df$GENE_ID, sum)
        c(N_eQTL = sum(N_eSNP_g), N_eGene = sum(N_eSNP_g>0), pi1=1-qobj$pi0)
      }))
      Summary_df <- reshape2::melt(FDR) %>% 
        rename(c("method" = "Var1", "parameter" = "Var2")) %>% 
        mutate(Cell.Type = CT, seed = seed)
      return(Summary_df)
    }, mc.cores = 50)
  })
  Permute_df <- do.call(rbind, unlist(Permute_res, recursive = F))
  saveRDS(Permute_df, file = sprintf("%s/Permute_df.rds", DIR_summary))
}

# ===== 7. Summarize runtime benchmarks =====
{
  # Collect scarf-QTL runtime profiles
  {
    fns_RT <- list.files(DIR_runtime, pattern = "time_RT")
    Time_RT_list <- lapply(fns_RT, function(fn){
      readRDS(sprintf("%s/%s", DIR_runtime, fn))
    })
    names(Time_RT_list) <- fns_RT
    fn_RT_df <- apply(str_split_fixed(fns_RT, "[_.]", 7)[, 3:6], 2, as.numeric)
    colnames(fn_RT_df) <- c("N_indi", "N_cell", "N_gene", "N_snp")
    fn_RT_df <- data.frame(fn_RT_df, fn = fns_RT)
    
    Time_RT_df <- do.call(rbind, lapply(fns_RT, function(fn){
      tl <- Time_RT_list[[fn]]
      df <- as.data.frame(tl$step_times)
      df$total <- sum(unlist(tl$step_times))
      df$fn <- fn
      return(df)
    }))
    Time_RT_long_df <- reshape2::melt(Time_RT_df) %>%
      rename(step = variable, seconds = value) %>%
      filter(step %in% c("eigendecomp", "null", "retro", "total")) %>%
      left_join(fn_RT_df, by = "fn") %>%
      mutate(Method = "scarf-QTL")
  }
  
  # Collect pseudobulk runtime profiles and merge summaries
  {
    fns_PB <- list.files(DIR_runtime, pattern = "time_PB")
    Time_PB_list <- lapply(fns_PB, function(fn){
      readRDS(sprintf("%s/%s", DIR_runtime, fn))
    })
    names(Time_PB_list) <- fns_PB
    fn_PB_df <- apply(str_split_fixed(fns_PB, "[_.]", 7)[, 3:6], 2, as.numeric)
    colnames(fn_PB_df) <- c("N_indi", "N_cell", "N_gene", "N_snp")
    fn_PB_df <- data.frame(fn_PB_df, fn = fns_PB)
    
    Time_PB_df <- do.call(rbind, lapply(fns_PB, function(fn){
      tl <- Time_PB_list[[fn]]
      df <- as.data.frame(tl$step_times)
      df$fn <- fn
      return(df)
    }))
    Time_PB_long_df <- reshape2::melt(Time_PB_df) %>%
      rename(step = variable, seconds = value) %>%
      filter(step %in% c("test")) %>%
      left_join(fn_PB_df, by = "fn") %>%
      mutate(Method = "Pseudobulk")
  }
  
  Time_long_df <- rbind(Time_RT_long_df, Time_PB_long_df) %>%
    mutate(Method = factor(Method, levels = c("scarf-QTL", "Pseudobulk")),
           Step = factor(step, levels = c("eigendecomp", "null", "retro", "total", "test"),
                         labels = c("Eigen transformation", "Null model fitting", 
                                  "Retrospective association test", "Total", 
                                  "Spearman correlation test")))
  
  saveRDS(Time_long_df, file = sprintf("%s/Runtime_df.rds", DIR_summary))
}