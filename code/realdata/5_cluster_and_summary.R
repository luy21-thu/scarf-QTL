# ================================================================
# This script summarizes effect-estimation and testing outputs across
# cell types and chromosomes, and prepares downstream clustering and
# visualization tables.
# Main tasks:
#   (i) collect estimated effect-size coefficients across chromosomes,
#   (ii) cluster dynamic eQTL trajectories within each cell type,
#   (iii) collect significant eQTL testing results and top SNP summaries,
#   (iv) merge summary tables with cell-type and rsID annotations,
#   (v) prepare QQ-plot input tables for each cell type.
# ================================================================

setwd("~/sceQTL")
library(dplyr)
library(stringr)

{
  DIR_res <- "results/realdata/intermediate"
  DIR_summary <- "results/realdata/summary"
  dir.create(DIR_summary, recursive = TRUE, showWarnings = FALSE)
  corder <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", 
              "NK", "NK R", "Plasma","B Mem", "B IN",
              "Mono C", "Mono NC","DC")
}

# ===== 1. Collect estimated effect sizes =====
{
  eQTL_est <- do.call(rbind, lapply(corder, function(CT){
    do.call(rbind, lapply(seq(22), function(chr){
      cat(CT, chr, "\n")
      ct <- str_replace(tolower(CT), " ", "")
      res_Est <- NULL
      try(res_Est <- readRDS(sprintf("%s/Estimation/%s_chr%s.rds", DIR_res, ct, chr)))
      res_Estimation <- lapply(res_Est, "[[", "1")
      # res_Estimation <- res_Est
      if(length(res_Estimation)>0){
        Beta <- as.matrix(do.call(rbind, lapply(res_Estimation, t)))
        colnames(Beta) <- paste0("beta_", seq(10))
        Beta_df <- data.frame(
          CELL_ID = ct,
          GENE_ID = rep(names(res_Estimation), sapply(res_Estimation, ncol)),
          SNPID = unlist(lapply(res_Estimation, colnames)),
          Beta
        )
      }else{
        Beta_df <- NULL
      }
      return(Beta_df)
    }))
  }))
  saveRDS(eQTL_est, file = sprintf("%s/eQTL_est.rds", DIR_summary))
}

# ===== 2. Cluster effect-size trajectories =====
{
  # build eQTL_Cluster with empty column "Cluster"
  cellid <- data.frame(Cell.type = corder, CELL_ID = str_replace(tolower(corder), " ", ""))
  eQTL_Cluster <- left_join(eQTL_est, cellid)
  rownames(eQTL_Cluster) <- paste0(eQTL_Cluster$Cell.type, "|", 
                                   eQTL_Cluster$GENE_ID, "|", 
                                   eQTL_Cluster$SNPID)
  eQTL_Cluster$Cluster <- NA
  
  # build spline basis
  {
    Time <- seq(0.05,0.95,0.01); names(Time) <- Time
    K_Phi <- 10
    phi <- splines::bs(Time, knots = seq(0,1,length.out = K_Phi-3+1), 
                       Boundary.knots = c(0,1))
    phi <- phi[,-ncol(phi)]
  }
  
  # run kmeans per cell type
  for(CT in corder){
    eQTL_df <- eQTL_Cluster[(eQTL_Cluster$Cell.type == CT) & (!is.na(eQTL_Cluster$beta_1)),]
    coef <- as.matrix(eQTL_df[, sprintf("beta_%s", seq(K_Phi))])
    fit <- coef %*% t(phi)
    fit_normalized <- fit/rowMeans(abs(fit))
    if(F){
      set.seed(1024)
      maxclunum <- 8 #20
      rss <- lapply(seq_len(maxclunum), function(clunum) {
        cat("=")
        tmp <- kmeans(fit_normalized, clunum, iter.max = 1000)
        tmp$betweenss / tmp$totss
      })
      rss <- unlist(rss)
      x <- 2:maxclunum
      number.cluster <-
        x[which.min(sapply(seq_len(length(x)), function(i) {
          x2 <- pmax(0, x - i)
          sum(lm(rss[-1] ~ x + x2)$residuals ^ 2)  ## check this
        }))]
      cat(CT, number.cluster, "\n")
    }else{
      # we recorded the chosen K here
      Nclu <- c(6,5,6, 5,5,7, 5,4, 6,6,6, 5,5,5)
      names(Nclu) <- corder
      number.cluster <- Nclu[CT]
      cat(CT,"")
    }
    {
      set.seed(1024)
      Clu <- kmeans(fit_normalized, number.cluster)
      clu <- Clu$cluster
      eQTL_Cluster[rownames(coef), "Cluster"] <- paste0(CT, "_clu", clu)
    }
  }
  
  # save clustered table
  saveRDS(eQTL_Cluster, file = sprintf("%s/eQTL_Cluster.rds", DIR_summary))
}

# ===== 3. Collect eQTL testing summaries =====
{
  eQTL_test <- do.call(rbind, lapply(corder, function(CT){
    do.call(rbind, lapply(seq(22), function(chr){
      cat(CT, chr, "\n")
      ct <- str_replace(tolower(CT), " ", "")
      P_df <- readRDS(sprintf("%s/Test/eQTL/eSNP_%s_chr%s.rds", DIR_res, ct, chr))
      Sig <- rowSums(as.matrix(P_df[ , str_detect(colnames(P_df), "FDR")])<0.05)
      P_df <- P_df[Sig>0,]
      return(P_df)
    }))
  }))
  saveRDS(eQTL_test, file = sprintf("%s/eQTL_test.rds", DIR_summary))
  
  eQTL_top <- do.call(rbind, lapply(corder, function(CT){
    df_top <- do.call(rbind, lapply(seq(22), function(chr){
      cat(CT, chr, "\n")
      ct <- str_replace(tolower(CT), " ", "")
      readRDS(sprintf("%s/Test/eQTL/topSNP_%s_chr%s.rds", DIR_res, ct, chr))
    }))
    return(df_top)
  }))
  saveRDS(eQTL_top, file = sprintf("%s/eQTL_top.rds", DIR_summary))
}

# ===== 4. Build merged summary table =====
{
  cTeX <- sprintf("$%s%s$", str_replace(corder, " ", "_{"), ifelse(str_detect(corder, " "), "}", ""))
  tol14rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
                    "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", 
                    "#F1932D", "#E8601C", "#DC050C")
  cellplot <- data.frame(CELL_ID = str_replace(tolower(corder), " ", ""), 
                         Cell.type = corder, cTeX = cTeX, color = tol14rainbow)
  eQTL_summary <- left_join(full_join(eQTL_test, select(eQTL_Cluster, -"Cell.type")), cellplot)
  
  rsid <- readRDS("data/genotype/rsid.rds")
  eQTL_summary <- left_join(eQTL_summary, rsid)
  saveRDS(eQTL_summary, file = sprintf("%s/eQTL_summary.rds", DIR_summary))
}

# ===== 5. Prepare QQ-plot data =====
{
  n_quant <- 1e5
  parallel::mclapply(corder, function(CT){
    ct <- stringr::str_replace(tolower(CT), " ", "")
    
    # load pseudobulk reference p-values
    eQTL_Yazar_CT <- readr::read_tsv(sprintf("data/eQTL_Yazar/%s_eqtl_table.tsv.gz", ct))
    eQTL_Yazar_CT_ <- (eQTL_Yazar_CT %>% arrange(P_VALUE) %>% 
                         mutate(POS = as.numeric(str_split_fixed(SNPID, "[:_]", 3)[,2]),
                                eQTL = paste(GENE, SNPID)))[,c("eQTL", "P_VALUE")]
    
    # load scarf-QTL p-values
    P_RT <- do.call(rbind, parallel::mclapply(seq(22), function(chr){
      P_df <- readRDS(sprintf("%s/Test/P_df/%s_chr%s.rds", DIR_res, ct, chr))
      P_df$eQTL <- paste(P_df$GENE_ID, P_df$SNPID)
      P_df <- P_df[,c("eQTL", "P_static", "P_dynamic", "P_combined")]
    }, mc.cores = 22))
    P_df <- left_join(eQTL_Yazar_CT_, P_RT, by = c("eQTL"))
    P_df <- P_df[,c("eQTL", "P_VALUE", "P_static", "P_dynamic", "P_combined")]
    colnames(P_df) <- c("eQTL", "Pseudobulk", 
                        "scar-QTL: static", "scar-QTL: dynamic", "scar-QTL: combined")

    # construct QQ sampled table
    qq_dfs <- do.call(rbind, lapply(colnames(P_df)[-1], function(m){
      p <- as.vector(P_df[[m]]); p <- p[!is.na(p)]
      p_sorted <- sort(p)
      
      max_logP_exp <- -log10((1 - 0.5) / length(p_sorted))
      min_logP_exp <- -log10((length(p_sorted) - 0.5) / length(p_sorted))
      seq_logP_exp <- seq(max_logP_exp, min_logP_exp, length.out = n_quant)
      idx <- round(10^(-seq_logP_exp) * length(p_sorted)+0.5)
      idx <- sort(unique(pmax(pmin(idx, length(p_sorted)), 1)))
      p_sampled <- p_sorted[idx]
      
      exp_p <- (idx - 0.5) / length(p_sorted)
      qq_df <- data.frame(
        Cell.Type = CT, 
        method = m,
        observed = -log10(p_sampled),
        expected = -log10(exp_p)
      )
      return(qq_df)
    }))
    
    # save QQ table
    saveRDS(qq_dfs, file = sprintf("%s/qqplot_df/qqplot_sample1e5_%s.rds", DIR_summary, ct))
  }, mc.cores = 14)
}