# ================================================================
# This script tests eQTL effect (retrospective association test).
# For each chromosome (chr = 1..22) and each of 14 immune cell types, this script:
#   (i) loads pre-split SCT expression (genes × cells), covariates,
#       cis-SNP lists, per-gene null-fit eta_hat, and per-chr genotypes;
#   (ii) runs scarfQTL with Test_type = "RT" (retrospective score tests);
#   (iii) writes per cell type outputs under scarf_QTL/Test/, including:
#         - raw p-values (P_dynamic / P_static / P_combined) and qvalue/lfdr summaries,
#         - per-gene top SNPs, significant eSNP tables, and per-gene significant SNP lists.
# Usage: Rscript 3_retrospective_test.R <chr>
# ================================================================

setwd("~/sceQTL")
source("code/scarf_QTL_function.R")
library(data.table)
library(dplyr)
library(qvalue)
library(stringr)

# Parse command-line arguments and cell-type labels
{
  args <- commandArgs(trailingOnly = TRUE)
  chr <- as.integer(args[1])
  stopifnot(!is.na(chr), chr >= 1, chr <= 22)
  
  # Example
  # chr <- 1
  CT_14 <- c("B IN", "B Mem", "CD4 ET", "CD4 NC", "CD4 SOX4", 
             "CD8 ET", "CD8 NC", "CD8 S100B", "DC", 
             "Mono C", "Mono NC", "NK", "NK R", "Plasma")
  
  DIR_MLE0 <- "results/realdata/intermediate/MLE0"
  DIR_Test <- "results/realdata/intermediate/Test"
  dir.create(DIR_Test, recursive = TRUE, showWarnings = FALSE)
}



# Read genotype data
G_SNP_mat <- readRDS(sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))

for (CT_idx in seq(14)) {
  CT <- CT_14[CT_idx]
  ct <- str_replace(tolower(CT), " ", "")
  
  out_file <- sprintf("%s/%s_chr%s.rds", DIR_Test, ct, chr)
  
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
        Gene_SNP_list <- readRDS(sprintf(
          "data/Gene_cis_SNP/%s_chr%s.rds", ct, chr))
        eta_hat0 <- readRDS(sprintf("%s/%s_chr%s.rds", DIR_MLE0, ct, chr))
        if (!identical(names(Gene_SNP_list), names(eta_hat0))) {
          eta_hat0 <- eta_hat0[names(Gene_SNP_list)]
          names(eta_hat0) <- names(Gene_SNP_list)
        }
        G_SNP_mat_ <- G_SNP_mat[, INDI, drop = F]
      }
    }
    
    # Run retrospective association tests
    cat(sprintf("[%s | chr%d] running retrospective eQTL tests for %d genes\n",
                ct, chr, length(Gene_SNP_list)))
    res_Score <- scarfQTL(
      Y_mat = Y_mat, INDI = INDI, INDI_C = INDI_C, 
      W_INDI_ = W_INDI_, W_CELL = W_CELL, TIME = TIME, 
      Test_type = "RT", Gene_SNP_list = Gene_SNP_list, 
      G_SNP_mat = G_SNP_mat_, GRM_INDI = GRM_INDI, 
      eta_hat0 = eta_hat0, max_iter = 0)
    
    # Save results
    {
      tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
      saveRDS(res_Score, file = tmp_file)
      file.rename(tmp_file, out_file)
    }
  }else{
    res_Score <- readRDS(out_file)
  }
  
  # Summarize p-values, q-values, and significant eQTL results
  {
    if (length(unlist(res_Score$P_dynamic))>0) {
      P <- res_Score$P_dynamic
      P_df <- cbind(data.frame(
        CELL_ID = ct,
        Chromosome = chr,
        GENE_ID = rep(names(P), sapply(P, length)),
        SNPID = (unlist(lapply(P, names)))),
        as.data.frame(lapply(res_Score[str_detect(names(res_Score), "P_")], unlist)))
    }else{
      # No tested SNPs for this cell type / chromosome
      P_df <- NULL
    }
    
    sapply(sprintf("%s/%s", DIR_Test, c("eQTL", "P_df", "P_Significant")), dir.create,
           recursive = TRUE, showWarnings = FALSE)
    if (nrow(P_df)>0) {
      Qs <- lapply(P_df[, c("P_dynamic", "P_static", "P_combined")], function(P) {
        P[is.na(P)] <- 1
        qobj <- NA
        try(qobj <- qvalue(P, lambda = 0.5), F)
        if (is.na(qobj)[1]) {
          qobj <- list(qvalues = NA, lfdr = p.adjust(P, "BY"))
          cat(ct, "fail")
        }else{
          cat(qobj$pi0, "")
        }
        return(qobj)
      })
      P_df <- cbind(
        P_df,
        Q_dynamic = Qs$P_dynamic$qvalues, 
        Q_static = Qs$P_static$qvalues, 
        Q_combined = Qs$P_combined$qvalues, 
        FDR_dynamic = Qs$P_dynamic$lfdr,
        FDR_static = Qs$P_static$lfdr,
        FDR_combined = Qs$P_combined$lfdr
      )
      NNA <- sapply(P_df, function(i) sum(is.na(i)))
      if (sum(NNA)>0) cat(sprintf("NA_%s_chr%s.rds", ct, chr))
      saveRDS(P_df, file = sprintf("%s/P_df/%s_chr%s.rds", DIR_Test, ct, chr))
      
      df_top <- P_df %>% group_by(GENE_ID) %>% 
        summarise(topSNP_combined = SNPID[which.min(P_combined)], 
                  topP_combined = min(P_combined),
                  topSNP_static = SNPID[which.min(P_static)], 
                  topP_static = min(P_static)) %>% 
        mutate(Cell.type = CT)
      saveRDS(df_top, file = sprintf("%s/eQTL/topSNP_%s_chr%s.rds", DIR_Test, ct, chr))
      
      Sig <- rowSums(as.matrix(P_df[ , str_detect(colnames(P_df), "FDR")])<0.05)
      P_df_Significant <- P_df[Sig>0,]
      saveRDS(P_df_Significant, file = sprintf("%s/eQTL/eSNP_%s_chr%s.rds", DIR_Test, ct, chr))
      
      P_df_Significant <- P_df_Significant[P_df_Significant$FDR_combined<0.05,]
      P_Significant <- list()
      if (!is.null(P_df_Significant)) {
        if (nrow(P_df_Significant)>0) {
          P_Significant <- tapply(
            seq(nrow(P_df_Significant)), P_df_Significant$GENE_ID, function(eqtl) {
              p_df <- P_df_Significant[eqtl,]
              Ps <- p_df$P_combined
              names(Ps) <- p_df$SNPID
              Ps <- sort(Ps)
              return(Ps)
            }, simplify = F)
        }
      }
      saveRDS(P_Significant, file = sprintf("%s/P_Significant/%s_chr%s.rds", DIR_Test, ct, chr))
    } else {
      cat(sprintf("No eQTL tested in %s, chr%s\n", ct, chr))
      saveRDS(NULL, file = sprintf("%s/eQTL/topSNP_%s_chr%s.rds", DIR_Test, ct, chr))
      saveRDS(NULL, file = sprintf("%s/eQTL/eSNP_%s_chr%s.rds", DIR_Test, ct, chr))
      saveRDS(NULL, file = sprintf("%s/P_df/%s_chr%s.rds", DIR_Test, ct, chr))
      saveRDS(NULL, file = sprintf("%s/P_Significant/%s_chr%s.rds", DIR_Test, ct, chr))
    }
  }
}
