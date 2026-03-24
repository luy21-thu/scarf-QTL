# ================================================================
# This script generates main-text and supplementary figures for the
# real-data analyses, including eGene comparisons, QQ plots, example
# boxplots, heatmaps, and dendritic-cell subset visualizations.
# ================================================================

setwd("~/sceQTL")

library(cowplot)
library(dplyr)
library(grDevices)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(latex2exp)
library(pheatmap)
library(RColorBrewer)
library(slingshot)
library(stringr)
library(tibble)
library(tidyr)

DIR_fig <- "results/realdata/figures"
dir.create(DIR_fig, recursive = TRUE, showWarnings = FALSE)

# Load summary tables and reference results
{
  eQTL_summary <- readRDS("results/realdata/summary/eQTL_summary.rds")
  eQTL_Yazar <- read.csv("data/eQTL_Yazar/science.abf3041_tables_s10.csv", 
                         header = T, skip = 2)
  eQTL_summary$Yazar_eGene <- ifelse(paste(eQTL_summary$Cell.type, eQTL_summary$GENE_ID) %in% 
                                       paste(eQTL_Yazar$Cell.type, eQTL_Yazar$Gene.ID), T, F)
  rownames(eQTL_summary) <- paste0(eQTL_summary$Cell.type, "|", 
                                   eQTL_summary$GENE_ID, "|", eQTL_summary$SNPID)
  eQTL_topSNP <- readRDS("results/realdata/summary/eQTL_top.rds")
  corder <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", 
              "NK", "NK R", "Plasma","B Mem", "B IN",
              "Mono C", "Mono NC","DC")
}

# Define plotting palettes and spline basis
{
  tol14rainbow <- setNames(c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
                             "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", 
                             "#F1932D", "#E8601C", "#DC050C"), corder)
  method_colors_short <- c(
    static = "deepskyblue3",  
    combine = "mediumpurple3",
    Yazar = "orange"
  )
  method_colors <- c(
    "scar-QTL: static" = "deepskyblue3",  
    "scar-QTL: dynamic" = "brown3", 
    "scar-QTL: combined" = "mediumpurple3",
    "Pseudobulk" = "orange"
  )
  col_fun <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  col_fun2 <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100)
  {
    Time <- seq(0.05,0.95,0.01); names(Time) <- Time
    K_Phi <- 10
    phi <- splines::bs(Time, knots = seq(0,1,length.out = K_Phi-3+1), 
                       Boundary.knots = c(0,1))
    phi <- phi[,-ncol(phi)]
  }
}

# ===== Figure 3a. Venn diagrams of detected eGenes =====
{
  eGene <- lapply(corder, function(CT) {
    eGene_dynamic <- unique(eQTL_summary$GENE_ID[
      eQTL_summary$Cell.type == CT & eQTL_summary$FDR_combined<0.05])
    eGene_static <- unique(eQTL_summary$GENE_ID[
      eQTL_summary$Cell.type == CT & eQTL_summary$FDR_static<0.05])
    eGene_Yazar <- eQTL_Yazar$Gene.ID[
      eQTL_Yazar$Cell.type == CT & eQTL_Yazar$eSNP.rank == "eSNP1"]
    list(combine = eGene_dynamic, static = eGene_static, Yazar = eGene_Yazar)
  }); names(eGene) <- corder
  N_Gene <- table(eQTL_topSNP$Cell.type)[corder]
  
  gVenn <- lapply(seq(14), function(CT_idx) {
    cat(CT_idx)
    eGene_ct <- eGene[[CT_idx]][c(2,1,3)]
    
    category.names2 <- sprintf("%s(%s)", c(
      Yazar = "Yazar\n", static = "Static\n", combine = "scarf-QTL:combined "
    )[names(eGene_ct)], sapply(eGene_ct, length))
    g <- ggVennDiagram(eGene_ct, category.names = category.names2,
                       label = "count", set_size = 4, set_color = method_colors_short) +
      ggtitle(TeX(sprintf("$%s%s\\ (%s)$", str_replace(names(eGene)[CT_idx], " ", "_{"), 
                          ifelse(str_detect(names(eGene)[CT_idx], " "), "}", ""), 
                          N_Gene[CT_idx]))) +
      scale_x_continuous(expand = expansion(mult = 0.1))+
      scale_y_continuous(expand = expansion(mult = 0.1))+
      scale_fill_gradient(low = "white", high = tol14rainbow[CT_idx], guide = "none")+
      scale_color_manual(values = unname(method_colors_short))+
      theme(plot.title = element_text(hjust = 0.5))
  })
  pdf(sprintf("%s/Venn.pdf", DIR_fig), width = 18, height = 9)
  grid.arrange(grobs = gVenn, layout_matrix =
                 rbind(seq(6), c(seq(7,11), NA), c(seq(12,14), rep(NA,3))))
  dev.off()
}

# ===== Figure S4. Number of eGenes versus cell counts =====
{
  N_cell <- setNames(c(463496, 61777, 4065,  205008, 133470, 34524,
                       159762, 9674,  3625, 47648, 81284, 
                       38218, 15161, 8689), corder)
  ggplot(left_join(reshape2::melt(sapply(eGene, sapply, length)),
                   data.frame(Var2 = corder, N_cell = N_cell)), 
         aes(x = N_cell, y = value,
             shape = factor(Var1, levels = c("Yazar", "static", "combine"), 
                            labels = c("Pseudobulk", "scarf-QTL:static", "scarf-QTL:combined")), 
             col = Var2, group = Var2))+
    geom_line(linetype = "dotted") + geom_point()+
    scale_color_manual(values = tol14rainbow)+
    theme_classic()+
    scale_x_continuous(transform = "log10")+
    scale_shape(name = "Method")+
    labs(x = "Total number of cells", y = "number of eGenes", col = "Cell type")
  ggsave(sprintf("%s/N_egene.pdf", DIR_fig), width = 9, height = 6)
}

# ===== Figure S5. Overlap with MAGMA disease-associated genes =====
{
  MAGMA_Traits <- read.csv("data/MAGMA/scDRS_trait.csv")
  MAGMA_Traits_select <- filter(MAGMA_Traits, Category == "blood/immune", 
                                !str_detect(Trait_Identifier, "UKB_460K.blood"))
  MAGMA_P <- as.matrix(read.table("data/MAGMA/MAGMA_v108_GENE_10_PSTAT_for_scDRS.txt"))
  MAGMA_P <- MAGMA_P[, MAGMA_Traits_select$Trait_Identifier]
  
  MAGMA_long <- reshape2::melt(MAGMA_P)
  colnames(MAGMA_long) <- c("Gene_MAGMA", "Trait_Identifier", "p_MAGMA")
  MAGMA_df <- MAGMA_long %>% filter(!is.na(p_MAGMA)) %>% 
    group_by(Trait_Identifier) %>% 
    mutate(FDR_significant = (p.adjust(p_MAGMA, method = "BH"))<0.05,
           GENE = toupper(Gene_MAGMA)) %>% 
    filter(GENE %in% toupper(unique(eQTL_topSNP$GENE_ID)))
  #print(MAGMA_df %>% summarise(N_MAGMA = sum(FDR_significant)))
  
  Overlap_Methods <- setNames(lapply(names(eGene[[1]]), function(Method) {
    df_all <- MAGMA_df %>% 
      mutate(eGENE = GENE %in% toupper(unlist(lapply(eGene, "[[", Method)))) %>% 
      summarise(N_MAGMA = sum(FDR_significant), 
                N_eGENE = sum(eGENE),
                N_overlap = sum(eGENE & FDR_significant),
                OR = N_overlap/N_eGENE/N_MAGMA*n(),
                Fisher_p = fisher.test(
                  table(eGENE, FDR_significant), alternative = "greater")$p.value) %>% 
      mutate(Cell.Type = "All", method = Method)
    
    df_CT <- do.call(rbind, lapply(corder, function(CT) {
      MAGMA_df %>% mutate(eGENE = GENE %in% toupper(unlist(eGene[[CT]][[Method]]))) %>% 
        summarise(N_MAGMA = sum(FDR_significant), 
                  N_eGENE = sum(eGENE),
                  N_overlap = sum(eGENE & FDR_significant),
                  OR = N_overlap/N_eGENE/N_MAGMA*n(),
                  Fisher_p = fisher.test(
                    table(eGENE, FDR_significant), alternative = "greater")$p.value) %>% 
        mutate(Cell.Type = CT, method = Method)
    }))
    return(rbind(df_all, df_CT))
  }), names(eGene[[1]]))
  Overlap_improvement <- do.call(rbind, Overlap_Methods) %>% 
    select(c("Trait_Identifier", "N_overlap", "method", "Cell.Type")) %>%
    group_by(Cell.Type) %>% 
    pivot_wider(names_from = "method", values_from = "N_overlap") %>% 
    left_join(MAGMA_Traits_select[,1:2], by = "Trait_Identifier") %>% 
    mutate(Combine_rel = (combine-Yazar)/Yazar,
           Static_rel = (static-Yazar)/Yazar)
  Overlap_improvement_df <- rbind(
    cbind(Overlap_improvement, rel = Overlap_improvement$Combine_rel, met = "scarf-QTL:combined"),
    cbind(Overlap_improvement, rel = Overlap_improvement$Static_rel, met = "scarf-QTL:static"))
  
  ggplot(filter(Overlap_improvement_df, Cell.Type != "All", rel != Inf), 
         aes(x = Trait.Name, y = rel))+
    geom_boxplot(aes(col = factor(met, levels = c("scarf-QTL:static", "scarf-QTL:combined"))), 
                 outliers = F)+
    geom_point(aes(x = as.numeric(factor(Trait.Name)) + 0.18*ifelse(met == "scarf-QTL:combined", 1, -1), 
                   y = rel, col = factor(Cell.Type, levels = corder)))+
    geom_point(data = filter(Overlap_improvement_df, Cell.Type == "All"), 
               aes(x = as.numeric(factor(Trait.Name)) + 0.18*ifelse(met == "scarf-QTL:combined", 1, -1), 
                   y = rel),
               shape = "*", size = 10)+
    geom_hline(yintercept = 0, col = "grey", linetype = "dotted")+
    geom_hline(yintercept = 0.3, col = "grey", linetype = "dotted")+
    geom_hline(yintercept = 1, col = "grey", linetype = "dotted")+
    scale_color_manual(name = "", values = c(
      tol14rainbow, `scarf-QTL:static` = "deepskyblue3",
      `scarf-QTL:combined` = "mediumpurple3"))+
    scale_y_continuous(breaks = c(seq(-1,3),0.3), labels = c(seq(-1,3),0.3)*100)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Immune diseases", y = "Relative improvement over Pseudobulk (%)")
  ggsave(sprintf("%s/More_magma_gene.pdf", DIR_fig), width = 9, height = 6)
}

# ===== Figure S6. QQ plots across cell types =====
{
  QQplots <- lapply(corder, function(CT) {
    ct <- stringr::str_replace(tolower(CT), " ", "")
    CT_TeX <- sprintf("$%s%s$", str_replace(CT, " ", "_{"), ifelse(str_detect(CT, " "), "}", ""))
    qq_df <- readRDS(sprintf("results/realdata/summary/qqplot_df/qqplot_sample1e5_%s.rds", ct))
    cutoff <- 200
    ggplot(filter(qq_df, method != "scar-QTL: dynamic"), 
           aes(expected, pmin(observed, cutoff), col = method)) +
      geom_point(size = 1, alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      theme_classic() +
      scale_color_manual(values = method_colors, guide = "none")+
      scale_y_continuous(breaks = 50*seq(0,4), labels = c(50*seq(0,3), ">200"), 
                         limits = c(0,200), expand = c(0,0,0,1))+
      labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = TeX(CT_TeX))
  })
  layout_mat3 <- matrix(c(seq(6), c(seq(7,11), NA), c(seq(12,14))),
                        byrow = T, ncol = 3)
  ggsave(plot = grid.arrange(grobs = QQplots, layout_matrix = layout_mat3), 
         height = 15, width = 12,
         filename = sprintf("%s/realdata_1e5.pdf", DIR_fig))
}

# ===== Figure 3b. Example dynamic eQTL boxplot =====
{
  if (FALSE) {
    parallel::mclapply(CT_14, function(CT) {
      ct <- str_replace(tolower(CT), " ", "")
      cell_meta <- readRDS(sprintf("data/covariate/cell_pt_%s.rds", ct))
      INDI_C <- cell_meta$individual; INDI <- sort(unique(INDI_C))
      TIME <- cell_meta$pseudotime
      TIME_quantile <- sapply(quantile(TIME, seq(5)/6), function(qt) TIME > qt)
      TIME_quantile <- factor(sprintf("Q%s", rowSums(TIME_quantile)+1), 
                              levels = sprintf("Q%s", seq(6)))
      PB_Q_mat <- do.call(cbind, parallel::mclapply(seq(22), function(chr) {
        Y_mat <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
        PB_mat <- sapply(INDI, function(indi) {
          Y_indi <- Y_mat[,INDI_C == indi, drop = F]
          #Y_indi <- log1p(Y_indi)
          PB_indi <- rowMeans(Y_indi)
        })
        PB_Q <- apply(Y_mat, 1, function(y) {
          tapply(y, interaction(TIME_quantile, INDI_C), mean)
        })
        cat(dim(PB_Q))
        return(PB_Q)
      }, mc.cores = 11))
      #table(rowSums(is.na(PB_Q_mat)))
      PB_Q_mat <- PB_Q_mat[rowSums(!is.na(PB_Q_mat))>0,]
      saveRDS(PB_Q_mat, file = sprintf("data/scRNA_SCT/PB_Q6_%s.rds", ct))
      return()
    }, mc.cores = 14)
  }
  if (FALSE) {
    eQTL_summary <- readRDS("scarf_QTL/Summary/eQTL_summary.rds")
    SNPs <- tapply(eQTL_summary$SNPID, eQTL_summary$Chromosome, unique)
    G_eSNP_mat <- parallel::mclapply(seq(22), function(chr) {
      G_SNP_mat <- readRDS(sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))
      G_SNP_mat_chr <- G_SNP_mat[SNPs[[chr]], ]
      cat(chr, "")
      return(G_SNP_mat_chr)
    }, mc.cores = 22)
    G_eSNP_mat <- do.call(rbind, G_eSNP_mat)
    saveRDS(G_eSNP_mat, file = "data/genotype/genotype_eSNPs.rds")
  }else{
    G_eSNP_mat  <- readRDS("data/genotype/genotype_eSNPs.rds")
  }
  # Helper function for dynamic eQTL boxplots
  plot_boxplot_eQTL <- function(GENE_ID, Cell.type, SNPID = NULL, RSID = NULL, Y_PB = NULL) {
    if (!is.null(SNPID)) {
      idx <- which(eQTL_summary$GENE_ID == GENE_ID & 
                     eQTL_summary$Cell.type == Cell.type & eQTL_summary$SNPID == SNPID)
    }else if (!is.null(RSID)) {
      idx <- which(eQTL_summary$GENE_ID == GENE_ID & 
                     eQTL_summary$Cell.type == Cell.type & eQTL_summary$RSID == RSID)
    }else{
      idx <- c()
    }
    if (length(idx)>0) {
      df <- eQTL_summary[idx,]
      ct <- df$CELL_ID; cTeX_ <- df$cTeX
      chr <- df$Chromosome; SNPID <- df$SNPID; RSID <- df$RSID
      if (is.null(Y_PB)) Y_PB <- readRDS(sprintf("data/scRNA_SCT/PB_Q6_%s.rds", ct))
      Y_pb <- Y_PB[,GENE_ID]
      G_SNP_vec <- G_eSNP_mat[SNPID,]
      
      deqtl_df <- data.frame(Expression = Y_pb, Q_I = rownames(Y_PB))
      deqtl_df <- cbind(str_split_fixed(deqtl_df$Q_I, "\\.", 2), deqtl_df)
      colnames(deqtl_df)[1:2] <- c("pseudotime_quantile", "INDI")
      deqtl_df <- left_join(deqtl_df, data.frame(INDI = colnames(G_eSNP_mat), Genotype = G_SNP_vec),
                            by = "INDI")
      {
        SNPinfo <- readRDS(sprintf("data/genotype/chr%s_SNP_info.rds", chr))
        SNPinfo_idx <- SNPinfo[SNPinfo$POS == str_split_fixed(SNPID, "[:_]", 3)[,2],]
        SNPinfo_idx$REF
        genotype_actg <- c(paste0(SNPinfo_idx$REF, SNPinfo_idx$REF),
                           paste0(SNPinfo_idx$REF, SNPinfo_idx$ALT),
                           paste0(SNPinfo_idx$ALT, SNPinfo_idx$ALT))
      }
      g <- ggplot(deqtl_df, aes(x = as.factor(pseudotime_quantile), 
                                fill = as.factor(Genotype), y = Expression))+
        geom_boxplot()+
        labs(x = TeX(sprintf("Pseudotime bin (%s cells)", cTeX_)), 
             y = sprintf("%s Expression", GENE_ID),
             fill = sprintf("%s\ngenotype", RSID))+
        scale_fill_manual(values = c("#4575b4", "#906bb1", "#d7301f"), labels = genotype_actg)+
        theme_classic()+
        theme(axis.ticks.y = element_blank(),        # 去掉刻度线
              axis.text.y = element_blank())
      attr(g, "name") <- paste0(df$Cell.type, "|", df$GENE_ID, "|", df$SNPID)
    }else{
      g <- NA
    }
    return(g)
  }
  
  b1 <- plot_boxplot_eQTL("FCER1A", "DC", "1:159223509_A")
  ggsave(b1, filename = sprintf("%s/%s.pdf", DIR_fig, attr(b1, "name")), 
         width = 9, height = 4.5)
}

# ===== Figure S7. DC PHATE and subset-marker visualization =====
{
  CT <- CT_ <- "DC"; ct <- "dc"
  {
    sds <- readRDS(sprintf("data/pseudotime/sds_%s.rds", CT_))
    pt <- slingCurves(sds)[[1]]$lambda
    em <- readRDS(sprintf("data/pseudotime/em_%s.rds", CT_))
    if (FALSE) {
      SCT <- lapply(seq(22), function(chr) {
        readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
      })
      SCT <- do.call(rbind, SCT)
      saveRDS(SCT, file = sprintf("data/scRNA_SCT/scRNA_SCT_%s_WG.rds", ct))
    }else{
      SCT <- readRDS(sprintf("data/scRNA_SCT/scRNA_SCT_%s_WG.rds", ct))
    }
    
    dc_marker_all <- xlsx::read.xlsx(file = "data/DC_subset/aah4573_supplementary_tables_1-16.xlsx", 
                                     sheetIndex = 2, startRow = 3)
    dc_marker_all <- dc_marker_all[,seq(11)] %>% filter(!is.na(AUC.value))
    dc_marker_top5 <- xlsx::read.xlsx(file = "data/DC_subset/aah4573_supplementary_tables_1-16.xlsx", 
                                      sheetIndex = 1, startRow = 3)
    dc_marker_top5 <- dc_marker_top5 %>% filter(!is.na(AUC.value))
    Subsets <- c("CLEC9A+", "CD1C_A", "CD1C_B", "CD141–CD1C–", 
                 "AXL+SIGLEC6+", "pDC")
  }
  
  {
    df <- data.frame(
      cell = colnames(SCT),
      sds@elementMetadata@listData[["reducedDim"]], 
      pseudotime = pt
    )
    cor_gene_pt <- apply(t(SCT), 2, cor, pt, method = "spearman")
    
    dc_marker_all <- dc_marker_all %>% group_by(Gene.ID) %>% 
      arrange(desc(AUC.value)) %>% slice_head(n = 1) %>% 
      ungroup() %>% arrange(as.integer(Rank)) %>% 
      mutate(Subset = factor(Associated.Cell.Population, 
                             levels = unique(dc_marker_all$Associated.Cell.Population),
                             labels = Subsets))
    dc_marker_sct <- inner_join(dc_marker_all, data.frame(
      Gene.ID = names(cor_gene_pt), cor = cor_gene_pt))
    exp_marker <- as.matrix(SCT)[dc_marker_sct$Gene.ID, ] %>%
      as.data.frame() %>%
      rownames_to_column("Gene.ID") %>%
      pivot_longer(-Gene.ID, names_to = "cell", values_to = "expr") %>%
      group_by(Gene.ID) %>%
      mutate(Expression = (expr - min(expr)) / (max(expr) - min(expr))) %>%
      ungroup()
    exp_marker <- left_join(exp_marker, dc_marker_all, by = ("Gene.ID"))
    exp_marker <- exp_marker %>%
      left_join(df, by = "cell")
    df_group_mean5 <- exp_marker %>%
      filter(Gene.ID%in%dc_marker_top5$Gene.ID) %>%
      group_by(Subset, cell) %>%
      summarise(mean_expression = mean(Expression, na.rm = TRUE), .groups = "drop")
    col.marker <- setNames(
      (colorRampPalette(brewer.pal(n = 6, name = "Set2"))(6))[c(6,3,4,2,1,5)],
      Subsets)
  }
  
  {
    g_pt <- ggplot(df, aes(x = PHATE1, y = PHATE2, col = pseudotime)) +
      geom_point(alpha = 1, shape = 20) +
      scale_color_gradientn(colours = col_fun) +
      geom_path(data = as.data.frame(slingCurves(sds)[[1]]$s), 
                aes(x = PHATE1, y = PHATE2), color = "red", linewidth = 1)+
      theme_minimal()
    
    g_cor <- ggplot(dc_marker_sct, aes(x = Subset, y = cor))+
      geom_hline(yintercept = 0, linetype = "dotted")+
      geom_violin(aes(fill = Subset)) + geom_boxplot(width = 0.2)+
      scale_fill_manual(values = col.marker)+
      guides(fill = "none")+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))+
      labs(x = "Discriminative gene of dendritic cell subsets",
           y = "Spearman correlation")
    
    g_top5 <- ggplot(df_group_mean5 %>% left_join(df, by = "cell"), 
                     aes(x = PHATE1, y = PHATE2, col = mean_expression)) +
      geom_point(alpha = 0.2, shape = 20) +
      scale_color_gradientn(colours = col_fun2,
                            breaks = c(0),
                            labels = c("0")) +
      theme_minimal() +
      guides(color = guide_colourbar(title = "expression", label = FALSE)) +
      facet_wrap(~Subset, nrow = 3)
    
    plot_grid(plot_grid(g_pt, g_cor, ncol = 1, labels = c("a","c")), 
              g_top5, nrow = 1, labels = c("","b"))
    ggsave(sprintf("%s/DC_phate.pdf", DIR_fig), width = 12, height = 8)
  }
}

# ===== Figure 3c. Heatmap of dynamic eQTL trajectories in DC =====
{
  CT <- "DC"
  {
    eQTL_df <- eQTL_summary[(eQTL_summary$Cell.type == CT) & (!is.na(eQTL_summary$beta_1)),]
    clu <- setNames(as.numeric(factor(eQTL_df$Cluster, levels = paste0("DC_clu", seq(5)))),
                    rownames(eQTL_df))
    
    marker_orig <- as.character(left_join(eQTL_df, dc_marker_all, 
                                          multiple = "first", by = c("GENE_ID" = "Gene.ID"))$Subset)
    names(marker_orig) <- rownames(eQTL_df)
    marker_orig[is.na(marker_orig)] <- "na"
    marker_orig_int <- as.numeric(factor(marker_orig, levels = c(names(col.marker), "na")))
    Yazar_orig <- ifelse(eQTL_df$GENE_ID %in% eQTL_Yazar$Gene.ID[eQTL_Yazar$Cell.type == "DC"], 
                         "eGene", "not")
    names(Yazar_orig) <- rownames(eQTL_df)
    
    coef <- as.matrix(eQTL_df[, sprintf("beta_%s", seq(K_Phi))])
    fit <- coef %*% t(phi)
    fit_normalized <- fit/rowMeans(abs(fit))
    fit_normalized1 <- fit_normalized/apply(abs(fit_normalized), 1, max)
    fit_normalized_mean <- rowMeans(fit_normalized1)
    fit_normalized_centered <- fit_normalized1-fit_normalized_mean
    mat.scale_ordered <- fit_normalized_centered[order(
      clu*1e4 + marker_orig_int*100 + (Yazar_orig != "eGene")*1000-
        sign(fit_normalized_mean)*10-abs(fit_normalized_mean)),]
  }
  {
    annotation_colors <- list()
    col.clu <- setNames(brewer.pal(8, 'Set1'), seq(8))[seq_len(length(unique(clu)))]
    annotation_colors$Cluster <- col.clu
    col.pseudotime <- setNames(
      colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(ncol(fit_normalized)),
      colnames(fit_normalized))
    annotation_colors$pseudotime <- col.pseudotime
    annotation_colors$Subset <- c(col.marker, na = "grey90")
    annotation_colors$Yazar <- c(eGene = "red", not = "grey")
    
    anno.col <- data.frame(pseudotime = colnames(fit_normalized))
    rownames(anno.col) <- colnames(mat.scale_ordered)
    
    mat.clust <- sort(clu)
    anno.row <- data.frame(Cluster = as.character(mat.clust))
    rownames(anno.row) <- names(mat.clust)
    anno.row$Subset <- marker_orig[rownames(anno.row)]
    anno.row$Yazar <- Yazar_orig[rownames(anno.row)]
    
    mat.scale_ordered <- mat.scale_ordered[
      !duplicated(paste0(mat.clust, str_split_fixed(rownames(mat.scale_ordered), ":", 2)[,1])),]
  }
  {
    P <- pheatmap(
      fit_normalized1[rownames(mat.scale_ordered), colnames(mat.scale_ordered)],
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      show_colnames = FALSE,
      color = colorRampPalette(colors = c("blue", "white","red"))(100),
      annotation_col = anno.col,
      annotation_row = anno.row,
      annotation_colors = annotation_colors,
      annotation_legend = F,
      cellwidth = 350 / ncol(mat.scale_ordered),
      cellheight = 350 / nrow(mat.scale_ordered),
      border_color = NA,
      silent = FALSE
    )
    ggsave(plot = grid.arrange(P[[4]], nrow = 1, ncol = 1), height = 6, width = 7,
           filename = sprintf("%s/%s_heatmap.pdf", DIR_fig, CT))
    
    P <- pheatmap(
      #mat.scale_ordered,
      fit_normalized1[1:5, 1:5],
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      show_colnames = FALSE,
      color = colorRampPalette(colors = c("blue", "white","red"))(100),
      annotation_col = NULL,
      annotation_row = anno.row,
      annotation_colors = annotation_colors,
      annotation_legend = T,
      cellwidth = 350 / ncol(mat.scale_ordered),
      cellheight = 350 / nrow(mat.scale_ordered),
      border_color = NA,
      main = "Legend",
      silent = FALSE
    )
    ggsave(plot = grid.arrange(P[[4]], nrow = 1, ncol = 1), height = 8, width = 4,
           filename = sprintf("%s/%s_heatmap_legend.pdf", DIR_fig, CT))
  }
}

# ===== Figure S8. Additional heatmaps in other cell types =====
{
  # for (CT in corder) {
  for (CT in c("CD4 NC", "CD8 S100B", "Mono C", "Plasma")) {
    {
      eQTL_df <- eQTL_summary[(eQTL_summary$Cell.type == CT) & (!is.na(eQTL_summary$beta_1)),]
      clu_levels <- paste0(CT, "_clu", seq(length(unique(eQTL_df$Cluster))))
      clu <- setNames(as.numeric(factor(eQTL_df$Cluster, levels = clu_levels)),
                      rownames(eQTL_df))
      
      Yazar_orig <- ifelse(eQTL_df$GENE_ID %in% eQTL_Yazar$Gene.ID[eQTL_Yazar$Cell.type == CT], 
                           "eGene", "not")
      names(Yazar_orig) <- rownames(eQTL_df)
      
      coef <- as.matrix(eQTL_df[, sprintf("beta_%s", seq(K_Phi))])
      fit <- coef %*% t(phi)
      fit_normalized <- fit/rowMeans(abs(fit))
      fit_normalized1 <- fit_normalized/apply(abs(fit_normalized), 1, max)
      fit_normalized_mean <- rowMeans(fit_normalized1)
      fit_normalized_centered <- fit_normalized1-fit_normalized_mean
      mat.scale_ordered <- fit_normalized_centered[order(
        clu*1e4 + (Yazar_orig != "eGene")*1000-
          sign(fit_normalized_mean)*10-abs(fit_normalized_mean)),]
    }
    {
      annotation_colors <- list()
      col.clu <- setNames(brewer.pal(8, 'Set1'), seq(8))[seq_len(length(unique(clu)))]
      annotation_colors$Cluster <- col.clu
      col.pseudotime <- setNames(
        colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(ncol(fit_normalized)),
        colnames(fit_normalized))
      annotation_colors$pseudotime <- col.pseudotime
      annotation_colors$Yazar <- c(eGene = "red", not = "grey")
      
      anno.col <- data.frame(pseudotime = colnames(fit_normalized))
      rownames(anno.col) <- colnames(mat.scale_ordered)
      
      mat.clust <- sort(clu)
      anno.row <- data.frame(Cluster = as.character(mat.clust))
      rownames(anno.row) <- names(mat.clust)
      anno.row$Yazar <- Yazar_orig[rownames(anno.row)]
      
      mat.scale_ordered <- mat.scale_ordered[
        !duplicated(paste0(mat.clust, str_split_fixed(rownames(mat.scale_ordered), ":", 2)[,1])),]
    }
    {
      P <- pheatmap(
        fit_normalized1[rownames(mat.scale_ordered), colnames(mat.scale_ordered)],
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = FALSE,
        show_colnames = FALSE,
        color = colorRampPalette(colors = c("blue", "white","red"))(100),
        annotation_row = anno.row,
        annotation_col = anno.col,
        annotation_colors = annotation_colors,
        annotation_legend = F,
        cellwidth = 350 / ncol(mat.scale_ordered),
        cellheight = 350 / nrow(mat.scale_ordered),
        border_color = NA,
        main = CT,
        silent = FALSE
      )
      ggsave(plot = grid.arrange(P[[4]], nrow = 1, ncol = 1), height = 6, width = 7,
             filename = sprintf("%s/more_heatmap_%s.pdf", DIR_fig, CT))
    }
  }
}

# ===== Figure S9. Example boxplots in CD4 NC =====
{
  eQTL_CT_cd4nc_dynamic <- eQTL_summary %>%
    filter(Cluster %in% c("CD4 NC_clu3", "CD4 NC_clu6")) %>%
    group_by(GENE_ID) %>% arrange(P_dynamic) %>% slice_head(n = 1) %>%
    ungroup() %>% arrange(P_dynamic)
  eQTL_CT_cd4nc_d <- eQTL_summary[eQTL_summary$Cluster %in% c("CD4 NC_clu3", "CD4 NC_clu6"),]
  eQTL_CT_cd4nc_d <- eQTL_CT_cd4nc_d[order(eQTL_CT_cd4nc_d$P_dynamic),]
  eQTL_CT_cd4nc_d <- eQTL_CT_cd4nc_d[!duplicated(eQTL_CT_cd4nc_d$GENE_ID),]
  
  Y_PB_cd4nc <- readRDS("data/scRNA_SCT/PB_Q6_cd4nc.rds")
  boxplot_cd4nc <- lapply(seq(10), function(i) {
    g <- plot_boxplot_eQTL(eQTL_CT_cd4nc_d$GENE_ID[i], "CD4 NC", eQTL_CT_cd4nc_d$SNPID[i],
                           Y_PB = Y_PB_cd4nc)
  })
  
  ggsave(grid.arrange(grobs = boxplot_cd4nc, ncol = 2), 
         filename = sprintf("%s/boxplot_cd4nc.pdf", DIR_fig),
         width = 12, height = 14)
}
