# ================================================================
# This script generates simulation figures and LaTeX-ready summaries
# for scarf-QTL evaluations, including:
#   (i) type I error tables,
#   (ii) power curves across dynamic effect patterns,
#   (iii) effect-size estimation summaries.
#   (iv) permutation-based null summaries,
#   (v) runtime benchmark figures.
# ================================================================

setwd("~/sceQTL")
library(cowplot)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(latex2exp)
library(stringr)

# ===== 1. Define directories, labels, and plotting constants =====
{
  DIR_summary <- "results/simulation/intermediate/summary"
  DIR_fig <- "results/simulation/figures"
  dir.create(DIR_fig, recursive = TRUE, showWarnings = FALSE)
  
  CT_simu <- c("CD4 NC", "B IN", "Plasma")
  names(CT_simu) <- str_replace(tolower(CT_simu), " ", "")
  cTeX_simu <- sprintf("$%s%s$", str_replace(CT_simu, " ", "_{"), 
                       ifelse(str_detect(CT_simu, " "), "}", ""))
  models <- c("normal", "poisson", "nbinom")
}

# ===== Table 1, S1 and S2. Type I error control across simulation settings =====
{
  T1E_df <- readRDS(sprintf("%s/T1E_df.rds", DIR_summary)) 
  T1E_df <- T1E_df %>% mutate(inflate = value>CR)
  
  print_method <- c("P_static", "P_dynamic", "P_combined", 
                    "P_static_Standard", "P_dynamic_Standard")
  alphas <- c(0.05, 0.1^(2:4))
  print(formatC(alphas+qnorm(0.99)*(sqrt(alphas*(1-alphas)/1e6)), digits = 3))
  sci_tex <- function(x, bf = F, digits = 2) {
    s <- formatC(x, format = "e", digits = digits)  
    if(!bf){
      tex <- sub("^(.*)e([+-]?)(0*)(\\d+)$", "$\\1\\\\times 10^{\\2\\4}$", s)
    }else{
      tex <- sub("^(.*)e([+-]?)(0*)(\\d+)$", "$\\\\mathbf{\\1\\\\times 10^{\\2\\4}}$", s)
    }
    return(tex)
  }
  model_table <- c("Normal", "Poisson", "\\makecell[c]{Negative \\\\ Binomial}")
  names(model_table) <- models
  for(CT in CT_simu[c(2,1,3)]){
    cat("\n", CT, "\n")
    for(model_ in models){
      cat("\\midrule\n\\multirow{", length(alphas), "}{*}{", model_table[model_], "}")
      for(alph in alphas){
        df <- left_join(data.frame(method = print_method),
                        filter(T1E_df, model == model_ & Cell.Type == CT & alpha == alph),
                        by = "method")
        t1e_print <- mapply(sci_tex, df$value, !df$inflate)
        cat("&", formatC(alph), paste("&", t1e_print), "\\\\", "\n")
      }
    }
  }
}

# ===== Figure 2. False positive eGenes under permutation-based calibration =====
{
  method_levels <- c("P_static", "P_dynamic", "P_combined", 
                     "P_static_Standard", "P_dynamic_Standard")
  method_labels <- c("scarf-QTL: static", "scarf-QTL: dynamic", "scarf-QTL: combined",
                     "Prospective: static", "Prospective: dynamic")
  method_colors <- setNames(
    c("deepskyblue3", "brown3", "mediumpurple3",
      "lightblue", "pink3"), method_labels)
  Permute_df <- readRDS(sprintf("%s/Permute_df.rds", DIR_summary))
  N_eGene_df <- Permute_df %>% 
    filter(parameter == "N_eGene") %>% 
    mutate(Method = factor(method, levels = method_levels, labels = method_labels)) %>% 
    filter(!is.na(Method))
  
  ggplot(N_eGene_df, aes(x = Method, y = log1p(value), col = Method))+
    geom_boxplot(outliers = T)+
    scale_color_manual(values = method_colors)+
    scale_y_continuous(breaks = log1p(c(0, 5^seq(0,4))), labels = c(0, 5^seq(0,4)))+
    facet_grid( ~ factor(Cell.Type, levels = CT_simu), scales = "free_y")+
    geom_hline(yintercept = 0, linetype = "dotted", col = "grey")+
    theme_classic()+
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(y = "Number of eGenes")
  ggsave(filename = sprintf("%s/Figure_permut_negene.pdf", DIR_fig), 
         width = 6, height = 5)
}

# ===== Figure 3, S1 and S2. Power comparison under dynamic genetic effect patterns =====
{
  # Define labels, methods, and plotting styles
  {
    Patterns_raw <- c("Static", "Linear", "TwoStage_onORoff", 
                      "Peak", "TwoStage_switch", 
                      "Cyclic_positive", "Cyclic_2period", "Cyclic_switch")
    Patterns_label <- c("Static", "Linear", "Two~stage~(On/Off)", 
                        "Peak", "Two~stage~(reverse)", 
                        "Cyclic", "Cyclic~(2~period)", "Cyclic~(reverse)")
    Patterns_label_abc <- c("Static", "Linear", "Two stage (On/Off)", 
                            "Peak", "Two stage (reverse)", 
                            "Cyclic", "Cyclic (2 period)", "Cyclic (reverse)")
    
    method_raw <- c("P_static", "P_dynamic", "P_combined", "P_Spearman")
    method_label <- c("scar-QTL: static", "scar-QTL: dynamic", "scar-QTL: combined", "Pseudobulk")
    method_colors4 <- c(
      "P_static"   = "deepskyblue3",  
      "P_dynamic"  = "brown3", 
      "P_combined" = "mediumpurple3",
      "P_Spearman" = "orange"
    )
    method_lty <- c(rep("solid", 3), "dashed")
    method_shap <- c(rep(21, 3), 23)
    names(method_colors4) <- names(method_lty) <- names(method_shap) <- method_label
  }
  
  # Load gamma-pattern and power summaries
  {
    gamma_df <- readRDS(sprintf("%s/Gamma_df.rds", DIR_summary))
    gamma_df <- gamma_df %>% mutate(
      label = "Pattern", 
      ABC = "\\textbf{A}", 
      pattern = factor(gamma_pattern, levels = Patterns_raw, labels = Patterns_label),
      abc = factor(gamma_pattern, levels = Patterns_raw, 
                   labels = sprintf("(\\textbf{%s}) %s", c("a", "b", "c", "d", "e", "f", "g", "h"),
                                    Patterns_label_abc))
    )
    
    gamma_df_text <- data.frame(
      pseudotime = 0, 
      value = max(gamma_df$value[gamma_df$pattern %in% Patterns_label[seq(4)]]),
      ABC = "\\textbf{A}",
      pattern = factor(Patterns_label[1], levels = Patterns_label),
      Ylab = "Effect~size",
      abc = "(\\textbf{a}) Static")
    
    Power_df <- readRDS(sprintf("%s/Power_df.rds", DIR_summary)) 
    Power_df <- Power_df %>% filter(alpha == 0.05) %>% 
      mutate(
        P_method = factor(method, levels = method_raw, labels = method_label),
        effect_size_x = factor(paste0(effect_size, "x"),
                               levels = paste0(sort(unique(Power_df$effect_size)), "x")),
        effect_size_break = ifelse(effect_size>5, effect_size/10+5, effect_size),
        effect_size_x10 = factor(
          paste0(ifelse(effect_size>5, effect_size/10, effect_size), "x"),
          levels = paste0(seq(5), "x")),
        CT = factor(Cell.Type, levels = CT_simu, labels = cTeX_simu),
        CT10 = factor(paste0(Cell.Type, ifelse(effect_size>5, "_10x", "")), 
                      levels = c(CT_simu, paste0(CT_simu, "_10x")), 
                      labels = c(cTeX_simu, paste0(cTeX_simu, " (10x)"))),
        ABC = c("\\textbf{B}", "\\textbf{C}", "\\textbf{D}")[as.numeric(CT)],
        ABC2 = c("\\textbf{A}", "\\textbf{B}", "\\textbf{C}")[as.numeric(CT)],
        pattern = factor(gamma_pattern, levels = Patterns_raw, labels = Patterns_label),
        abc = factor(gamma_pattern, levels = Patterns_raw, 
                     labels = sprintf("(\\textbf{%s}) %s", c("a", "b", "c", "d", "e", "f", "g", "h"),
                                      Patterns_label_abc))
      )
    
    Power_df_text <- data.frame(
      P_method = "",
      value = 1, effect_size_x = "1x", 
      CT_label = c("CD4[NC]", "B[IN]", "Plasma"),
      ABC = c("\\textbf{B}", "\\textbf{C}", "\\textbf{D}"),
      ABC2 = c("\\textbf{A}", "\\textbf{B}", "\\textbf{C}"),
      pattern = factor(Patterns_label[1], levels = Patterns_label),
      abc = "(\\textbf{a}) Static")
  }
  
  # Generate power figure for the normal setting
  {
    model_ <- "normal"
    p_pattern <- ggplot(gamma_df, aes(x = pseudotime, y = value)) +
      geom_area(aes(y = pmax(0, value)), fill = "lightblue", alpha = 0.5) +
      geom_area(aes(y = pmin(0, value)), fill = "pink2", alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
      geom_line() +
      facet_grid(abc ~ ABC, switch = "y",
                 labeller = as_labeller(TeX, default = label_parsed)) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(breaks = 0, position = "right") +
      scale_linetype_manual(values = c("solid", "dashed")) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.ticks.x = element_line(),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 10)) +
      labs(y = "") +
      guides(color = "none", linetype = "none") +
      geom_text(aes(label = Ylab), data = gamma_df_text, 
                parse = T, hjust = 0, vjust = 1, col = "black")
    
    p_power1 <- ggplot(filter(Power_df, model == model_ & effect_size < 10),
                       aes(x = effect_size_x, y = value,
                           group = P_method, col = P_method, 
                           linetype = P_method, shape = P_method)) +
      geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey50") +
      geom_line() + 
      geom_point(fill = "white") +
      scale_color_manual(values = method_colors4) +
      scale_linetype_manual(values = method_lty) +
      scale_shape_manual(values = method_shap) +
      theme_minimal() +
      facet_grid(pattern ~ ABC,
                 labeller = (CT = as_labeller(TeX, default = label_parsed)),
                 axis.labels = "margins") +
      theme(panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_rect(color = NA),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_blank(),
            legend.title = element_blank()
      ) +
      labs(x = "Effect size", y = TeX("% p-value < 0.05")) +
      geom_text(aes(label = CT_label), data = Power_df_text,
                parse = T, hjust = 0, vjust = 1,col = "black")
    
    plot_grid(p_pattern, p_power1, nrow = 1, 
              align = "hv", axis = "lr", rel_widths = c(0.8,3))
    ggsave(sprintf("%s/Figure_power_%s_all.pdf", DIR_fig, model_), 
           width = 12, height = 13)
  }
  
  # Generate power figures for count-based settings
  for (model_ in c("poisson", "nbinom")) {
    p_power2 <- ggplot(filter(Power_df, model == model_),
                       aes(x = effect_size_x, y = value,
                           group = P_method, col = P_method,
                           linetype = P_method, shape = P_method)) +
      geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey50") +
      geom_vline(xintercept = 5.5, linewidth = 0.2) +
      geom_line(data = filter(Power_df, model == model_ & effect_size < 10)) +
      geom_line(data = filter(Power_df, model == model_ & effect_size > 5)) +
      geom_point(fill = "white") +
      scale_color_manual(values = method_colors4) +
      scale_linetype_manual(values = method_lty) +
      scale_shape_manual(values = method_shap) +
      scale_x_discrete(expand = expansion(add = 0.4)) +
      scale_y_continuous(position = "right") +
      theme_minimal() +
      facet_grid(abc ~ ABC2, scales = "free_x", space = "free_x",
                 switch = "y", labeller = (CT = as_labeller(TeX, default = label_parsed)),
                 axis.labels = "margins") +
      theme(panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_rect(color = NA),
            strip.text.x = element_text(size = 15),
            legend.title = element_blank()) +
      guides(color = "none", linetype = "none", shape = "none") +
      labs(x = "Effect size", y = TeX("% p-value < 0.05")) +
      geom_text(aes(label = CT_label), data = Power_df_text,
                parse = T, hjust = 0, vjust = 1, col = "black")
    
    ggsave(plot = p_power2, 
           filename = sprintf("%s/Figure_power_%s.pdf", DIR_fig, model_), 
           width = 10, height = 12)
  }
}

# ===== Figure S3. Effect-size estimation accuracy in simulation studies =====
{
  Est_mean_curve <- readRDS(sprintf("%s/Estimation_df.rds", DIR_summary))
  Est_mean_curve <- Est_mean_curve %>% mutate(
    pseudotime = as.numeric(str_split_fixed(Est_mean_curve$coef, "_", 2)[,2]),
    CT = factor(Cell.Type, levels = CT_simu, labels = cTeX_simu),
    pattern = factor(gamma_pattern, levels = Patterns_raw, labels = Patterns_label),
    lambda = paste0(ifelse(model == "normal", "", paste0(model,"; ")), "\\lambda=", ridge_lamb)
  )
  Est_mean_curve_plot <- filter(
    Est_mean_curve, Cell.Type == "B IN", 
    (model == "normal" & effect_size == 5) |
      (model %in% c("nbinom", "poisson") & effect_size %in% c(5, 50) & ridge_lamb == 1))
  
  ggplot(Est_mean_curve_plot, aes(x = pseudotime, y = mean_est_scale)) +
    geom_hline(aes(yintercept = y), col = "white", 
               data = data.frame(
                 y = c(-0.5,2,3), # to adjust ylim
                 pattern = factor(Patterns_label[c(1,1,3)], levels = Patterns_label)))+
    geom_ribbon(aes(ymin = mean_est_scale - sd_est_scale,
                    ymax = mean_est_scale + sd_est_scale), 
                alpha = 0.2, fill = "brown3", 
                data = filter(Est_mean_curve_plot, model == "normal")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
    geom_line(aes(col = effect_size == 5), alpha = 0.9) +
    geom_line(aes(y = value), data = gamma_df, linetype = "dotted") +
    scale_color_manual(values = c("deepskyblue3", "brown3"), guide = NULL, na.value = NA) +
    facet_grid(pattern ~ lambda, switch = "y", 
               scales = "free_y", space = "free_y",
               labeller = (as_labeller(TeX, default = label_parsed))) +
    scale_x_continuous(breaks = seq(4) * 0.2, limits = c(0.05, 0.95)) +
    scale_y_continuous(breaks = seq(-5, 5), position = "right") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.ticks.y.right = element_line()) +
    labs(y = "Rescaled effect size", x = "Pseudotime")
  ggsave(filename = sprintf("%s/Figure_est_bin2.pdf", DIR_fig), 
         width = 9, height = 11)
}

# ===== Figure 4 and S10. Runtime benchmarking across simulation settings =====
{
  # Load and preprocess runtime summaries
  {
    Time_long_df <- readRDS(sprintf("%s/Runtime_df.rds", DIR_summary))
    
    Time_df_nsample <- filter(Time_long_df, N_gene==100, N_snp==3000) %>% 
      mutate(minutes = seconds/60)
    Time_df_nsample_RT <- filter(Time_df_nsample, Method != "Pseudobulk")
    df_ex <- do.call(rbind, lapply(unique(Time_df_nsample_RT$Step), function(i){
      D <- Time_df_nsample_RT[Time_df_nsample_RT$Step == i,]
      LM <- lm(log10(minutes)~log10(N_cell)*log10(N_indi), data = D)
      newD <- as.data.frame(t(rep(NA, ncol(D))))
      colnames(newD) <- colnames(D)
      newD$N_indi <- newD$N_cell <- 1e4
      newD$Step <- i
      newD$minutes <- 10^predict(LM, newD)
      return(newD)
    }))
    
    Time_df_ngene <- filter(Time_long_df, N_indi==1000, N_cell==200, N_snp==3000)
    Time_df_nsnp <- filter(Time_long_df, N_indi==1000, N_cell==200, N_gene==100)
    
    Time_df_nsample_comp <- Time_df_nsample %>%
      mutate(pointshape = ifelse(Method == "Pseudobulk", 23, 21)) %>%
      filter(step %in% c("total", "test"))
  }
  
  # Define runtime plotting scales and labels
  {
    min_breaks <- c(1/60, 10/60, 1, 10, 60, 60*10, 
                    60*24*seq(10))
    time_labels <- c("1s", "10s", "1min", "10min", "1hour", "10hour", 
                     paste0(seq(10), "day"))
    n_to_label <- function(x) {
      ifelse(x < 1e3, x, ifelse(x < 1e6, paste0(x / 1e3, "K"), paste0(x / 1e6, "M")))
    }
    linetyp <- c(rep("solid", length(levels(Time_long_df$Step))-1), "dashed")
    names(linetyp) <- levels(Time_long_df$Step)
    shap <- c(rep(21, length(levels(Time_long_df$Step))-1), 23)
    names(shap) <- levels(Time_long_df$Step)
  }
  
  # Generate sample-size contour plot
  {
    ggplot(rbind(Time_df_nsample_RT, df_ex), 
           aes(x = N_indi, y = N_cell, z = minutes)) +
      geom_contour_filled(breaks = c(0, min_breaks), alpha = 0.8) +
      scale_x_log10(breaks = unique(Time_df_nsample_$N_indi),
                    labels = n_to_label(unique(Time_df_nsample_$N_indi)))+
      scale_y_log10(breaks = unique(Time_df_nsample_$N_cell),
                    labels = n_to_label(unique(Time_df_nsample_$N_cell)))+
      scale_fill_brewer(labels = paste(c("0", time_labels[-length(time_labels)]), "-",
                                       time_labels),
                        name = "Running time", palette = "Spectral", direction = -1)+
      facet_grid( ~ Step)+
      new_scale_fill() +
      geom_point(aes(x = N_indi, y = N_cell, fill = log10(minutes),
                     shape = ifelse(!is.na(N_gene), "actual", "extrapolated")), 
                 col = "black", size = 3)+
      scale_fill_gradientn(
        breaks = log10(min_breaks)[seq(sum(max(Time_df_nsample_$minutes)>min_breaks)+1)],
        colors = rev(RColorBrewer::brewer.pal(
          sum(max(Time_df_nsample_$minutes)>min_breaks)+1, "Spectral")), 
        guide = "none")+
      scale_shape_manual(values = c(actual = 21, extrapolated = 24), 
                         guide = "none")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      labs(x = TeX("$n_{indi}$"), y = TeX("$n_{cell}$"), shape = NULL)+
      geom_point(
        data = filter(Time_df_nsample_, N_cell<1e3 & N_indi == 1e3 & step == "total"), 
        colour = "red",
        size = 5, shape=0,
      )
    ggsave(sprintf("%s/Figure_runtime_samplesize.pdf", DIR_fig), 
           width = 12, height = 3.3)
  }
  
  # Generate runtime scaling plots
  {
    p_ngene_pb <- ggplot(Time_df_ngene, aes(x = N_gene, y = seconds, col = Step))+
      geom_line(aes(linetype = Step))+
      geom_point(aes(shape = Step), fill = "white")+
      scale_x_log10(breaks = unique(Time_df_ngene$N_gene))+
      scale_y_log10(breaks = min_breaks*60, labels = time_labels)+
      scale_color_brewer(palette = "Set1")+
      scale_linetype_manual(values = linetyp)+scale_shape_manual(values = shap)+
      theme_classic()+
      guides(col = "none", linetype = "none", shape = "none")+
      labs(x = TeX("$n_{gene}$"), y = "Running time")
    
    p_nsnp_pb <- ggplot(Time_df_nsnp, aes(x = N_snp, y = seconds, col = Step))+
      geom_line(aes(linetype = Step))+
      geom_point(aes(shape = Step), fill = "white")+
      scale_y_log10(breaks = min_breaks*60, labels = time_labels)+
      scale_x_log10(breaks = unique(Time_df_nsnp$N_snp), 
                    labels = n_to_label(unique(Time_df_nsnp$N_snp)))+
      scale_color_brewer(palette = "Set1")+
      scale_linetype_manual(values = linetyp) + scale_shape_manual(values = shap)+
      theme_classic()+
      labs(x = TeX("$n_{snp}$"), y = "Running time")
    
    S3 <- ggplot(Time_df_nsample_comp, aes(x = N_indi, y = minutes))+
      geom_line(aes(col = log10(N_cell), 
                    linetype = Method, 
                    group = interaction(N_cell, Method)))+
      geom_point(aes(fill = log10(N_indi*N_cell), 
                     shape = pointshape), size = 3)+
      scale_x_log10(breaks = unique(Time_df_nsample_comp$N_indi),
                    labels = n_to_label(unique(Time_df_nsample_comp$N_indi)))+
      scale_y_log10(breaks = min_breaks, labels = time_labels)+
      scale_shape_identity(name = NULL)+
      scale_linetype_discrete(
        name = "Step", 
        labels = c("scarf-QTL (Total)", "Pseudobulk (Spearman\ncorrelation test)"))+
      scale_colour_gradient2(
        name = TeX("$n_{cell}$"),
        high = "Purple", low = "darkgreen", mid = "deepskyblue3", midpoint = 3, 
        breaks = log10(unique(Time_df_nsample_comp$N_cell)),
        labels = n_to_label(unique(Time_df_nsample_comp$N_cell)))+
      scale_fill_gradient2(
        name = "Total number of cells", 
        high = "red", low = "olivedrab", mid = "yellow",
        midpoint = 6, breaks = c(seq(4,7), log10(5)+7), 
        labels = TeX(c(sprintf("$10^%s$", seq(4,7)), "$5\\times 10^7$")))+
      theme_classic()+
      labs(x = TeX("$n_{indi}$"), y = "Running time")
    
    S12 <- plot_grid(p_ngene_pb, p_nsnp_pb, 
                     rel_widths = c(1,1.65), labels = c("a", "b"))
    S123 <- plot_grid(S12, S3, rel_widths = c(1,1), 
                      labels = c("", "c"), nrow = 2)
    ggsave(printf("%s/runtime_samplesize_pb.pdf", DIR_fig), 
           width = 9, height = 8)
  } 
}