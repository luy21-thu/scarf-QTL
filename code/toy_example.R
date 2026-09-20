# ================================================================
# Minimal reproducible example for scarf-QTL
#
# This script generates a fully synthetic single-cell eQTL dataset
# and demonstrates the basic scarf-QTL workflow.
#
# The synthetic data are intended only to illustrate the required
# input format and software usage. They are independent of the
# simulation studies and real-data analyses reported in the paper.
# ================================================================

set.seed(20260920)

source("code/scarf_QTL_function.R")

if (file.exists("data/toy_example_data.RData")) {
  load("data/toy_example_data.RData")
} else {
  # ================================================================
  # 1. Generate individuals and cells
  # ================================================================
  
  n_ind <- 1000
  n_gene <- 5
  n_snp_per_gene <- 3000
  
  INDI <- sprintf("ID%04d", seq_len(n_ind))
  
  # Number of cells per individual.
  # The parameters can later be calibrated to roughly reflect
  # the scale observed in the OneK1K dendritic-cell data.
  n_cell_ind <- pmax(
    20,
    rnbinom(n_ind, mu = 200, size = 5)
  )
  
  INDI_C <- rep(INDI, times = n_cell_ind)
  n_cell <- length(INDI_C)
  
  cat("Number of individuals:", n_ind, "\n")
  cat("Number of cells:", n_cell, "\n")
  
  
  # ================================================================
  # 2. Generate individual-level covariates
  # ================================================================
  
  # Ten synthetic individual-level covariates, matching the dimension
  # used in the real-data analysis.
  #
  # cov1 is binary; the remaining covariates are continuous.
  # These variables do not correspond to actual OneK1K observations.
  
  W_INDI_raw <- cbind(
    rbinom(n_ind, 1, 0.5),
    rnorm(n_ind, mean = 0, sd = 1),
    matrix(rnorm(n_ind * 8), nrow = n_ind, ncol = 8)
  )
  
  colnames(W_INDI_raw) <- paste0("cov", seq_len(ncol(W_INDI_raw)))
  rownames(W_INDI_raw) <- INDI
  
  W_INDI_ <- t(scale(W_INDI_raw))
  
  
  # ================================================================
  # 3. Generate cell-level covariates and cell-state coordinates
  # ================================================================
  
  # A one-dimensional cell-state coordinate on [0, 1].
  TIME <- runif(n_cell, min = 0, max = 1)
  
  # Two synthetic cell-level covariates.
  W_CELL_raw <- cbind(
    rnorm(n_cell, mean = 0, sd = 1),
    rnorm(n_cell, mean = 0, sd = 1)
  )
  
  colnames(W_CELL_raw) <- c("cell_cov1", "cell_cov2")
  
  W_CELL <- t(scale(W_CELL_raw))
  
  
  # ================================================================
  # 4. Generate synthetic genotypes
  # ================================================================
  
  n_snp <- n_gene * n_snp_per_gene
  
  snp_names <- sprintf("SNP%04d", seq_len(n_snp))
  
  # Minor allele frequencies are generated independently for each SNP.
  maf <- runif(n_snp, min = 0.10, max = 0.50)
  
  G_SNP_mat <- t(
    vapply(
      maf,
      function(p) rbinom(n_ind, size = 2, prob = p),
      numeric(n_ind)
    )
  )
  
  rownames(G_SNP_mat) <- snp_names
  colnames(G_SNP_mat) <- INDI
  
  
  # Assign a separate set of cis-SNPs to each example gene.
  gene_names <- paste0("Gene", seq_len(n_gene))
  
  Gene_SNP_list <- lapply(
    seq_len(n_gene),
    function(g) {
      first <- (g - 1) * n_snp_per_gene + 1
      last  <- g * n_snp_per_gene
      snp_names[first:last]
    }
  )
  
  names(Gene_SNP_list) <- gene_names
  
  
  # ================================================================
  # 5. Genetic relationship matrix
  # ================================================================
  
  # For this minimal example, individuals are assumed to be unrelated.
  GRM_INDI <- diag(n_ind)
  
  rownames(GRM_INDI) <- INDI
  colnames(GRM_INDI) <- INDI
  
  
  # ================================================================
  # 6. Generate synthetic single-cell expression
  # ================================================================
  
  # Expression contains:
  #   - gene-specific baseline expression;
  #   - individual-level covariate effects;
  #   - cell-level covariate effects;
  #   - smooth cell-state effects;
  #   - individual-specific random effects;
  #   - cell-level noise.
  #
  # A cell-state-dependent genetic effect is introduced for
  # Gene1-SNP0001 to provide a known positive association.
  # No direct genetic effects are introduced for the remaining
  # gene-SNP pairs.
  
  Y_mat <- matrix(
    NA_real_,
    nrow = n_gene,
    ncol = n_cell,
    dimnames = list(gene_names, paste0("Cell", seq_len(n_cell)))
  )
  
  # Map each cell to its individual's row in W_INDI_raw.
  ind_index <- match(INDI_C, INDI)
  
  # Introduce one cell-state-dependent genetic effect for demonstration.
  causal_snp <- Gene_SNP_list[[1]][1]
  causal_genotype <- G_SNP_mat[causal_snp, ]
  
  # Standardize genotype to keep the effect size interpretable.
  causal_genotype <- as.numeric(scale(causal_genotype))
  
  for (g in seq_len(n_gene)) {
    
    # Gene-specific coefficients
    beta_ind <- rnorm(ncol(W_INDI_raw), mean = 0, sd = 0.15)
    beta_cell <- rnorm(ncol(W_CELL_raw), mean = 0, sd = 0.15)
    
    # Individual-specific random effects
    b_ind <- rnorm(n_ind, mean = 0, sd = 0.4)
    
    # Smooth cell-state effect
    state_effect <-
      0.5 * sin(2 * pi * TIME) +
      0.3 * TIME
    
    mu <-
      2 +
      as.numeric(W_INDI_raw[ind_index, , drop = FALSE] %*% beta_ind) +
      as.numeric(W_CELL_raw %*% beta_cell) +
      state_effect +
      b_ind[ind_index]
    
    # Add a cell-state-dependent genetic effect to Gene1-SNP0001.
    if (g == 1) {
      genetic_effect <-
        0.04 * causal_genotype[ind_index] * (2 * TIME - 0.5)
      
      mu <- mu + genetic_effect
    }
    
    Y_mat[g, ] <- mu + rnorm(n_cell, mean = 0, sd = 1)
  }
  
  
  # ================================================================
  # 7. Basic checks
  # ================================================================
  
  stopifnot(
    ncol(Y_mat) == length(INDI_C),
    length(TIME) == length(INDI_C),
    ncol(W_CELL) == length(INDI_C),
    ncol(W_INDI_) == length(INDI),
    ncol(G_SNP_mat) == length(INDI),
    all(colnames(G_SNP_mat) == INDI),
    all(rownames(GRM_INDI) == INDI),
    all(colnames(GRM_INDI) == INDI),
    all(names(Gene_SNP_list) == rownames(Y_mat))
  )
  
  save(Y_mat,
       INDI,
       INDI_C,
       W_INDI_,
       W_CELL,
       TIME,
       Gene_SNP_list,
       G_SNP_mat,
       GRM_INDI,
       file = "data/toy_example_data.RData")
}


# ================================================================
#  Run scarf-QTL
# ================================================================

system.time(
  res_Score <- scarfQTL(
    Y_mat = Y_mat,
    INDI = INDI,
    INDI_C = INDI_C,
    W_INDI_ = W_INDI_,
    W_CELL = W_CELL,
    TIME = TIME,
    Test_type = "RT",
    Gene_SNP_list = Gene_SNP_list,
    G_SNP_mat = G_SNP_mat,
    GRM_INDI = GRM_INDI,
    eta_hat0 = NULL,
    max_iter = 1000
  )
)

# ================================================================
#  Examine output
# ================================================================

names(res_Score)

example_gene <- names(res_Score$P_combined)[1]

head(data.frame(
  SNP = names(res_Score$P_combined[[example_gene]]),
  P_dynamic = res_Score$P_dynamic[[example_gene]],
  P_static = res_Score$P_static[[example_gene]],
  P_combined = res_Score$P_combined[[example_gene]]
))
