# ================================================================
# scarf-QTL: core functions for single-cell eQTL analysis
#
# This file contains two main functions:
#   1) scarfQTL():  retrospective single-cell eQTL analysis with static,
#                   functional, and combined tests
#   2) PB_runtime(): pseudobulk-based eQTL analysis using Spearman correlation
#
# Required packages:
#   Matrix, stringr, qvalue, SparseM, data.table, RNOmni
# ================================================================

library(MASS)
library(Matrix)
library(splines)
library(SparseM)
library(RNOmni)    

# ------------------------------------------------
# INPUTS:
#   Y_mat      : gene × cell expression matrix (C columns aligned with INDI_C)
#   INDI       : vector of I unique individuals
#   INDI_C     : length-C vector mapping each cell to an individual
#
#   W_INDI_    : static individual-level covariates (covariate × individual)
#   W_CELL     : static cell-level covariates (covariate × cell)
#   Z_INDI_    : dynamic individual-level covariates (covariate × individual)
#   Z_CELL     : dynamic cell-level covariates (covariate × cell)
#
#   TIME       : pseudotime vector aligned with INDI_C
#   K_Phi      : number of B-spline basis functions
#
#   Gene_SNP_list :
#       list mapping each gene to its cis-SNPs;
#       names must be a subset of rownames(Y_mat)
#
#   G_SNP_mat  : SNP × I genotype matrix
#   GRM_INDI   : I × I genetic relatedness matrix
#
# OPTIONS:
#   Test_type     : character vector indicating the type of association test.
#                      NULL fits the null model only;
#                     "RT" performs the tests of scarf-QTL (retrospective score test);
#                     "Standard" performs the prospective score test for comparison.
#   Estimation    : whether to estimate smooth eQTL effect curves
#   ridge_lamb    : ridge penalty for effect-size estimation
#
#   Y_RINT        : whether to apply rank-based inverse normal transform
#   nc_min        : minimum number of cells per individual
#   eta_hat0      : initial values of eta_hat
#   lazy          : whether to use lazy version of IDUL^† updates
#   eps           : relative convergence tolerance
#   max_iter      : the maximum number of IDUL iterations
#   LM_cholesky   : whether to use SparseM::slm.fit for sparse WLS fitting
#   print_step    : whether to print algorithm steps
#   profile_steps : whether to return runtime of steps
# ------------------------------------------------

scarfQTL <- function(
    Y_mat, INDI, INDI_C,
    W_INDI_, W_CELL,
    Z_INDI_ = matrix(1, ncol = length(INDI), dimnames = list("baseline", INDI)), 
    Z_CELL = NULL,
    TIME, K_Phi = 10,
    Test_type = "RT", Gene_SNP_list, G_SNP_mat = NULL,
    GRM_INDI = diag(nrow = length(INDI)),
    Estimation = FALSE, ridge_lamb = NULL,
    Y_RINT = TRUE, nc_min = 1,
    eta_hat0 = NULL, lazy = TRUE, 
    eps = 1e-6, max_iter = 100, LM_cholesky = FALSE,
    print_step = FALSE, profile_steps = FALSE){
  # ===== 1. Internal helpers =====
  {
    Eigen_transform_indi <- function(indi){
      cell_indi <- which(INDI_C==indi); C_indi <- length(cell_indi)
      Phi_t_indi <- Phi_t[,cell_indi, drop=F]
      # eigen decomposition
      {
        n_U <- min(K_Phi, C_indi) # n_W <- C_indi-n_U
        svd_indi <- svd(t(Phi_t_indi), nu=C_indi, nv=0)
        U_indi <- svd_indi$u[,seq(n_U), drop=F]
        V_indi <- svd_indi$u[,-seq(n_U), drop=F]
        Lamb_indi <- (svd_indi$d)^2
      }
      # projection
      {
        Y_indi <- Y_mat[,cell_indi, drop=F]
        Y_U_indi <- t(as.matrix(Y_indi %*% U_indi))
        Y_V_indi <- t(as.matrix(Y_indi %*% V_indi))
        
        if (is.null(W_CELL)){
          W_U_CELL_indi <- W_V_CELL_indi <- NULL
        }else{
          W_CELL_indi <- W_CELL[,cell_indi, drop=F]
          W_U_CELL_indi <- t(as.matrix(W_CELL_indi %*% U_indi))
          W_V_CELL_indi <- t(as.matrix(W_CELL_indi %*% V_indi))
        }
        if (is.null(Z_CELL)){
          Z_U_CELL_indi <- Z_V_CELL_indi <- NULL
        }else{
          Z_CELL_indi <- Z_CELL[,cell_indi, drop=F]
          Phi_Z_CELL_indi <- mapply(outer, as.data.frame(Phi_t_indi), 
                                    as.data.frame(Z_CELL_indi))
          Z_U_CELL_indi <- t(Phi_Z_CELL_indi %*% U_indi)
          Z_V_CELL_indi <- t(Phi_Z_CELL_indi %*% V_indi)
          colnames(Z_U_CELL_indi) <- colnames(Z_V_CELL_indi) <- paste0(
            rep(rownames(Z_CELL), each=K_Phi), ".", seq(K_Phi))
        }
        X_U_CELL_indi <- cbind(W_U_CELL_indi, Z_U_CELL_indi)
        X_V_CELL_indi <- cbind(W_V_CELL_indi, Z_V_CELL_indi)
        if (is.null(Z_V_CELL_indi)) X_V_CELL_indi <- W_V_CELL_indi
        Phi_U_indi <- t(Phi_t_indi %*% U_indi)
        # Phi_U_indi <- svd_indi$d * t(svd_indi$v)
        colSum_U_indi <- colSums(U_indi)
      }
      return(list(Lamb = Lamb_indi, Y_U = Y_U_indi, Y_V = Y_V_indi, 
                  X_U_CELL = X_U_CELL_indi, X_V_CELL = X_V_CELL_indi, 
                  Phi_U = Phi_U_indi, colSum_U = colSum_U_indi))
    }
    IDUL_iterations <- function(Y_, X_, D, eta_0, lazy=T, eps=1e-6, max_iter=100, 
                                g=NULL, LM_cholesky=F){
      if (max_iter!=0){
        eta_t <- eta_0
        H_t <- eta_t * D + 1
        if (LM_cholesky){
          WLS1_s <- H_t * (slm.wfit(X_, Y_, weights = 1/H_t)$residuals)^2
        }else{
          WLS1_s <- (lm(Y_~0+X_, weights = 1/H_t)$residuals)^2
        }
        tau_inv_t <- mean(WLS1_s/H_t)
        LL_t <- -sum(log(H_t))-length(D)*log(tau_inv_t)
        for(t in 1:max_iter){
          WLS2 <- lm(WLS1_s~D, weights = 1/H_t^2)
          eta_t1 <- WLS2$coefficients[2] / tau_inv_t +
            (1 - WLS2$coefficients[1] / tau_inv_t)*eta_t/(1+lazy)
          # eta_t1 <- WLS2$coefficients[2]/WLS2$coefficients[1]
          eta_t1 <- max(eta_t1, 0)
          while(abs(eta_t1-eta_t)/(abs(eta_t)+eps)>eps){
            H_t1 <- eta_t1 * D + 1
            if (LM_cholesky){
              WLS1_t1 <- slm.wfit(X_, Y_, weights = 1/H_t)
              WLS1_s_t1 <- H_t * (WLS1_t1$residuals)^2
            }else{
              WLS1_t1 <- lm(Y_~0+X_, weights = 1/H_t)
              WLS1_s_t1 <- (WLS1_t1$residuals)^2
            }
            tau_inv_t1 <- mean(WLS1_s_t1/H_t1)
            LL_t1 <- -sum(log(H_t1))-length(D)*log(tau_inv_t1)
            if (LL_t1>=LL_t) break
            eta_t1 <- (eta_t1+eta_t)/2
          }
          if (abs(eta_t1-eta_t)/(abs(eta_t)+eps)<=eps) break
          eta_t <- eta_t1; H_t <- H_t1
          WLS1_s <- WLS1_s_t1
          tau_inv_t <- tau_inv_t1; LL_t <- LL_t1
        }
        if (t==max_iter){
          mle_wls <- mle_wls_NA
          cat(paste("Fail to converge for gene", g, "!\n"))
        }else{
          if (print_step) cat(sprintf("t=%s; ", t))
          eta_hat <- eta_t1
          H_hat <- eta_hat * D + 1
          if (LM_cholesky){
            MLE_wls <- slm.wfit(X_, Y_, weights = 1/H_hat)
            resids <- sqrt(H_hat) * MLE_wls$residuals
            sigma_hat <- sqrt(sum(MLE_wls$residuals^2)/MLE_wls$df.residual)
            coef <- MLE_wls$coefficients; names(coef) <- colnames(X_)
          }else{
            MLE_wls <- lm(Y_~0+X_, weights = 1/H_hat)
            resids <- MLE_wls$residuals
            sigma_hat <- summary(MLE_wls)$sigma
            coef <- MLE_wls$coefficients
          }
          mle_wls <- list(
            eta_hat = eta_hat, weights = 1/H_hat, residuals = resids, 
            sigma_hat = sigma_hat, coef = coef)
        }
      }else{
        if (is.na(eta_0)){
          mle_wls <- mle_wls_NA
        }else{
          eta_hat <- eta_0
          H_hat <- eta_hat * D + 1
          if (LM_cholesky){
            MLE_wls <- slm.wfit(X_, Y_, weights = 1/H_hat)
            resids <- sqrt(H_hat) * MLE_wls$residuals
            sigma_hat <- sqrt(sum(MLE_wls$residuals^2)/MLE_wls$df.residual)
            coef <- MLE_wls$coefficients; names(coef) <- colnames(X_)
          }else{
            MLE_wls <- lm(Y_~0+X_, weights = 1/H_hat)
            resids <- MLE_wls$residuals
            sigma_hat <- summary(MLE_wls)$sigma
            coef <- MLE_wls$coefficients
          }
          mle_wls <- list(eta_hat = unname(eta_hat), weights = 1/H_hat, residuals = resids, 
                          sigma_hat = sigma_hat, coef = coef)
        }
      }
      return(mle_wls)
    }
    CCT_2vec <- function(P1, P2, name = names(P1)){
      if (length(P1) != length((P2))) stop("Length not matched!")
      P1[is.na(P1)] <- P2[is.na(P2)] <- 1
      P_combin <- rep(-1, length(P1))
      P_combin[P1==1] <- P2[P1==1]; P_combin[P2==1] <- P1[P2==1]
      P_combin[P1==0 | P2==0] <- 0
      
      p1 <- P1[P_combin==-1]; p2 <- P2[P_combin==-1]
      q1 <- ifelse(p1 < 1e-16, 1/p1/pi, tan((0.5-p1)*pi))
      q2 <- ifelse(p2 < 1e-16, 1/p2/pi, tan((0.5-p2)*pi))
      cct.stat <- (q1+q2)/2 # weighted mean
      pval <- ifelse(cct.stat > 1e+15, (1/cct.stat)/pi, 1-pcauchy(cct.stat))
      P_combin[P_combin==-1] <- pval
      names(P_combin) <- name
      return(P_combin)
    }
  }
  
  # ===== 2. Preparation and input reformatting =====
  {
    if (print_step) cat("Preprocessing...\n")
    tick <- function() proc.time()[["elapsed"]]
    tock <- function(t0) proc.time()[["elapsed"]] - t0
    step_times <- list()    
    t0_step <- if (profile_steps) tick() else NULL
  }
  {
    if (is.vector(Y_mat)) Y_mat <- matrix(Y_mat, nrow=1)
    if (!identical(names(Gene_SNP_list), rownames(Y_mat))){
      Y_mat <- Y_mat[names(Gene_SNP_list), , drop=F]
    } 
    if (is.vector(G_SNP_mat)) G_SNP_mat <- matrix(G_SNP_mat, nrow=1)
    {
      INDI <- INDI[table(INDI_C)>=nc_min]
      if (sum(table(INDI_C)<nc_min)){
        W_INDI_ <- W_INDI_[,INDI, drop=F]
        Z_INDI_ <- Z_INDI_[,INDI, drop=F]
        if (!is.null(G_SNP_mat)) G_SNP_mat <- G_SNP_mat[,INDI, drop=F]
        GRM_INDI <- GRM_INDI[INDI, INDI, drop=F]
      }
    }
    if (Y_RINT){
      Y_mat <- t(apply(Y_mat, 1, function(Y){
        if (var(Y)==0){
          Y_ <- Y*0
        }else{
          Y_ <- RankNorm(Y)
          Y_ <- scale(Y_)*sd(Y) + mean(Y)
        }
        return(Y_)
      }))
    }
    if (sum((diag(GRM_INDI)-1)^2)>1e-10){
      std_devs <- sqrt(diag(GRM_INDI))
      GRM_INDI <- GRM_INDI / (std_devs %*% t(std_devs))
      warning("Input GRM not a correlation matrix! \n")
    } 
  }
  
  # ===== 3. Construct B-spline basis =====
  {
    TIME_ <- rank(TIME)/(length(TIME)+1)
    if (K_Phi>3){
      Phi <- bs(TIME_, knots=seq(0,1,length.out=K_Phi-3+1), Boundary.knots=c(0,1))
      Phi <- Phi[,-ncol(Phi), drop=F]; colnames(Phi) <- sprintf("bs%s", seq(K_Phi))
    }else{
      if (K_Phi=="linear"){
        Phi <- cbind(TIME_, 1-TIME_)
      }else{
        stop("Number of splines should be more than 3!")
      }
    }
    Phi_t <- t(Phi)
    rm(TIME, TIME_, Phi)
  }
  if (profile_steps) step_times$prepare <- tock(t0_step)
  
  # ===== 4. Eigen transformation =====
  {
    if (print_step) cat("Eigen decomposition...\n")
    t0_step <- if (profile_steps) tick() else NULL
    eigen_transform <- lapply(INDI, Eigen_transform_indi)
    rm(Y_mat, Z_CELL, W_CELL, Phi_t)
    # rbind of individuals
    {
      Lamb_list <- lapply(eigen_transform, "[[", "Lamb")
      Lamb_raw <- do.call(c, Lamb_list)
      N_U <- length(Lamb_raw); N_W <- sum(INDI_C %in% INDI)-N_U
      D_raw <- c(Lamb_raw, rep(0, N_W))
      D_scale <- max(D_raw); D <- D_raw/D_scale
      if (!is.null(Test_type)){
        f_U <- unlist(mapply(function(L, i){
          rep(i, length(L))
        }, Lamb_list, INDI, SIMPLIFY=F))
        f_Q <- c(f_U, unlist(mapply(function(L, i){
          rep(i, nrow(L))
        }, lapply(eigen_transform, "[[", "X_V_CELL"), INDI, SIMPLIFY=F)))
      }
      Y_Q_mat <- rbind(do.call(rbind, lapply(eigen_transform, "[[", "Y_U")), 
                       do.call(rbind, lapply(eigen_transform, "[[", "Y_V")))
      
      Phi_U_list <- lapply(eigen_transform, "[[", "Phi_U")
      colSum_U_list <- lapply(eigen_transform, "[[", "colSum_U")
      Z_U_INDI <- do.call(rbind, mapply(
        function(z, pu) kronecker(t(z), pu), 
        as.data.frame(Z_INDI_), Phi_U_list, SIMPLIFY=F))
      colnames(Z_U_INDI) <- paste0(
        rep(rownames(Z_INDI_), each=K_Phi), ".", seq(K_Phi))
      W_U_INDI <- do.call(rbind, mapply(
        function(z, pu) kronecker(t(z), pu), 
        as.data.frame(W_INDI_), 
        colSum_U_list, SIMPLIFY=F))
      colnames(W_U_INDI) <- rownames(W_INDI_)
      X_U_INDI <- cbind(W_U_INDI, Z_U_INDI)
      
      X_U_COV <- cbind(do.call(rbind, lapply(eigen_transform, "[[", "X_U_CELL")), 
                       X_U_INDI)
      X_V_CELL <- do.call(rbind, lapply(eigen_transform, "[[", "X_V_CELL"))
      if (!is.null(Test_type)){
        XtX_V_CELL <- t(X_V_CELL)%*%X_V_CELL / N_U
      }
      X_Q_COV <- rbind(X_U_COV, cbind(X_V_CELL, 
                                      matrix(0, nrow=N_W, ncol=ncol(X_U_INDI))))
      if (LM_cholesky){
        X_Q_COV_s <- as(X_Q_COV, "matrix.csr")
        rm(X_Q_COV)
      }
    }
    rm(Z_INDI_, W_INDI_, W_U_INDI, Z_U_INDI, X_U_INDI, X_V_CELL, eigen_transform); gc()
    if (profile_steps) step_times$eigendecomp <- tock(t0_step)
  }
  
  # ===== 5. Null model fitting =====
  {
    if (print_step) cat("Fitting null MLE for", length(Gene_SNP_list), "genes...\n")
    t0_step <- if (profile_steps) tick() else NULL
    MLE_0_list <- vector("list", length(Gene_SNP_list))
    names(MLE_0_list) <- names(Gene_SNP_list)
    mle_wls_NA <- list(eta_hat=NA, weights=NA, residuals=NA, sigma_hat=NA, coef=NA)
    for(g in names(Gene_SNP_list)){
      Y_Q <- Y_Q_mat[,g]
      if (var(Y_Q)!=0){
        if (print_step) cat(g, " ")
        if (!is.null(eta_hat0)){
          eta_0 <- eta_hat0[g]
        }else{
          eta_0 <- 1
        }
        if (!LM_cholesky){
          MLE_0 <- IDUL_iterations(
            Y_Q, X_Q_COV, D, eta_0, lazy, eps, max_iter, g=g)
        }else{
          MLE_0 <- IDUL_iterations(
            Y_Q, X_Q_COV_s, D, eta_0, lazy, eps, max_iter, g=g, LM_cholesky=T)
        }
        MLE_0_list[[g]] <- MLE_0
      }else{
        MLE_0_list[[g]] <- mle_wls_NA
      }
    }
    eta_hat0 <- sapply(MLE_0_list, "[[", "eta_hat")
    res <- list(eta_hat0=eta_hat0)
    if (profile_steps) step_times$null <- tock(t0_step)
    if (!LM_cholesky) rm(X_Q_COV) else rm(X_Q_COV_s)
    rm(Y_Q_mat, MLE_0, Y_Q); gc()
  }
  
  # ===== 6. Retrospective association test =====
  if ("RT" %in% Test_type){
    if (print_step) cat("\nRetrospective score test...\n")
    t0_step <- if (profile_steps) tick() else NULL
    P_dynamic <- P_static <- vector("list", length(Gene_SNP_list))
    names(P_dynamic) <- names(P_static) <- names(Gene_SNP_list)
    if (Estimation) Estimation_beta <- P_dynamic
    for(g in names(Gene_SNP_list)){
      if (print_step) cat(g, " ")
      MLE_0 <- MLE_0_list[[g]]
      if (is.na(MLE_0)[1]){
        P_dynamic[[g]] <- P_static[[g]] <- rep(NA, length(Gene_SNP_list[[g]]))
        names(P_dynamic[[g]]) <- names(P_static[[g]]) <- Gene_SNP_list[[g]]
      }else{
        Weight_U_hat0 <- MLE_0$weights[1:N_U]
        Gamma_U <- split(Weight_U_hat0 * MLE_0$residuals[1:N_U], f_U)[INDI]
        Phi_Gamma_U <- mapply("%*%", Gamma_U, Phi_U_list)
        G_SNP_mat_g <- G_SNP_mat[Gene_SNP_list[[g]], , drop=F]
        Score_U_mat <- G_SNP_mat_g %*% t(Phi_Gamma_U)# / sqrt(N_U)
        
        {
          MAF <- rowMeans(G_SNP_mat_g)/2
          var_MAF <- 2*MAF*(1-MAF)
          # V <- Phi_Gamma_U %*% t(Phi_Gamma_U) / N_U
          V <- Phi_Gamma_U %*% GRM_INDI %*% t(Phi_Gamma_U)# / N_U
          V_inv <- ginv(V)
          # Score_Stat_noV <- apply(Score_U_mat, 1, function(U) U %*% V_inv  %*% U)
          Score_Stat_noV <- rowSums((Score_U_mat%*%V_inv) * Score_U_mat)
          Score_Stat <- 1/var_MAF * Score_Stat_noV
          Score_P <- pchisq(Score_Stat, K_Phi, lower.tail=F)
        }
        {
          Score_U_constant <- rowSums(Score_U_mat)
          V_constant <- sum(V)
          Score_Stat_constant <- Score_U_constant^2/V_constant/var_MAF
          Score_P_constant <- pchisq(Score_Stat_constant, 1, lower.tail=F)
        }
        P_dynamic[[g]] <- Score_P
        P_static[[g]] <- Score_P_constant
        names(P_dynamic[[g]]) <- names(P_static[[g]]) <- rownames(G_SNP_mat_g)
        if (Estimation){
          Weight_U <- split(Weight_U_hat0, f_U)[INDI]
          Phi_W_Phi <- mapply(function(a, b){
            t(a) %*% (b*a)
          }, Phi_U_list, Weight_U)
          XtVX <- G_SNP_mat_g^2 %*% t(Phi_W_Phi)
          if (length(ridge_lamb)>0){
            ridge_D <- diff(diag(K_Phi), differences = 2)
            ridge_P <- t(ridge_D) %*% ridge_D
            # XtVX <- XtVX + ridge_P
            beta_hat_SNP <- lapply(ridge_lamb, function(lamb){
              sapply(rownames(Score_U_mat), function(s){
                ginv(matrix(XtVX[s,], nrow=K_Phi) + lamb * ridge_P) %*% Score_U_mat[s,]
              })
            })
            names(beta_hat_SNP) <- ridge_lamb
          }else{
            beta_hat_SNP <- sapply(rownames(Score_U_mat), function(s){
              ginv(matrix(XtVX[s,], nrow=K_Phi)) %*% Score_U_mat[s,]
            })
          }
          Estimation_beta[[g]] <- beta_hat_SNP
        }
      }
    }
    P_combined <- mapply(CCT_2vec, P_dynamic, P_static, SIMPLIFY = F)
    res <- c(res, list(P_dynamic = P_dynamic, P_static = P_static, P_combined = P_combined))
    if (profile_steps) step_times$retro <- tock(t0_step)
    if (Estimation)  res <- c(res, list(Estimation = Estimation_beta))
  }
  
  # ===== 7. Prospective score test =====
  if ("Standard" %in% Test_type){
    if (print_step) cat("\nStandard score test...\n")
    t0_step <- if (profile_steps) tick() else NULL
    P_dynamic_Standard <- P_static_Standard <- vector("list", length(Gene_SNP_list))
    names(P_dynamic_Standard) <- names(P_static_Standard) <- names(Gene_SNP_list)
    for(g in names(Gene_SNP_list)){
      if (print_step) cat(g, " ")
      MLE_0 <- MLE_0_list[[g]]
      if (is.na(MLE_0)[1]){
        P_dynamic_Standard[[g]] <- P_static_Standard[[g]] <- rep(NA, length(Gene_SNP_list[[g]]))
        names(P_dynamic_Standard[[g]]) <- names(P_static_Standard[[g]]) <- Gene_SNP_list[[g]]
      }else{
        Weight_U_hat0 <- MLE_0$weights[1:N_U]
        Gamma_U <- split(Weight_U_hat0 * MLE_0$residuals[1:N_U], f_U)[INDI]
        Phi_Gamma_U <- mapply("%*%", Gamma_U, Phi_U_list)
        G_SNP_mat_g <- G_SNP_mat[Gene_SNP_list[[g]], , drop=F]
        Score_U_mat <- G_SNP_mat_g %*% t(Phi_Gamma_U) / 
          (sqrt(N_U) * MLE_0$sigma_hat)
        
        I_11 <- t(X_U_COV) %*% (Weight_U_hat0 * X_U_COV) / N_U
        C_idx <- seq(ncol(XtX_V_CELL))
        I_11[C_idx, C_idx] <- I_11[C_idx, C_idx] + XtX_V_CELL
        
        I_21_noG <- mapply(function(w_x_cov, phi_u){
          t(phi_u) %*% matrix(w_x_cov, ncol=ncol(X_U_COV))
        }, split(Weight_U_hat0*X_U_COV, f_U)[INDI], Phi_U_list)
        I_21_mat <- G_SNP_mat_g %*% t(I_21_noG) / N_U
        
        I_22_noG <- mapply(function(w, phi_u){
          t(phi_u) %*% (w * phi_u)
        }, split(Weight_U_hat0, f_U)[INDI], Phi_U_list)
        I_22_mat <- G_SNP_mat_g^2 %*% t(I_22_noG) / N_U
        
        P_dynamic_Standard[[g]] <- sapply(seq(nrow(G_SNP_mat_g)), function(SNP_idx){
          I_21=matrix(I_21_mat[SNP_idx,], ncol=ncol(X_U_COV))
          I_22=matrix(I_22_mat[SNP_idx,], ncol=K_Phi)
          I <- rbind(cbind(I_11, t(I_21)),
                     cbind(I_21, I_22))
          S_idx <- seq(ncol(I_11)+1, ncol(I))
          {
            I_inv <- NA
            try(I_inv <- ginv(I))
            if (is.na(I_inv)[1]) I_inv <- solve(I)
          }
          Score_I <- I_inv[S_idx,S_idx]
          Score_U <- Score_U_mat[SNP_idx,]
          Score_Stat <- as.vector(t(Score_U) %*% Score_I %*% Score_U)
          Score_P <- pchisq(Score_Stat, K_Phi, lower.tail=F)
        })
        names(P_dynamic_Standard[[g]]) <- rownames(G_SNP_mat_g)
        
        P_static_Standard[[g]] <- sapply(seq(nrow(G_SNP_mat_g)), function(SNP_idx){
          I_21=colSums(matrix(I_21_mat[SNP_idx,], ncol=ncol(X_U_COV)))
          I_22=sum(I_22_mat[SNP_idx,])
          I <- rbind(cbind(I_11, I_21), c(I_21, I_22))
          S_idx <- ncol(I_11)+1
          {
            I_inv <- NA
            try(I_inv <- ginv(I))
            if (is.na(I_inv)[1]) I_inv <- solve(I)
          }
          Score_I <- I_inv[S_idx,S_idx]
          Score_U <- sum(Score_U_mat[SNP_idx,])
          Score_Stat <- Score_U^2 * Score_I
          Score_P <- pchisq(Score_Stat, 1, lower.tail=F)
        })
        names(P_static_Standard[[g]]) <- rownames(G_SNP_mat_g)
      }
    }
    P_combined_Standard <- mapply(CCT_2vec, P_dynamic_Standard, P_static_Standard, SIMPLIFY = F)
    if (profile_steps) step_times$prosp <- tock(t0_step)
    res <- c(res, list(P_dynamic_Standard = P_dynamic_Standard, 
                       P_static_Standard = P_static_Standard,
                       P_combined_Standard = P_combined_Standard))
  }
  
  attr(res, "profile") <- list(
    step_times = step_times
  )
  return(res)
}


# ------------------------------------------------
# PB_runtime
#
# Perform pseudobulk-based eQTL testing using
# Spearman correlation between genotype and
# covariate-adjusted pseudobulk expression.
#
# INPUTS:
#   Y_mat          : gene × cell expression matrix
#   INDI           : vector of individual IDs
#   INDI_C         : individual index for each cell
#   W_INDI_        : individual-level covariates (covariate × individual)
#   Gene_SNP_list  : named list of cis-SNPs for each gene
#   G_SNP_mat      : SNP × individual genotype matrix
#
# OPTIONS:
#   Log2p          : whether to apply log1p to expression before aggregation
#   profile_steps  : whether to record runtime of each step
#   return_corr    : whether to return correlation coefficients
#
# OUTPUT:
#   A named list of p-values (or p-values + correlations),
#   one element per gene, with SNP-level results.
# ------------------------------------------------

PB_runtime <-
  function(Y_mat, INDI, INDI_C, W_INDI_, Gene_SNP_list, G_SNP_mat, # same as scarfQTL
           Log2p = TRUE, # whether to do log transform 
           profile_steps = T, # whether to record runtime
           return_corr = F # whether to return correlation coefficient
  ){
    {
      tick <- function()
        proc.time()[["elapsed"]]
      tock <- function(t0)
        proc.time()[["elapsed"]] - t0
      step_times <- list()
    }
    
    {
      t0_step <- if (profile_steps) tick() else NULL
      PB_mat <- sapply(INDI, function(indi) {
        Y_indi <- Y_mat[, INDI_C == indi, drop = F]
        if (Log2p)
          Y_indi <- log1p(Y_indi)
        PB_indi <- rowMeans(Y_indi)
      })
      if (profile_steps)
        step_times$prepare <- tock(t0_step)
    }
    {
      t0_step <- if (profile_steps) tick() else NULL
      PB_res_mat <-
        lm.fit(x = cbind(1, t(W_INDI_)), y = t(PB_mat))$residuals
      if (profile_steps)
        step_times$null <- tock(t0_step)
    }
    {
      t0_step <- if (profile_steps) tick() else NULL
      res <- lapply(names(Gene_SNP_list), function(G) {
        df <- nrow(PB_res_mat) - 2
        PB_res <- PB_res_mat[, G]
        SNP_g <- Gene_SNP_list[[G]]
        G_SNP_g <- G_SNP_mat[SNP_g, , drop = F]
        
        PB_rank <- rank(PB_res, ties.method = "average")
        G_rank <- apply(G_SNP_g, 1, rank, ties.method = "average")
        
        r <- suppressWarnings(cor(G_rank, PB_rank, method = "pearson"))
        r <- as.vector(r)
        tstat <- r * sqrt(df / pmax(1e-12, 1 - r ^ 2))
        pval  <- 2 * pt(abs(tstat), df = df, lower.tail = FALSE)
        names(pval) <- SNP_g
        
        if (return_corr)
          ret <- list(pval = pval, r = r)
        else
          ret <- pval
      })
      names(res) <- names(Gene_SNP_list)
      
      if (profile_steps)
        step_times$test <- tock(t0_step)
    }
    if (profile_steps) {
      attr(res, "profile") <- list(step_times = step_times)
    }
    return(res)
  }