### Metrics_store_ftn for storing paramters, given the true \Sigma_{true}.

library(infinitefactor) # SBIFM
library(BayesianGLasso) # BGLasso
library(spcov) # SPCOV
library(glasso) # GLASSO
library(MASS)
library(Matrix)

Metrics_store_ftn <- function(B, p, Sigma_true) {
  
  Uns <- Sigma_true

  ### Store parameters ###
  # Ent_mat 
  Ent_mat <- matrix(NA, nrow = B, ncol = 5)
  # Quad_mat
  Quad_mat <- matrix(NA, nrow = B, ncol = 5)
  Frob_mat <- matrix(NA, nrow = B, ncol = 5)
  # MCC_mat
  MCC_mat <- matrix(NA, nrow = B, ncol = 5)
  # covrge_mat: report for selected indices (5,5) and (7,8); one diagonal term, and one off-diagonal term
  covrge_diag_mat <- matrix(NA, nrow = B, ncol = 2)
  covrge_offdiag_mat <- matrix(NA, nrow = B, ncol = 2)
  # avg_numFacts 
  avg_numFacts <- matrix(NA, nrow = B, ncol = 1)
  com_time <- matrix(NA, nrow = B, ncol = 5)
  
  for (b in 1:B) {
    
    mu <- rep(0, p)
    Y <- MASS::mvrnorm(n=n, mu = mu, Sigma = Uns)
    
    # 1. SCOV
    start <- Sys.time()
    SCOV <- cov(Y)
    end <- Sys.time()
    
    com_time[b,1] <- as.numeric(end - start)
    Frob_mat[b,1] <- Frob_norm(SCOV, Uns)
    Ent_mat[b,1] <- Ent_loss(SCOV, Uns)
    Quad_mat[b,1] <- Quad_loss(SCOV, Uns)
    MCC_mat[b,1] <- MCC_index(SCOV, Uns)
    
    # 2. SPCOV
    # 2. GLasso
    # start <- Sys.time()
    # P <- matrix(1, p, p)
    # diag(P) <- 0
    # # select tuning parameter via 5-fold cross-validation w.r.t entropy loss
    # # SPCOV estimate with final tuning parameter. 
    # lam <- sqrt(log(p)/n)
    # S <- SCOV + 0.001*diag(p)
    # mm <- spcov::spcov(Sigma=S, # Initial guess for \Sigma
    #                    S=SCOV, # Empirical p.d. covariance matrix.  
    #                    lambda=lam * P, # Penalty parameter, p x p matrix. Penalizes off-diagonal elements.
    #                    step.size=100, n.inner.steps=200,
    #                    thr.inner=0, tol.outer=1e-3, trace=1) # majorize-minimize algorithm tuning parameters
    # SPCOV <- mm$Sigma
    # end <- Sys.time()
    
    start <- Sys.time()
    mod_glasso <- glasso::glasso(s=SCOV,rho=sqrt(log(p)/n))
    GLasso <- mod_glasso$wi
    GLasso_inv <- solve(GLasso)
    end <- Sys.time()
    
    com_time[b,2] <- as.numeric(end - start)
    Frob_mat[b,2] <- Frob_norm(GLasso_inv, Uns)
    Ent_mat[b,2] <- Ent_loss(GLasso_inv, Uns)
    Quad_mat[b,2] <- Quad_loss(GLasso_inv, Uns)
    MCC_mat[b,2] <- MCC_index(GLasso_inv, Uns)
    
    # 3. SBIFM
    start <- Sys.time()
    #ags <- infinitefactor::linearMGSP(Y, nrun = 10000, burn = 5000, thin = 1)
    ags <- infinitefactor::linearMGSP(Y, nrun = 2000, burn = 1000, thin = 1)
    SBIFM <- ags$omegaSamps # posterior \Sigma samples
    numFacts <- ags$numFacts # estimated number of factors
    SBIFM_mean <- ags$covMean # posterior mean of \Sigma, apply(SBIFM, c(1,2), mean)
    end <- Sys.time()
    
    com_time[b,3] <- as.numeric(end - start)
    Frob_mat[b,3] <- Frob_norm(SBIFM_mean, Uns)
    Ent_mat[b,3] <- Ent_loss(SBIFM_mean, Uns)
    Quad_mat[b,3] <- Quad_loss(SBIFM_mean, Uns)
    MCC_mat[b,3] <- MCC_index(SBIFM_mean, Uns)
    covrge_diag_mat[b,1] <- {Uns[5,5] >= quantile(SBIFM[5,5,], c(0.025))} & {Uns[5,5] <= quantile(SBIFM[5,5,], c(0.975))}
    covrge_offdiag_mat[b,1] <- {Uns[7,8] >= quantile(SBIFM[7,8,], c(0.025))} & {Uns[7,8] <= quantile(SBIFM[7,8,], c(0.975))}
    avg_numFacts[b] <- mean(numFacts)
    
    # 4. STEP
    # Initial value of parameters
    start <- Sys.time()
    mod_glasso <- glasso::glasso(s=SCOV,rho=sqrt(log(p)/n))
    sigma_inv_lasso <- mod_glasso$wi
    fit_StepOmega <- StepOmega(S=SCOV,temp=sigma_inv_lasso,nmc=1000,burnin=2000,n=n) # one step of the STEPWISE algorithm
    STEP <- fit_StepOmega$Omegahat
    STEP_inv <- solve(STEP)
    end <- Sys.time()
    
    com_time[b,4] <- as.numeric(end - start)
    Frob_mat[b,4] <- Frob_norm(STEP_inv, Uns)
    Ent_mat[b,4] <- Ent_loss(STEP_inv, Uns)
    Quad_mat[b,4] <- Quad_loss(STEP_inv, Uns)
    MCC_mat[b,4] <- MCC_index(STEP_inv, Uns)
    
    # 5. BGLasso
    start <- Sys.time()
    #bg <- blockGLasso(X=Y,iterations=5000,burnIn=5000)
    # bg <- blockGLasso(X=Y,iterations=1000,burnIn=100)
    # BGLasso <- bg$Sigmas
    # BGLasso_mean <- Reduce("+", BGLasso) / length(BGLasso)
    # end <- Sys.time()
    # 
    # niter = 1000
    # com_time[b,5] <- as.numeric(end - start)
    # Frob_mat[b,5] <- Frob_norm(BGLasso_mean, Uns)
    # Ent_mat[b,5] <- Ent_loss(BGLasso_mean, Uns)
    # Quad_mat[b,5] <- Quad_loss(BGLasso_mean, Uns)
    # MCC_mat[b,5] <- MCC_index(BGLasso_mean, Uns)
    # covrge_diag_mat[b,2] <- {Uns[5,5] >= quantile(sapply(1:niter, function(i) {BGLasso[[i]][5,5]}), 0.025)} & {Uns[5,5] <= quantile(sapply(1:niter, function(i) {BGLasso[[i]][5,5]}), 0.975)}
    # covrge_offdiag_mat[b,2] <- {Uns[7,8] >= quantile(sapply(1:niter, function(i) {BGLasso[[i]][7,8]}), 0.025)} & {Uns[7,8] <= quantile(sapply(1:niter, function(i) {BGLasso[[i]][7,8]}), 0.975)}

  }
  
  
  return(list(
    com_time = com_time,
    Frob_mat = Frob_mat,
    Ent_mat = Ent_mat, 
    Quad_mat = Quad_mat,
    MCC_mat = MCC_mat,
    covrge_diag_mat = covrge_diag_mat,
    covrge_offdiag_mat = covrge_offdiag_mat,
    avg_numFacts = avg_numFacts
  ))
  
}