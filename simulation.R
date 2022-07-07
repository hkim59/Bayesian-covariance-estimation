###--------------------------------------------------------------------------------------------###
### A Simulation study: Bayesian approaches to structured covariance estimation                ###
### Author: Hyoshin Kim                                                                        ###
### Date: July 5th, 2022                                                                       ###
### Goal: Covariance matrix estimation performance comparison under structure misspecification ###
###--------------------------------------------------------------------------------------------###

library(Rcpp)
library(RcppArmadillo)

# JRNS 
sourceCpp("JRNS/StepOmega.cpp")

source("utils.R") # Load performance evaluation metrics
source("code.R") # Load Metrics_store_ftn for storing paramters

### Models to compare: 
# SCOV: Basic sample covariance
# SPCOV: A Frequentist approach. Sparse Estimation of a Covariance Matrix, Bien & Tibshirani 2012. https://statweb.stanford.edu/~jbien/biometrika2011spcov.pdf
# SBIFM: Sparse Bayesian infinite factor models, Bhattacharya & Dunson 2011, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3419391/pdf/asr013.pdf
# STEP: Precision matrix estimation. Simultaneous Regression and Covariance Selection, Samanta et al. 2022. https://arxiv.org/abs/2201.05653
# BGLasso: Precision matrix estimation. Bayesian Graphical Lasso, Wang 2012. https://projecteuclid.org/journals/bayesian-analysis/volume-7/issue-4/Bayesian-Graphical-Lasso-Models-and-Efficient-Posterior-Computation/10.1214/12-BA729.full

### Sample size n, dimensionality p: 
# moderate n, low p, p < n: n = 100, p = 10
# moderate n, moderate p, p > n: n = 100, p = 100
# moderate n, large p, p > n: n = 100, p = 200
# moderate n, very large p, p >> n: n = 100, p = 1000

### Consider independent samples with a common covariance structure: Y_i \sim^{iid} N_p(0, \Sigma), i=1,...,n.
### For each n,p generate the true covariance \Sigma_{true} from the following.

# Unstructured, random: Independently assign \Sigma_{ij} = \Sigma_{ji} to be non-zero with some probability. Generate \Sigma from the Cholesky decomposition.
# Structured, AR(1): Assign \sigma_{ij} = 0.7^{|i-j|} 
# Structured, exponential correlation model: Assign \sigma_{ij} = exp(-\alpha |i - j|), i < j.
# Structured, Block model: \Sigma = diag(\Sigma_1, ..., \Sigma_5). \Sigma_k are dense. For each \Sigma_k, \sigma_{ii} = 1, \sigma_{ij} = 0.5. 
# Structured, Factor model: \Sigma = \Lambda\Lambda^{\top} + D. D: diagonal, independently sample entries of D^{-1} from Gamma(shape = 1, rate = 0.25) with mean 4.


### Note. The JRNS method estimates the precision matrix. 
###       AR(1) covariance matrix => Tridiagonal precision matrix
###       Tridigonal covariance matrix => AR(1) precision matrix
###       Block diagonal covariance matrix => Block diagonal precision matrix
###       Factor models: Assume a lower dimensional structure, k < p
###       We may estimate the precision matrix => Invert to obtain the covariance matrix estimate. 

### Data dimension tuple
n <- 100
p_list <- c(10, 100, 200, 1000)

### Tuning parameters for generating the true covariance structure \Sigma_{true}, for each p=10,100,200,1000
alpha_list <- c(0.4, 0.15, 0.06, 0.02) # Uns
gam_list <- c(1, 0.5, 0.3, 0.1)        # Exp
K_list <- c(2,5,5,10)                  # Blk
k_list <- c(3,10,10,15)                # Fct, factor dimensionality


for (id in 1:length(p_list)) {
  
  p <- p_list[id]
  alpha <- alpha_list[id]
  gam <- gam_list[id]
  K <- K_list[id]
  k <- k_list[id]
  
  set.seed(123456)
  
  ###--- STEP 1. Generate the true covariance, \Sigma_{true}---###
  # 1. Unstructured, random. (Khondker et al. 2013)
  L <- diag(p)
  entries <- runif(p*(p-1)/2,-0.5,0.5)
  idx <- rbinom(p*(p-1)/2, 1, prob = alpha)
  entries[!idx] <- 0
  L[lower.tri(L)] <- entries
  Uns <- L%*%t(L) # Generate unstructured \Sigma from Cholesky decomposition.
  #sum(Uns[lower.tri(Uns)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.73
  
  # 2. Structured, AR(1). (Wang 2012)
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
  AR1 <- 0.7^exponent
  #sum(AR1[lower.tri(AR1)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0. 
  #sum(eigen(AR1)$values > 0) # Check: Positive-definiteness
  
  # 3. Structured, exponential correlation model (Khondker et al. 2013)
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
  Exp <- exp(-gam * exponent)
  #sum(Exp[lower.tri(Exp)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.33
  #sum(eigen(Exp)$values > 0) # Check: Positive-definiteness
  
  # 4. Structured, block diagonal (Wang 2012, Bien & Tibshirani 2010)
  block <- matrix(0.5, nrow = p/K, ncol = p/K)
  diag(block) <- 1
  Blist <- rep(list(block),K)
  Blk <- as.matrix(Matrix::bdiag(Blist))
  #sum(Blk[lower.tri(Blk)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.89
  #sum(eigen(Blk)$values > 0) # Check: Positive-definiteness
  
  # 5. Factor model (Bhattacharya & Dunson 2011)
  D <- matrix(0, nrow = p, ncol = p)
  entries <- rgamma(p, shape = 1, rate = 0.25)
  diag(D) <- entries
  
  Lambda <- matrix(0, nrow = p, ncol = k)
  nonzero_num <- seq(2*k, k+1, length.out = k)
  for (j in 1:k) {
    num <- nonzero_num[j]
    entries <- rnorm(num, mean = 0, sd = 1)
    idx <- sample(p, size = num, replace = F)
    Lambda[idx,j] <- entries
  }
  Fct <- Lambda %*% t(Lambda) + D
  #sum(Fct[lower.tri(Fct)] == 0)/(p*(p-1)/2) # Check: Sparsity index = 0.98
  #sum(eigen(Fct)$values > 0) # Check: Positive-definiteness
  
  
  
  ###--- STEP 2. Store the performance metrics ---###
  
  # Total number of Monte Carlo samples
  B <- 50
  set.seed(123456)
  start <- Sys.time()
  para_Fct <- Metrics_store_ftn(B, p, Sigma_true = Fct)
  
  end <- Sys.time()
  time <- end - start
  
  # para_Uns <- Metrics_store_ftn(B, p, Sigma_true = Uns)
  # para_AR1 <- Metrics_store_ftn(B, p, Sigma_true = AR1)
  # para_Exp <- Metrics_store_ftn(B, p, Sigma_true = Exp)
  # para_Blk <- Metrics_store_ftn(B, p, Sigma_true = Blk)
  # para_Fct <- Metrics_store_ftn(B, p, Sigma_true = Fct)

  
}


###--- Generate summary tables ---###
para <- para_Fct
round(colMeans(para$com_time),2)
round(colMeans(para$Ent_mat),2)
round(colMeans(para$Quad_mat),2)
round(colMeans(para$Frob_mat),2)
round(colMeans(para$MCC_mat),2)
round(colMeans(para$covrge_diag_mat),2)
round(colMeans(para$covrge_offdiag_mat),2)
round(colMeans(para$avg_numFacts),2)
















