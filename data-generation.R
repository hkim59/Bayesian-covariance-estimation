### Data generation settings. Check the sparsity level of \Sigma_{true} for each (n,p) combinations. 
# Check: Sparsity index, the percentage of zeros on the off-diagonals
# Check: Symmetry & Positive-definiteness.

# Check: Sparsity index, the percentage of zeros on the off-diagonals
# Check: Symmetry & Positive-definiteness.


n <- 100   # sample size
p <- 10    # data dimension
p <- 100
p <- 200
p <- 1000

# 1. Unstructured, random. (Khondker et al. 2013)
alpha <- 0.4
alpha <- 0.15
alpha <- 0.06
alpha <- 0.02
L <- diag(p)
entries <- runif(p*(p-1)/2,-0.5,0.5)
idx <- rbinom(p*(p-1)/2, 1, prob = alpha)
entries[!idx] <- 0
L[lower.tri(L)] <- entries
Uns <- L%*%t(L) # Generate unstructured \Sigma from Cholesky decomposition.
sum(Uns[lower.tri(Uns)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.73


# 2. Structured, AR(1). (Wang 2012)
exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
AR1 <- 0.7^exponent
sum(AR1[lower.tri(AR1)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0. 
sum(eigen(AR1)$values > 0) # Check: Positive-definiteness


# 3. Structured, exponential correlation model (Khondker et al. 2013)
gam <- 1
gam <- 0.5
gam <- 0.3
gam <- 0.1
exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
Exp <- exp(-gam * exponent)
sum(Exp[lower.tri(Exp)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.33
sum(eigen(Exp)$values > 0) # Check: Positive-definiteness

# 4. Structured, block diagonal (Wang 2012, Bien & Tibshirani 2010)
K <- 2
K <- 5
K <- 5
K <- 10
block <- matrix(0.5, nrow = p/K, ncol = p/K)
diag(block) <- 1
Blist <- rep(list(block),K)
Blk <- as.matrix(Matrix::bdiag(Blist))
sum(Blk[lower.tri(Blk)] < 0.01)/(p*(p-1)/2) # Check: Sparsity index = 0.89
sum(eigen(Blk)$values > 0) # Check: Positive-definiteness

# 5. Factor model (Bhattacharya & Dunson 2011)
k <- 3     # factor dimensionality
k <- 10
k <- 10
k <- 15
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
sum(Fct[lower.tri(Fct)] == 0)/(p*(p-1)/2) # Check: Sparsity index = 0.98
sum(eigen(Fct)$values > 0) # Check: Positive-definiteness

###--- STEP 2. Generate the data, Y_i \sim^{iid} N_p(0, \Sigma), i=1,...,n ---###
# Y : n x p matrix of responses
mu <- rep(0, p)
Y_Uns <- MASS::mvrnorm(n=n, mu = mu, Sigma = Uns)
Y_AR1 <- MASS::mvrnorm(n=n, mu = mu, Sigma = AR1)
Y_Exp <- MASS::mvrnorm(n=n, mu = mu, Sigma = Exp)
Y_Blk <- MASS::mvrnorm(n=n, mu = mu, Sigma = Blk)
Y_Fct <- MASS::mvrnorm(n=n, mu = mu, Sigma = Fct)
Y <- Y_Uns

###--- STEP 3. Compare the models ---###
# 1. SCOV
SCOV <- cov(Y)




# 2. SPCOV
P <- matrix(1, p, p)
diag(P) <- 0

# select tuning parameter via 5-fold cross-validation w.r.t entropy loss

lam_vec <- seq(0.01, 1.01, 0.2) # As noted in the SPCOV paper, CV(\lambda) is relatively flat around the maximum. 
idx <- sample(n, n, replace = F)
CV_lam <- vector(length = length(lam_vec))

for (l in 1:length(lam_vec)) {
  lam <- lam_vec[l]
  CV <- vector(length = 5)
  for (i in 1:5) {
    start <- 20*(i-1)+1
    end <- 20*i
    test_idx <- idx[start:end]
    train_idx <- idx[idx != test_idx]
    test_Y <- Y[test_idx,]
    train_Y <- Y[train_idx,]
    test_S <- cov(test_Y)
    train_S <- cov(train_Y)
    mm <- spcov::spcov(Sigma=train_S, # Initial guess for \Sigma
                       S=train_S, # Empirical p.d. covariance matrix.  
                       lambda=lam * P, # Penalty parameter, p x p matrix. Penalizes off-diagonal elements.
                       step.size=100, n.inner.steps=200,
                       thr.inner=0, tol.outer=1e-3, trace=1) # majorize-minimize algorithm tuning parameters
    Sigma_hat <- mm$Sigma
    CV[i] <- CV_loss(Sigma_hat, test_S)
  }
  CV_lam[l] <- mean(CV)
}
lam <- lam_vec[which.max(CV_lam)]

# SPCOV estimate with final tuning parameter. 
mm <- spcov::spcov(Sigma=SCOV, # Initial guess for \Sigma
                   S=SCOV, # Empirical p.d. covariance matrix.  
                   lambda=lam * P, # Penalty parameter, p x p matrix. Penalizes off-diagonal elements.
                   step.size=100, n.inner.steps=200,
                   thr.inner=0, tol.outer=1e-3, trace=1) # majorize-minimize algorithm tuning parameters
SPCOV <- mm$Sigma





# 3. SBIFM
ags <- infinitefactor::linearMGSP(Y, nrun = 10000, burn = 5000, thin = 1)
SBIFM <- ags$omegaSamps # posterior \Sigma samples
#ags$covMean # posterior mean of \Sigma, apply(SBIFM, c(1,2), mean)
numFacts <- ags$numFacts # estimated number of factors

# 4. STEP
# Initial value of parameters
mod_glasso <- glasso::glasso(s=SPCOV,rho=sqrt(log(p)/n))
sigma_inv_lasso <- mod_glasso$wi
fit_StepOmega <- StepOmega(S=SPCOV,temp=sigma_inv_lasso,nmc=1000,burnin=2000,n=n) # one step of the STEPWISE algorithm
STEP <- fit_StepOmega$Omegahat

# 6. BGLasso
bg <- blockGLasso(X=Y,iterations=5000,burnIn=5000)
bg$Sigmas[[1]]
bg$Omegas[[1]]



