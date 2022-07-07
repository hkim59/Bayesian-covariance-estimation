### Performance metrics

Ent_loss <- function(Sigma_hat, Sigma_true) {
  p <- nrow(Sigma_true)
  Mat <- Sigma_hat %*% solve(Sigma_true)
  sum(diag(Mat)) - log(det(Mat)) - p
}

Quad_loss <- function(Sigma_hat, Sigma_true) {
  p <- nrow(Sigma_true)
  Mat <- Sigma_hat %*% solve(Sigma_true)
  sum(diag({Mat - diag(p)}**2)) 
}

CV_loss <- function(Sigma1, Sigma2) {
  Mat <- Sigma2 %*% solve(Sigma1)
  -log(det(Sigma1)) -sum(diag(Mat))
}

MCC_index <- function(Sigma_hat, Sigma_true) {
  # TP: \sigma_{ij} = 0, \hat{\sigma}_{ij} = 0
  # TN: \sigma_{ij} \ne 0, \hat{\sigma}_{ij} \ne 0
  # FP: \sigma_{ij} \ne 0, \hat{\sigma}_{ij} = 0
  # FN: \sigma_{ij} = 0, \hat{\sigma}_{ij} \ne 0
  
  TP <- sum((Sigma_true <= 0.03) & (Sigma_hat <= 0.03))
  TN <- sum((Sigma_true > 0.03) & (Sigma_hat > 0.03))
  FP <- sum((Sigma_true > 0.03) & (Sigma_hat <= 0.03))
  FN <- sum((Sigma_true <= 0.03) & (Sigma_hat > 0.03))
  
  num <- (TP * TN) - (FP * FN)
  denom <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
  MCC <- num / sqrt(denom)
  return(MCC)
}