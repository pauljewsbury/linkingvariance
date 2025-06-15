
### supplement for "Standard error estimation for subpopulation non-invariance"
### accepted for publication in Applied Psychological Measurement in 2025

### questions can be directed to pjewsbury@gmail.com; pjewsbury@ets.org

library(countrycode)
library(dplyr)
library(haven)
library(tidyr)
library(ggplot2)

A <- c(538.07815,537.13672,536.30970,536.57072,536.93362)
B <- c(100.74315,101.42598,100.88021,101.72353,101.25662)
M <- 5 # number of PVs

weighted.sd <- function(x, w) sqrt(sum(w * (x - sum(w * x) / sum(w))^2) / (sum(w) - 1))
weighted.sd.p <- function(x, w) sqrt(sum(w * (x - sum(w * x) / sum(w))^2) / sum(w))
colWeightedMeans <- function(x, w) apply(x,2,function(k) weighted.mean(k, w))
colWeightedSDs <- function(x, w) apply(x,2,function(k) weighted.sd(k, w))
weighted.scale <- function(x, w) (x - weighted.mean(x, w))/weighted.sd.p(x, w)

weighted.msd <- function(x, w) weighted.mean(x, w) + weighted.sd(x,w)

constructJK <- function(JKREP,JKZONE,TOTWGT)
{
  #max_zone <- max(JKZONE)
  max_zone <- 125
  JKweights <- matrix(nrow=length(TOTWGT),ncol=(max_zone*2))
  
  for (i in 1:max_zone)
    for (j in 0:1)
      JKweights[,(i-1)*2 + j+1] <- TOTWGT*((i!=JKZONE)*1 + 
                                         (i==JKZONE)*(j==JKREP)*2 + 
                                         (i==JKZONE)*(j!=JKREP)*0)
  
  colnames(JKweights) <- sprintf("JKweights_%02d", 1:(max_zone * 2))
  JKweights
}

# apply fun for each level of group, return point estimate and BI variance
analyze_PVs <- function(mat, fun, group, weights) {
  if (length(unique(group)) == 1) {
    estimates <- apply(mat, 2, function(col) fun(col,weights))
    c(est = mean(estimates), BI = var(estimates))
  }
  else
  {
    estimates <- apply(mat, 2, function(col) tapply(seq_along(col), group, function(idx) fun(col[idx],weights[idx])))
    t(apply(estimates, 1, function(x) c(est = mean(x), BI = var(x))))
  }
}

# apply fun for each level of group, return point estimate and BI variance
# special function for standardized score (Z-score)
std_analyze_PVs <- function(mat, fun, group, weights) {
  if (length(unique(group)) == 1) {
    c(est = 0, BI = 0)
  }
  else
  {
    std_mat <- apply(mat, 2, function(x) weighted.scale(x, weights))
    estimates <- apply(std_mat, 2, function(col) tapply(seq_along(col), group, function(idx) fun(col[idx],weights[idx])))
    t(apply(estimates, 1, function(x) c(est = mean(x), BI = var(x))))
  }
}


# get WI variance
getWI <- function(scores, fun, group, weight, JKweights) {
  nrep <- ncol(JKweights)
  
  #orig_est <- tapply(seq_along(scores), group, function(idx) fun(scores[idx], weight[idx]))
  
  replicate_est <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(scores), group, function(idx) fun(scores[idx], JKweights[idx, r]))
  }))
  

  if (length(unique(group)) > 1) {
    orig_est <- colMeans(replicate_est)
    WI <- (1/2) * colSums((sweep(replicate_est, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est)
    WI <- (1/2) * sum((replicate_est - as.numeric(orig_est))^2)
  }
  
  WI
}

# average across WI estimate for each set of PVs
getWImat <- function(mat, fun, group, weight, JKweights) {
  if (length(unique(group)) > 1)
    rowMeans(sapply(1:ncol(mat), function(j) getWI(mat[, j], fun, group, weight, JKweights)))
  else
    mean(sapply(1:ncol(mat), function(j) getWI(mat[, j], fun, group, weight, JKweights)))
}

# get WI variance
# special function for Z-scores
std_getWI <- function(scores, fun, group, weight, JKweights) {
  nrep <- ncol(JKweights)
  
  if (length(unique(group)) == 1) return(0)
  
  replicate_est <- t(sapply(1:nrep, function(r) {
    std_scores <- weighted.scale(scores, JKweights[,r])
    tapply(seq_along(std_scores), group, function(idx) fun(std_scores[idx], JKweights[idx, r]))
  }))

  orig_est <- colMeans(replicate_est)
  (1/2) * colSums((sweep(replicate_est, 2, orig_est, "-"))^2)
}

# average across WI estimate for each set of PVs
std_getWImat <- function(mat, fun, group, weight, JKweights) {
  if (length(unique(group)) == 1) return(0)
  
  rowMeans(sapply(1:ncol(mat), function(j) std_getWI(mat[, j], fun, group, weight, JKweights)))
}

# links R to A, returns transformed R; matrix version
# m_R and m_A control which PVs (or all in case of 0) to use to calculate linking coefficients
linkmat <- function(mat_R,mat_A,m_R,m_A,weight_R,weight_A) {
  mu_theta <- if (m_R > 0) weighted.mean(mat_R[,m_R],weight_R)
                 else mean(apply(mat_R, 2, weighted.mean, w = weight_R))
  mu_X     <- if (m_A > 0) weighted.mean(mat_A[,m_A],weight_A)
                 else mean(apply(mat_A, 2, weighted.mean, w = weight_A))
  
  sigma_theta <- if (m_R > 0) weighted.sd(mat_R[,m_R],weight_R)
                 else mean(apply(mat_R, 2, weighted.sd, w = weight_R))
  sigma_X     <- if (m_A > 0) weighted.sd(mat_A[,m_A],weight_A)
                 else mean(apply(mat_A, 2, weighted.sd, w = weight_A))
  
  A_n  <- sigma_X / sigma_theta
  B_n  <- mu_X - A_n * mu_theta
  
  A_n * mat_R + B_n
}


# special code to get covariance between group standardized means and a sample stat
std_getCovariance <- function(mat, fun, group, weights, JKweights)
{
  
  # get within-imputation variance
  nrep <- ncol(JKweights)
  M <- ncol(mat)
  
  WI <- rowMeans(sapply(1:M,function(m) {
    
    scores <- mat[,m]
    
    std_replicate_est <- t(sapply(1:nrep, function(r) {
      std_scores <- weighted.scale(scores, JKweights[,r])
      tapply(seq_along(std_scores), group, function(idx) weighted.mean(std_scores[idx], JKweights[idx, r]))
    }))
    
    pop_replicate_est <- t(sapply(1:nrep, function(r) {
      fun(scores, JKweights[, r])
    }))
    
    std_orig_est <- colMeans(std_replicate_est)
    pop_orig_est <- mean(pop_replicate_est)
    pop_matrix <- matrix(pop_replicate_est, nrow = nrep, ncol = ncol(std_replicate_est))
    
    (1/2) * colSums((sweep(std_replicate_est, 2, std_orig_est, "-")) * (pop_matrix - pop_orig_est))
  }))
  
  
  # get between-imputation variance
  std_mat <- apply(mat, 2, function(x) weighted.scale(x, weights))
  std_estimates <- apply(std_mat, 2, function(col) tapply(seq_along(col), group, function(idx) weighted.mean(col[idx],weights[idx])))
  pop_estimates <- apply(mat, 2, function(col) fun(col, weights))
  BI <- t(apply(std_estimates, 1, function(x) cov(x,pop_estimates)))
  
  WI + (1 + 1/M) * BI
}


# links R to A, returns transformed R
link <- function(PV_R,PV_A,weight_R,weight_A) {
  A_n  <- weighted.sd(PV_A,weight_A) / weighted.sd(PV_R,weight_R)
  B_n  <- weighted.mean(PV_A,weight_A) - A_n * weighted.mean(PV_R,weight_R)
  PV_R <- A_n * PV_R + B_n
}


# R - A
# need to give it untransformed R
gettotalWI <- function(mat_R, mat_A, m_R, m_A, fun, group_R, group_A, 
                           weight_R, weight_A, JKweights_R, JKweights_A)
{
  scores_R <- mat_R[,m_R]
  scores_A <- mat_A[,m_A]
  
  # combine scores to simplify later code
  combined_scores <- c(scores_R,scores_A)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(scores_R)),rep(2,length(scores_A)))
  
  nrep <- ncol(JKweights_R)
  
  replicate_est <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(combined_scores), combined_group, function(idx) {
      trans_R <- linkmat(mat_R,mat_A,m_R,m_A,JKweights_R[,r],JKweights_A[,r])[,m_R]
      idx_R <- idx[sample_ind[idx] == 1]
      idx_A <- idx[sample_ind[idx] == 2] - length(scores_R)
      
      (fun(trans_R[idx_R],JKweights_R[idx_R,r]) -
          fun(scores_A[idx_A],JKweights_A[idx_A,r]))
    })}))
  
  if (length(unique(group_R)) > 1) {
    orig_est <- colMeans(replicate_est)
    WI_R <- (1/2) * colSums((sweep(replicate_est, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est)
    WI_R <- (1/2) * sum((replicate_est - as.numeric(orig_est))^2)
  }
}

# average across WI estimate for each set of PVs
gettotalWImat <- function(mat_R, mat_A, fun, group_R, group_A, 
                              weight_R, weight_A, JKweights_R, JKweights_A) {
  M <- ncol(mat_R)
  
  if (length(unique(group_R)) > 1)
    rowMeans(sapply(1:M,function(m) gettotalWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                                   weight_R, weight_A, JKweights_R, JKweights_A)))
  else
    mean(sapply(1:M,function(m) gettotalWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                               weight_R, weight_A, JKweights_R, JKweights_A)))
}

# R - A
# need untransformed R
gettotalBI <- function(mat_R, mat_A, fun, group_R, group_A, weight_R, weight_A)
{
  M <- ncol(mat_R)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(group_R)),rep(2,length(group_A)))
  
  # BI for R
  BI <- apply(sapply(1:M, function(m) tapply(seq_along(combined_group), combined_group, function(idx) {
    trans_R <- linkmat(mat_R, mat_A, m, m, weight_R, weight_A)
    idx_R <- idx[sample_ind[idx] == 1]
    idx_A <- idx[sample_ind[idx] == 2] - length(group_R)
    
    fun(trans_R[idx_R,m],weight_R[idx_R]) - fun(mat_A[idx_A,m],weight_A[idx_A])
  })),1,var)
  
  BI
}

# R - A
# need to give it untransformed R
getcomponentWI <- function(mat_R, mat_A, m_R, m_A, fun, group_R, group_A, 
                        weight_R, weight_A, JKweights_R, JKweights_A)
{
  scores_R <- mat_R[,m_R]
  scores_A <- mat_A[,m_A]

  # combine scores to simplify later code
  combined_scores <- c(scores_R,scores_A)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(scores_R)),rep(2,length(scores_A)))
  
  # WI for R
  nrep <- ncol(JKweights_R)
  
  replicate_est_R <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(combined_scores), combined_group, function(idx) {
      idx_R <- idx[sample_ind[idx] == 1]
      idx_A <- idx[sample_ind[idx] == 2] - length(scores_R)
      
      mean(sapply(1:M, function(i) {
        trans_R <- link(mat_R[,m_R],mat_A[,i],JKweights_R[,r],weight_A)
        fun(trans_R[idx_R],JKweights_R[idx_R,r]) - fun(mat_A[idx_A,i],weight_A[idx_A])
      })) 
    })}))
  
  if (length(unique(group_R)) > 1) {
    orig_est <- colMeans(replicate_est_R)
    WI_R <- (1/2) * colSums((sweep(replicate_est_R, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est_R)
    WI_R <- (1/2) * sum((replicate_est_R - as.numeric(orig_est))^2)
  }
  
  # WI for A
  nrep <- ncol(JKweights_A)
  
  replicate_est_A <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(combined_scores), combined_group, function(idx) {
      idx_R <- idx[sample_ind[idx] == 1]
      idx_A <- idx[sample_ind[idx] == 2] - length(scores_R)
      
      mean(sapply(1:M, function(i) {
        trans_R <- link(mat_R[,i],mat_A[,m_A],weight_R,JKweights_A[,r])
        fun(trans_R[idx_R],weight_R[idx_R]) - fun(mat_A[,m_A][idx_A],JKweights_A[idx_A,r])
      }))
      
    })}))
  
  if (length(unique(group_R)) > 1) {
    orig_est <- colMeans(replicate_est_A)
    WI_A <- (1/2) * colSums((sweep(replicate_est_A, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est_A)
    WI_A <- (1/2) * sum((replicate_est_A - as.numeric(orig_est))^2)
  }
  
  WI_R + WI_A
}

# average across WI estimate for each set of PVs
getcomponentWImat <- function(mat_R, mat_A, fun, group_R, group_A, 
                              weight_R, weight_A, JKweights_R, JKweights_A) {
  M <- ncol(mat_R)
  
  if (length(unique(group_R)) > 1)
    rowMeans(sapply(1:M,function(m) getcomponentWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                             weight_R, weight_A, JKweights_R, JKweights_A)))
  else
    mean(sapply(1:M,function(m) getcomponentWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                                   weight_R, weight_A, JKweights_R, JKweights_A)))
}

# R - A
# need untransformed R
getcomponentBI <- function(mat_R, mat_A, fun, group_R, group_A, weight_R, weight_A)
{
  M <- ncol(mat_R)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(group_R)),rep(2,length(group_A)))
  
  # BI for R
  BI_R <- apply(sapply(1:M, function(m) tapply(seq_along(combined_group), combined_group, function(idx) {
    idx_R <- idx[sample_ind[idx] == 1]
    idx_A <- idx[sample_ind[idx] == 2] - length(group_R)
    
    mean(sapply(1:M, function(i) {
      trans_R <- link(mat_R[,m], mat_A[,i], weight_R, weight_A)
      fun(trans_R[idx_R],weight_R[idx_R]) -  fun(mat_A[idx_A,i],weight_A[idx_A])
    }))
  })),1,var)
  
  # BI for A
  BI_A <- apply(sapply(1:M, function(m) tapply(seq_along(combined_group), combined_group, function(idx) {
    idx_R <- idx[sample_ind[idx] == 1]
    idx_A <- idx[sample_ind[idx] == 2] - length(group_R)
    
    mean(sapply(1:M, function(i) {
      trans_R <- link(mat_R[,i], mat_A[,m], weight_R, weight_A)
      fun(trans_R[idx_R],weight_R[idx_R])
    })) - fun(mat_A[idx_A,m],weight_A[idx_A])
  })),1,var)
  
  BI_R + BI_A
}


# linking component
getlinkWI <- function(mat_R, mat_A, m_R, m_A, fun, group_R, group_A, 
                           weight_R, weight_A, JKweights_R, JKweights_A)
{
  scores_R <- mat_R[,m_R]
  scores_A <- mat_A[,m_A]
  
  # combine scores to simplify later code
  combined_scores <- c(scores_R,scores_A)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(scores_R)),rep(2,length(scores_A)))
  
  # WI for link-R
  nrep <- ncol(JKweights_R)
  
  replicate_est_R <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(combined_scores), combined_group, function(idx) {
      trans_R <- linkmat(mat_R,mat_A,m_R,0,JKweights_R[,r],weight_A)[,m_R]
      idx_R <- idx[sample_ind[idx] == 1]
      idx_A <- idx[sample_ind[idx] == 2] - length(scores_R)
      
      (fun(trans_R[idx_R],weight_R[idx_R]) -
          fun(scores_A[idx_A],weight_A[idx_A]))
    })}))
  
  if (length(unique(group_R)) > 1) {
    orig_est <- colMeans(replicate_est_R)
    WI_R <- (1/2) * colSums((sweep(replicate_est_R, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est_R)
    WI_R <- (1/2) * sum((replicate_est_R - as.numeric(orig_est))^2)
  }
  
  # WI for link-A
  nrep <- ncol(JKweights_A)
  
  replicate_est_A <- t(sapply(1:nrep, function(r) {
    tapply(seq_along(combined_scores), combined_group, function(idx) {
      trans_R <- linkmat(mat_R,mat_A,0,m_A,weight_R,JKweights_A[,r])[,m_A]
      idx_R <- idx[sample_ind[idx] == 1]
      idx_A <- idx[sample_ind[idx] == 2] - length(scores_R)
      
      (fun(trans_R[idx_R],weight_R[idx_R]) -
          fun(scores_A[idx_A],weight_A[idx_A]))
    })}))
  
  if (length(unique(group_R)) > 1) {
    orig_est <- colMeans(replicate_est_A)
    WI_A <- (1/2) * colSums((sweep(replicate_est_A, 2, orig_est, "-"))^2)
  }
  else {
    orig_est <- mean(replicate_est_A)
    WI_A <- (1/2) * sum((replicate_est_A - as.numeric(orig_est))^2)
  }
  
  WI_R + WI_A
}

# average across WI estimate for each set of PVs
getlinkWImat <- function(mat_R, mat_A, fun, group_R, group_A, 
                              weight_R, weight_A, JKweights_R, JKweights_A) {
  M <- ncol(mat_R)
  
  if (length(unique(group_R)) > 1)
    rowMeans(sapply(1:M,function(m) getlinkWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                                   weight_R, weight_A, JKweights_R, JKweights_A)))
  else
    mean(sapply(1:M,function(m) getlinkWI(mat_R, mat_A, m, m, fun, group_R, group_A,
                                               weight_R, weight_A, JKweights_R, JKweights_A)))
}

# get link BI
getlinkBI <- function(mat_R, mat_A, fun, group_R, group_A, weight_R, weight_A)
{
  M <- ncol(mat_R)
  combined_group <- c(group_R,group_A)
  sample_ind <- c(rep(1,length(group_R)),rep(2,length(group_A)))
  
  # BI for link-R
  BI_R <- apply(sapply(1:M, function(m) tapply(seq_along(combined_group), combined_group, function(idx) {
    trans_R <- linkmat(mat_R, mat_A, m, 0, weight_R, weight_A)
    idx_R <- idx[sample_ind[idx] == 1]
    idx_A <- idx[sample_ind[idx] == 2] - length(group_R)
    
    analyze_PVs(trans_R, fun, group_R, weight_R)[1]
  })),1,var)
  
  # BI for link-A
  BI_A <- apply(sapply(1:M, function(m) tapply(seq_along(combined_group), combined_group, function(idx) {
    trans_R <- linkmat(mat_R, mat_A, 0, m, weight_R, weight_A)
    idx_R <- idx[sample_ind[idx] == 1]
    idx_A <- idx[sample_ind[idx] == 2] - length(group_R)
    
    analyze_PVs(trans_R, fun, group_R, weight_R)[1]
  })),1,var)

  
  BI_R + BI_A
}

# get conventional equation. receives transformed R
getconv_eq <- function(mat_R, mat_A, group_R, group_A, 
                       weight_R, weight_A, JKweights_R, JKweights_A)
{
  M <- ncol(mat_R)
  pop_R <- rep(1,length(weight_R))
  pop_A <- rep(1,length(weight_A))
  
  # get point estimates and BI
  PV_sigma_X <- analyze_PVs(mat_A, weighted.sd, pop_A, weight_A)
  PV_mu_X    <- analyze_PVs(mat_A, weighted.mean, pop_A, weight_A)
  
  PV_sigma_Y <- analyze_PVs(mat_R, weighted.sd, pop_R, weight_R)
  PV_mu_Y    <- analyze_PVs(mat_R, weighted.mean, pop_R, weight_R)
  
  PV_Z_X     <- std_analyze_PVs(mat_A, weighted.mean, group_A, weight_A)
  PV_Z_Y     <- std_analyze_PVs(mat_R, weighted.mean, group_R, weight_R)
  
  # point estimates
  sigma_X   <- PV_sigma_X[1]
  Z_X       <- PV_Z_X[,1]
  Z_Y       <- PV_Z_Y[,1]
  
  # total variances
  var_sigma_X <- (getWImat(mat_A, weighted.sd, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_sigma_X[2])
  var_mu_X <- (getWImat(mat_A, weighted.mean, pop_A, weight_A, JKweights_A)
               + (1+1/M)*PV_mu_X[2])
  
  var_sigma_Y <- (getWImat(mat_R, weighted.sd, pop_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_sigma_Y[2])
  var_mu_Y <- (getWImat(mat_R, weighted.mean, pop_R, weight_R, JKweights_R)
               + (1+1/M)*PV_mu_Y[2])  
  
  var_Z_X     <- (std_getWImat(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_Z_X[,2])
  var_Z_Y     <- (std_getWImat(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_Z_Y[,2])
  
  
  # covariances
  var_mu_sigma_X <- (getWImat(mat_A, weighted.msd, pop_A, weight_A, JKweights_A)
                     + (1+1/M)*analyze_PVs(mat_A, weighted.msd, pop_A, weight_A)[2])
  var_mu_sigma_Y <- (getWImat(mat_R, weighted.msd, pop_R, weight_R, JKweights_R)
                     + (1+1/M)*analyze_PVs(mat_R, weighted.msd, pop_R, weight_R)[2])
  cov_mu_sigma_X <- .5 * (var_mu_sigma_X - var_mu_X - var_sigma_X)
  cov_mu_sigma_Y <- .5 * (var_mu_sigma_Y - var_mu_Y - var_sigma_Y)
  
  cov_mu_Z_X    <- std_getCovariance(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
  cov_sigma_Z_X <- std_getCovariance(mat_A, weighted.sd, group_A, weight_A, JKweights_A)
  
  cov_mu_Z_Y    <- std_getCovariance(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
  cov_sigma_Z_Y <- std_getCovariance(mat_R, weighted.sd, group_R, weight_R, JKweights_R)
  
  
  (sigma_X^2 * var_Z_X + Z_X^2 * var_sigma_X + var_mu_X + sigma_X^2 * var_Z_Y + Z_Y^2 * var_sigma_Y + var_mu_Y
    + var_sigma_X * var_Z_X + var_sigma_Y * var_Z_Y
    + 2 * Z_X * cov_mu_sigma_X + 2 * Z_Y * cov_mu_sigma_Y
    + 2 * sigma_X * Z_X * cov_sigma_Z_X + 2 * sigma_X * cov_mu_Z_X
    + 2 * sigma_X * Z_Y * cov_sigma_Z_Y + 2 * sigma_X * cov_mu_Z_Y
  )
}

# get total equation
gettotal_eq <- function(mat_R, mat_A, group_R, group_A, 
                        weight_R, weight_A, JKweights_R, JKweights_A)
{
  # X is A. \theta is untransformed R. Y is transformed R
  # (var(sigma_X) + sigma_X^2)(var(Z_\theta) + var(Z_X))
  #   + (Z_\theta^2 + Z_X^2 -2 Z_\theta Z_\X) var(\sigma_X)
  
  M <- ncol(mat_R)
  pop_R <- rep(1,length(weight_R))
  pop_A <- rep(1,length(weight_A))
  
  # get point estimates and BI
  PV_sigma_X <- analyze_PVs(mat_A, weighted.sd, pop_A, weight_A)
  PV_Z_X     <- std_analyze_PVs(mat_A, weighted.mean, group_A, weight_A)
  PV_Z_Y     <- std_analyze_PVs(mat_R, weighted.mean, group_R, weight_R)
  
  # point estimates
  sigma_X   <- PV_sigma_X[1]
  Z_X       <- PV_Z_X[,1]
  Z_Y       <- PV_Z_Y[,1]
  
  # total variances
  var_sigma_X <- (getWImat(mat_A, weighted.sd, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_sigma_X[2])
  var_Z_X     <- (std_getWImat(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_Z_X[,2])
  var_Z_Y     <- (std_getWImat(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_Z_Y[,2])
  
  (var_sigma_X + sigma_X^2) * (var_Z_Y + var_Z_X) + (Z_X^2 + Z_Y^2 - 2 * Z_X * Z_Y) * var_sigma_X
}


# get component equation
getcomponent_eq <- function(mat_R, mat_A, group_R, group_A,
                            weight_R, weight_A, JKweights_R, JKweights_A)
{
  # X is A. \theta is untransformed R. Y is transformed R
  # (var(\sigma_X) + \sigma_X^2) var(Z_X) + (Z_\theta^2 + Z_X^2 -2 Z_\theta Z_\X) var(\sigma_X)
  # sigma_X^2 var(Z_\theta)
  
  M <- ncol(mat_R)
  pop_R <- rep(1,length(weight_R))
  pop_A <- rep(1,length(weight_A))
  
  # get point estimates and BI
  PV_sigma_X <- analyze_PVs(mat_A, weighted.sd, pop_A, weight_A)
  PV_Z_X     <- std_analyze_PVs(mat_A, weighted.mean, group_A, weight_A)
  PV_Z_Y     <- std_analyze_PVs(mat_R, weighted.mean, group_R, weight_R)
  
  # point estimates
  sigma_X   <- PV_sigma_X[1]
  Z_X       <- PV_Z_X[,1]
  Z_Y       <- PV_Z_Y[,1]
  
  # total variances
  var_sigma_X <- (getWImat(mat_A, weighted.sd, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_sigma_X[2])
  var_Z_X     <- (std_getWImat(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_Z_X[,2])
  var_Z_Y     <- (std_getWImat(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_Z_Y[,2])
  
  compX <- (var_sigma_X + sigma_X^2) * var_Z_X + (Z_Y^2 + Z_X^2 - 2 * Z_Y * Z_X) * var_sigma_X
  compY <- sigma_X^2 * var_Z_Y
  
  compX + compY
}


# takes transformed mat_R. takes transformed R
getCTcoefficients <- function(mat_R, mat_A, weight_R, weight_A, JKweights_R, JKweights_A)
{
  # X is A. \theta is untransformed R. Y is transformed R
  M <- ncol(mat_R)
  pop_R <- rep(1,length(weight_R))
  pop_A <- rep(1,length(weight_A))
  
  # get point estimates and BI
  PV_mu_Y <- analyze_PVs(mat_R, weighted.mean, pop_R, weight_R)
  PV_mu_X <- analyze_PVs(mat_A, weighted.mean, pop_A, weight_A)
  PV_sigma_Y <- analyze_PVs(mat_R, weighted.sd, pop_R, weight_R)
  PV_sigma_X <- analyze_PVs(mat_A, weighted.sd, pop_A, weight_A)
  
  # point estimates
  mu_X <- PV_mu_X[1]
  sigma_X <- PV_sigma_X[1]
  
  # total variances
  var_mu_X <- (getWImat(mat_A, weighted.mean, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_mu_X[2])
  var_mu_Y <- (getWImat(mat_R, weighted.mean, pop_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_mu_Y[2])
  var_sigma_X <- (getWImat(mat_A, weighted.sd, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_sigma_X[2])
  var_sigma_Y <- (getWImat(mat_R, weighted.sd, pop_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_sigma_Y[2])
  
  d1 <- var_sigma_X - var_sigma_Y
  d2 <- -2 * var_sigma_X
  
  var_mu_sigma_X <- (getWImat(mat_A, weighted.msd, pop_A, weight_A, JKweights_A)
                    + (1+1/M)*analyze_PVs(mat_A, weighted.msd, pop_A, weight_A)[2])
  var_mu_sigma_Y <- (getWImat(mat_R, weighted.msd, pop_R, weight_R, JKweights_R)
                    + (1+1/M)*analyze_PVs(mat_R, weighted.msd, pop_R, weight_R)[2])
  cov_mu_sigma_X <- .5 * (var_mu_sigma_X - var_mu_X - var_sigma_X)
  cov_mu_sigma_Y <- .5 * (var_mu_sigma_Y - var_mu_Y - var_sigma_Y)
  
  d3 <- 2 / sigma_X * (-1 * cov_mu_sigma_Y)
  d4 <- 2 / sigma_X * (-1 * cov_mu_sigma_X)
  
  # return coefficients, c(AY, AXY, BY, BX, C)
  c(d1 / (sigma_X^2),
    d2 / (sigma_X^2),
    -1 * mu_X / (sigma_X^2) * (2 * d1 + d2) + d3,
    -1 * d2 * mu_X / (sigma_X^2) + d4,
    -1 * var_mu_X - var_mu_Y + mu_X^2 / (sigma_X^2) * (d1 + d2) - mu_X * (d3 + d4)
  )
  
  # ## special code for 'complete' correction-term including covariances with Z
  # cov_mu_Z_X    <- std_getCovariance(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
  # cov_sigma_Z_X <- std_getCovariance(mat_A, weighted.sd, group_A, weight_A, JKweights_A)
  # 
  # cov_mu_Z_Y    <- std_getCovariance(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
  # cov_sigma_Z_Y <- std_getCovariance(mat_R, weighted.sd, group_R, weight_R, JKweights_R)
  # 
  # return coefficients, c(AY, AXY, BY, BX, C)
  # t(sapply(1:(length(cov_mu_Z_X)), function(i) {
  #   c(
  #     d1 / (sigma_X^2),
  #     d2 / (sigma_X^2),
  #     -1 * mu_X / (sigma_X^2) * (2 * d1 + d2) + d3 - 2 * cov_sigma_Z_Y[i],
  #     -1 * d2 * mu_X / (sigma_X^2) + d4 - 2 * cov_sigma_Z_X[i],
  #     (-1 * var_mu_X - var_mu_Y + mu_X^2 / (sigma_X^2) * (d1 + d2) - mu_X * (d3 + d4)
  #      + 2 * mu_X * cov_sigma_Z_Y[i] + 2 * mu_X * cov_sigma_Z_X[i]
  #      - 2 * sigma_X * cov_mu_Z_X[i] - 2 * sigma_X * cov_mu_Z_Y[i])
  #   )
  # }))
}



# takes transformed mat_R. takes transformed R
# special code for getting the 'complete' CT
getCTcoefficients_complete <- function(mat_R, mat_A, weight_R, weight_A, JKweights_R, JKweights_A)
{
  # X is A. \theta is untransformed R. Y is transformed R
  M <- ncol(mat_R)
  pop_R <- rep(1,length(weight_R))
  pop_A <- rep(1,length(weight_A))
  
  # get point estimates and BI
  PV_mu_Y <- analyze_PVs(mat_R, weighted.mean, pop_R, weight_R)
  PV_mu_X <- analyze_PVs(mat_A, weighted.mean, pop_A, weight_A)
  PV_sigma_Y <- analyze_PVs(mat_R, weighted.sd, pop_R, weight_R)
  PV_sigma_X <- analyze_PVs(mat_A, weighted.sd, pop_A, weight_A)
  
  # point estimates
  mu_X <- PV_mu_X[1]
  sigma_X <- PV_sigma_X[1]
  
  # total variances
  var_mu_X <- (getWImat(mat_A, weighted.mean, pop_A, weight_A, JKweights_A)
               + (1+1/M)*PV_mu_X[2])
  var_mu_Y <- (getWImat(mat_R, weighted.mean, pop_R, weight_R, JKweights_R)
               + (1+1/M)*PV_mu_Y[2])
  var_sigma_X <- (getWImat(mat_A, weighted.sd, pop_A, weight_A, JKweights_A)
                  + (1+1/M)*PV_sigma_X[2])
  var_sigma_Y <- (getWImat(mat_R, weighted.sd, pop_R, weight_R, JKweights_R)
                  + (1+1/M)*PV_sigma_Y[2])
  
  d1 <- var_sigma_X - var_sigma_Y
  d2 <- -2 * var_sigma_X
  

  # code for longer CT
  var_mu_sigma_X <- (getWImat(mat_A, weighted.msd, pop_A, weight_A, JKweights_A)
                     + (1+1/M)*analyze_PVs(mat_A, weighted.msd, pop_A, weight_A)[2])
  var_mu_sigma_Y <- (getWImat(mat_R, weighted.msd, pop_R, weight_R, JKweights_R)
                     + (1+1/M)*analyze_PVs(mat_R, weighted.msd, pop_R, weight_R)[2])
  cov_mu_sigma_X <- .5 * (var_mu_sigma_X - var_mu_X - var_sigma_X)
  cov_mu_sigma_Y <- .5 * (var_mu_sigma_Y - var_mu_Y - var_sigma_Y)
  
  d3 <- 2 / sigma_X * (-1 * cov_mu_sigma_Y)
  d4 <- 2 / sigma_X * (-1 * cov_mu_sigma_X)
  

  ## special code for 'complete' correction-term including covariances with Z
  cov_mu_Z_X    <- std_getCovariance(mat_A, weighted.mean, group_A, weight_A, JKweights_A)
  cov_sigma_Z_X <- std_getCovariance(mat_A, weighted.sd, group_A, weight_A, JKweights_A)

  cov_mu_Z_Y    <- std_getCovariance(mat_R, weighted.mean, group_R, weight_R, JKweights_R)
  cov_sigma_Z_Y <- std_getCovariance(mat_R, weighted.sd, group_R, weight_R, JKweights_R)

  #return coefficients, c(AY, AXY, BY, BX, C)
  t(sapply(1:(length(cov_mu_Z_X)), function(i) {
    c(
      d1 / (sigma_X^2),
      d2 / (sigma_X^2),
      -1 * mu_X / (sigma_X^2) * (2 * d1 + d2) + d3 - 2 * cov_sigma_Z_Y[i],
      -1 * d2 * mu_X / (sigma_X^2) + d4 - 2 * cov_sigma_Z_X[i],
      (-1 * var_mu_X - var_mu_Y + mu_X^2 / (sigma_X^2) * (d1 + d2) - mu_X * (d3 + d4)
       + 2 * mu_X * cov_sigma_Z_Y[i] + 2 * mu_X * cov_sigma_Z_X[i]
       - 2 * sigma_X * cov_mu_Z_X[i] - 2 * sigma_X * cov_mu_Z_Y[i])
    )
  }))
}

getdf <- function(BI, T, M, max_zones)
{
  FMI <- (1+1/M)*BI / T
  df_c <- max_zones
  
  # Bernard and Rubin
  (df_c + 1)/(df_c + 3) * df_c * (1-FMI)
  
  # PIRLS method
  # 1/(FMI^2/(M-1) + (1-FMI)^2 / d)
}

# getcombineddf <- function(df,T,k)
# {
#   # Welch–Satterthwaite equation
#   k <- as.numeric(k)
#   (rowSums(T * rep(k, each = nrow(T)))^2 
#   / rowSums((T * rep(k, each = nrow(T)))^2 / df))
# }


# Welch–Satterthwaite equation
getcombineddf <- function(df, T, k) { 
  k  <- as.numeric(k)
  
  # If df and T are vectors (not matrices), handle single-group case
  if (is.null(dim(T)) || is.null(dim(df))) {
    numerator <- (sum(T * k))^2
    denominator <- sum((T * k)^2 / df)
    return(numerator / denominator)
  }
  
  # If df and T are matrices, handle row-wise case (return vector of results)
  k_matrix <- matrix(rep(k, each = nrow(T)), nrow = nrow(T))
  numerator <- (rowSums(T * k_matrix))^2
  denominator <- rowSums((T * k_matrix)^2 / df)
  return(numerator / denominator)
}







# List only .Rdata files that start with ASA
files <- list.files(pattern = "^ASA.*\\.Rdata$")

# Filter to ensure filename is at least 8 characters (so 7th letter exists)
files <- files[nchar(files) >= 8]

# Extract the base key (everything except the 7th letter)
base_keys <- unique(substring(files, 1, 6))

# Initialize lists for R and A data
R_data_list <- list()
A_data_list <- list()

big_res_D <- data.frame()

for (base in base_keys) {
  file_R <- files[files == paste0(base, "R", substring(files[1], 8))]
  file_A <- files[files == paste0(base, "A", substring(files[1], 8))]
  
  if (length(file_R) == 1 && length(file_A) == 1) {
    # Load R file
    env_R <- new.env()
    load(file_R, envir = env_R)
    obj_R <- mget(ls(envir = env_R), envir = env_R)[[1]]
    
    # add contextual variables
    file_G <- paste0(substr(file_R, 1, 2), "G", substr(file_R, 4, nchar(file_R)))
    env_G <- new.env()
    load(file_G, envir = env_G)
    obj_G <- mget(ls(envir = env_G), envir = env_G)[[1]]
    new_cols <- setdiff(names(obj_G), names(obj_R))
    obj_R <- cbind(obj_R, obj_G[new_cols])
    
    # Load A file -- paper bridge
    env_A <- new.env()
    load(file_A, envir = env_A)
    obj_A <- mget(ls(envir = env_A), envir = env_A)[[1]]

    # add contextual variables
    file_G <- paste0(substr(file_A, 1, 2), "G", substr(file_A, 4, nchar(file_A)))
    env_G <- new.env()
    load(file_G, envir = env_G)
    obj_G <- mget(ls(envir = env_G), envir = env_G)[[1]]
    new_cols <- setdiff(names(obj_G), names(obj_A))
    obj_A <- cbind(obj_A, obj_G[new_cols])
    
    # get country
    countrycode <- countrycode(obj_R$IDCNTRY[1],"iso3n","iso3c")
    countryname <- countrycode(obj_R$IDCNTRY[1],"iso3n","country.name")
    
    # set weights
    weight_R <- obj_R$TOTWGT
    weight_A <- obj_A$TOTWGT
    
    # construct jackknife weights
    JKweights_R <- with(obj_R, constructJK(JKREP, JKZONE, TOTWGT))
    JKweights_A <- with(obj_A, constructJK(JKREP, JKZONE, TOTWGT))
    
    # extract PVs
    PVmat_R <- as.matrix(obj_R[, sprintf("ASRREA%02d", 1:M)],ncol=M)
    PVmat_A <- as.matrix(obj_A[, sprintf("ASRREA%02d", 1:M)],ncol=M)
    
    
    # untransform R PVs
    PVmat_R <- untrans_PVmat_R <- sweep(sweep(PVmat_R, 2, B, "-"), 2, A, "/")
    
    # retransform R PVs on country-level
    PVmat_R <- linkmat(PVmat_R,PVmat_A,0,0,obj_R$TOTWGT,obj_A$TOTWGT)
    
    
    # choose group
    group_R <- obj_R$ASBG05F
    group_A <- obj_A$ASBG05F
    
    # R statistics
    PVres_R <- with(obj_R,analyze_PVs(PVmat_R, weighted.mean, group_R, TOTWGT))
    WI_R    <- with(obj_R,getWImat(PVmat_R, weighted.mean, group_R, TOTWGT, JKweights_R))
    BI_R    <- PVres_R[,2]
    T_R     <- WI_R + (1+1/M)*BI_R
    df_R    <- getdf(BI_R, T_R, M, max(obj_R$JKZONE))
    
    res_R <- data.frame(countrycode = countrycode,
                        countryname = countryname,
                        means = PVres_R[,1],
                        WI = WI_R,
                        BI = BI_R,
                        T = T_R,
                        df = df_R)

    # A statistics
    PVres_A <- with(obj_A,analyze_PVs(PVmat_A, weighted.mean, group_A, TOTWGT))
    WI_A    <- with(obj_A,getWImat(PVmat_A, weighted.mean, group_A, TOTWGT, JKweights_A))
    BI_A    <- PVres_A[,2]
    T_A     <- WI_A + (1+1/M)*BI_A
    df_A <- getdf(BI_A, T_A, M, max(obj_A$JKZONE))
    
    res_A <- data.frame(countrycode = countrycode,
                        countryname = countryname,
                        means = PVres_A[,1],
                        WI = WI_A,
                        BI = BI_A,
                        T = T_A,
                        df = df_A)

    # conventional stats
    mean_diff <- res_A$means - res_R$means
    T_diff    <- res_R$T + res_A$T
    
    # total resampling
    TR_WI <- gettotalWImat(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A, 
                               obj_R$TOTWGT, obj_A$TOTWGT, JKweights_R, JKweights_A)
    
    TR_BI <- gettotalBI(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A,
                            obj_R$TOTWGT, obj_A$TOTWGT)
    
    T_diff_TR <- TR_WI + (1+1/M)*TR_BI
    
    # component resampling
    CR_WI <- getcomponentWImat(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A, 
                         obj_R$TOTWGT, obj_A$TOTWGT, JKweights_R, JKweights_A)
    
    CR_BI <- getcomponentBI(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A,
                         obj_R$TOTWGT, obj_A$TOTWGT)
    
    T_diff_CR <- CR_WI + (1+1/M)*CR_BI

    # correction-term (simple quad)
    CTcoeff <- getCTcoefficients(PVmat_R, PVmat_A, obj_R$TOTWGT, obj_A$TOTWGT,
                                        JKweights_R, JKweights_A)
    
    # R is Y, A is X
    CT <- (CTcoeff[1] * res_R$means^2 + CTcoeff[2] * res_R$means * res_A$means +
             CTcoeff[3] * res_R$means + CTcoeff[4] * res_A$means + CTcoeff[5])
    
    T_diff_CT <- res_R$T + res_A$T + CT
    T_diff_CTs <- T_diff_CT
    T_diff_CTs[T_diff_CT < 0] <- .Machine$double.eps
    
    # "complete" correction-term including covar with Z
    CTcoeff <- getCTcoefficients_complete(PVmat_R, PVmat_A, obj_R$TOTWGT, obj_A$TOTWGT,
                                          JKweights_R, JKweights_A)
    
    # for the complete correction term method
    cCT <- (CTcoeff[,1] * res_R$means^2 + CTcoeff[,2] * res_R$means * res_A$means +
              CTcoeff[,3] * res_R$means + CTcoeff[,4] * res_A$means + CTcoeff[,5])
    
    T_diff_cCT <- res_R$T + res_A$T + cCT
    

    # equation methods
    T_diff_TEQ <- gettotal_eq(PVmat_R, PVmat_A, group_R, group_A, weight_R, weight_A, JKweights_R, JKweights_A)
    T_diff_CEQ <- getcomponent_eq(PVmat_R, PVmat_A, group_R, group_A, weight_R, weight_A, JKweights_R, JKweights_A)
    
    size_R <- tapply(obj_R$TOTWGT, group_R, sum)/sum(obj_R$TOTWGT)*length(obj_R$TOTWGT)
    size_A <- tapply(obj_A$TOTWGT, group_A, sum)/sum(obj_A$TOTWGT)*length(obj_A$TOTWGT)
    
    
    # linking component
    LC_WI <- getlinkWImat(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A, 
                               obj_R$TOTWGT, obj_A$TOTWGT, JKweights_R, JKweights_A)
    
    LC_BI <- getlinkBI(untrans_PVmat_R, PVmat_A, weighted.mean, group_R, group_A,
                            obj_R$TOTWGT, obj_A$TOTWGT)
    
    T_diff_LC <- T_diff + LC_WI + (1+1/M)*LC_BI
    
    # df and t calculation
    df_diff <- getcombineddf(matrix(c(df_A,df_R),ncol=2),matrix(c(T_A,T_R),ncol=2),c(1,1))
    t <- qt(0.975, df = df_diff)
    
    # D = A - R; paper - digital
    res_D <- data.frame(countrycode = countrycode,
                        countryname = countryname,
                        category_R = as.numeric(names(tapply(group_R, group_R, function(x) x[1]))),
                        category_A = as.numeric(names(tapply(group_A, group_A, function(x) x[1]))),
                        mean_R = res_R$means,
                        mean_A = res_A$means,
                        size_R = size_R,
                        size_A = size_A,
                        pct_size_R = size_R/length(obj_R$TOTWGT),
                        pct_size_A = size_A/length(obj_A$TOTWGT),
                        mean_diff = mean_diff,
                        T_diff = T_diff,
                        lowerbound = mean_diff - t * sqrt(T_diff),
                        upperbound = mean_diff + t * sqrt(T_diff),
                        T_diff_TR = T_diff_TR,
                        TR_lowerbound = mean_diff - t * sqrt(T_diff_TR),
                        TR_upperbound = mean_diff + t * sqrt(T_diff_TR),
                        T_diff_CR = T_diff_CR,
                        CR_lowerbound = mean_diff - t * sqrt(T_diff_CR),
                        CR_upperbound = mean_diff + t * sqrt(T_diff_CR),
                        T_diff_CT = T_diff_CT,
                        CT_lowerbound = mean_diff - t * sqrt(T_diff_CTs),
                        CT_upperbound = mean_diff + t * sqrt(T_diff_CTs),
                        T_diff_TEQ = T_diff_TEQ,
                        TEQ_lowerbound = mean_diff - t * sqrt(T_diff_TEQ),
                        TEQ_upperbound = mean_diff + t * sqrt(T_diff_TEQ),
                        T_diff_CEQ = T_diff_CEQ,
                        CEQ_lowerbound = mean_diff - t * sqrt(T_diff_CEQ),
                        CEQ_upperbound = mean_diff + t * sqrt(T_diff_CEQ),
                        T_diff_LC = T_diff_LC,
                        LC_lowerbound = mean_diff - t * sqrt(T_diff_LC),
                        LC_upperbound = mean_diff + t * sqrt(T_diff_LC),
                        T_diff_cCT = T_diff_cCT,
                        cCT_lowerbound = mean_diff - t * sqrt(T_diff_cCT),
                        cCT_upperbound = mean_diff + t * sqrt(T_diff_cCT),
                        T_diff_CTs = T_diff_CTs,
                        df_diff = df_diff,
                        pval = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff)), df = df_diff)),
                        pval_TR = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_TR)), df = df_diff)),
                        pval_CR = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_CR)), df = df_diff)),
                        pval_CT = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_CTs)), df = df_diff)),
                        pval_TEQ = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_TEQ)), df = df_diff)),
                        pval_CEQ = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_CEQ)), df = df_diff)),
                        pval_LC = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_LC)), df = df_diff)),
                        pval_cCT = 2 * (1 - pt(abs(mean_diff / sqrt(T_diff_cCT)), df = df_diff))
    )
    
    big_res_D <- rbind(big_res_D,res_D)
    
    print(countryname)
    write.csv(big_res_D,"eight_methods_XmY.csv")
  }
}


