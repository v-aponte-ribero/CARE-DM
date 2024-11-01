pool_cindx <- function(cindx,se,n,transform=T){
  
  est.logit <- log(cindx/(1-cindx))
  
  pool <- pool.scalar(Q=est.logit,
                      U=var.logit,
                      n=n)
  
  est <- pool$qbar
  lo <- pool$qbar-qnorm(0.975)*sqrt(pool$t)
  hi <-pool$qbar+qnorm(0.975)*sqrt(pool$t)
  
  if(transform==TRUE){
    est <- exp(pool$qbar)/(1+exp(pool$qbar))
    lo <- exp(lo)/(1+exp(lo))
    hi <- exp(hi)/(1+exp(hi))
  }
  
  
  ret <- c(Estimate=est,
           logitSE=sqrt(pool$t),
           lower=lo,
           upper=hi)
  
  return(ret)
  
}

pool_cal <- function(est,se,n){
  
  var <- se^2
  
  pool <- pool.scalar(Q=est,
                      U=var,
                      n=n)
  
  ret <- c(Estimate=pool$qbar,
           SE=sqrt(pool$t),
           lower=pool$qbar-qnorm(0.975)*sqrt(pool$t),
           upper=pool$qbar+qnorm(0.975)*sqrt(pool$t))
  
  return(ret)
  
}

pool_oe <- function(est,se,n,transform=T){
  
  var <- se^2
  
  pool <- pool.scalar(Q=log(est),
                      U=var,
                      n=n)
  
  est <- pool$qbar
  lo <- pool$qbar-qnorm(0.975)*sqrt(pool$t)
  hi <- pool$qbar+qnorm(0.975)*sqrt(pool$t)
  
  if(transform==TRUE){
    est <- exp(pool$qbar)
    lo <- exp(lo)
    hi <- exp(hi)
  }
  
  ret <- c(Estimate=est,
           logSE=sqrt(pool$t),
           lower=lo,
           upper=hi)
  
  return(ret)
  
}
