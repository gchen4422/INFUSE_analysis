# Algebraically equivalent spectral implementation of the existing diagnostics.
# No dense block precision matrix is formed. Each ancestry eigensystem is reused.
# Finite-N z adjustment is performed ONCE, as in the prior revision analysis.
adjust_z <- function(z,n) {
  stopifnot(all(is.finite(z)),all(n>1))
  z*sqrt((n-1)/(z*z+n-2))
}
spectral_input <- function(e,z,lambda=0,permutation=seq_along(z)) {
  stopifnot(lambda>=0,lambda<1,length(e$values)==length(z))
  # Floor AFTER transforming eigenvalues of the supplied LD, not before.
  d <- (1-lambda)*e$values+lambda
  d[d<1e-8] <- 0
  inv <- order(permutation)
  y <- drop(crossprod(e$vectors,z[inv]))
  list(e=e,z=z,d=d,y=y,permutation=permutation)
}
fit_s <- function(inputs) {
  d <- unlist(lapply(inputs,`[[`,"d"),use.names=FALSE)
  y <- unlist(lapply(inputs,`[[`,"y"),use.names=FALSE)
  neglog <- function(s) {
    v <- (1-s)*d+s
    .5*sum(log(v))+.5*sum(y*y/v)
  }
  optim(.5,neglog,method="Brent",lower=0,upper=1)$par
}
conditional_spectral <- function(x,s) {
  dinv <- 1/((1-s)*x$d+s)
  dinv[is.infinite(dinv)] <- 0
  # Diagonal of V diag(dinv) V' and its product with z.
  pd <- drop(x$e$squared %*% dinv)[x$permutation]
  pz <- drop(x$e$vectors %*% (dinv*x$y))[x$permutation]
  stopifnot(all(is.finite(pd)),all(pd>0))
  residual <- pz/pd
  data.frame(z=x$z,condmean=x$z-residual,condvar=1/pd,
             z_std_diff=pz/sqrt(pd))
}
mixture_lr <- function(cd) {
  amax <- if(max(cd$z_std_diff^2)<1) 2 else 2*sqrt(max(cd$z_std_diff^2))
  npoint <- ceiling(log2(amax/.8)/log2(1.05))
  scales <- 1.05^(seq(-npoint,0))*amax
  sds <- outer(sqrt(cd$condvar),scales)
  ll <- dnorm(cd$z-cd$condmean,sd=sds,log=TRUE)
  offsets <- apply(ll,1,max)
  ll <- ll-offsets
  w <- mixsqp::mixsqp(ll,log=TRUE,control=list(verbose=FALSE))$x
  l0 <- drop(log(exp(ll)%*%(w+1e-15)))+offsets
  ll <- dnorm(cd$z+cd$condmean,sd=sds,log=TRUE)
  offsets <- apply(ll,1,max)
  ll <- ll-offsets
  l1 <- drop(log(exp(ll)%*%(w+1e-15)))+offsets
  cd$logLR <- l1-l0
  stopifnot(all(is.finite(cd$logLR)))
  cd
}
joint_from_conditionals <- function(cd1,cd2) {
  p <- nrow(cd1)
  cd <- mixture_lr(rbind(cd1,cd2))
  lr1 <- cd$logLR[seq_len(p)];lr2 <- cd$logLR[p+seq_len(p)]
  r1 <- abs(cd$z_std_diff[seq_len(p)])>4
  r2 <- abs(cd$z_std_diff[p+seq_len(p)])>4
  f1 <- which(r1 & (lr1>2 | lr1+lr2>2))
  f2 <- which(r2 & (lr2>2 | lr1+lr2>2))
  list(flags=list(WB=f1,BB=f2,Either=sort(union(f1,f2))),scores=cd)
}
run_diagnostics <- function(x1,x2,include_union=TRUE,include_separate=FALSE) {
  sh <- fit_s(list(x1,x2))
  ans <- list(joint_shared=joint_from_conditionals(
    conditional_spectral(x1,sh),conditional_spectral(x2,sh)))
  sep <- c(NA_real_,NA_real_)
  if(include_union || include_separate) {
    sep <- c(fit_s(list(x1)),fit_s(list(x2)))
    c1 <- conditional_spectral(x1,sep[1]);c2 <- conditional_spectral(x2,sep[2])
    if(include_union) {
      a <- mixture_lr(c1);b <- mixture_lr(c2)
      f1 <- which(a$logLR>2 & abs(a$z)>2)
      f2 <- which(b$logLR>2 & abs(b$z)>2)
      ans$separate_union <- list(flags=list(WB=f1,BB=f2,Either=sort(union(f1,f2))),
                                 scores=rbind(a,b))
    }
    if(include_separate) ans$joint_separate <- joint_from_conditionals(c1,c2)
  }
  list(methods=ans,s=c(shared=sh,WB=sep[1],BB=sep[2]))
}
count_metrics <- function(flags,truth,p) {
  stopifnot(all(flags %in% seq_len(p)),all(truth %in% seq_len(p)))
  tp <- length(intersect(flags,truth));fp <- length(setdiff(flags,truth))
  data.frame(n_snp=p,n_target=length(truth),n_non_target=p-length(truth),
             n_flag=length(flags),tp=tp,fp=fp,fn=length(truth)-tp)
}
