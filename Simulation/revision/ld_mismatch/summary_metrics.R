# Figure-3-style evaluation averages; shared by both revision summary scripts.
# `boot` must contain common locus-cluster bootstrap indices (100 rows).
add_evaluation_metrics <- function(d) {
  d[, `:=`(sensitivity_eval=NA_real_, fdp_conditional=NA_real_,
           fdp=fp/pmax(n_flag,1L))]
  d[n_target>0, sensitivity_eval:=tp/n_target]
  d[n_flag>0, fdp_conditional:=fp/n_flag]
  stopifnot(all(is.na(d$sensitivity_eval)==(d$n_target==0)),
            all(is.na(d$fdp_conditional)==(d$n_flag==0)),
            all(d$fdp>=0 & d$fdp<=1))
  invisible(d)
}

cluster_summary <- function(d) {
  a <- d[,.(tp=sum(tp),fp=sum(fp),n_target=sum(n_target),n_non_target=sum(n_non_target),
            n_flag=sum(n_flag),n_snp=sum(n_snp),n_any=sum(n_flag>0),
            sum_sensitivity=sum(sensitivity_eval,na.rm=TRUE),
            n_sensitivity=sum(!is.na(sensitivity_eval)),
            sum_fdp=sum(fdp),n_eval=.N),by=locus]
  setorder(a,locus);stopifnot(identical(a$locus,1:100))
  ratio <- function(num,den,scale=1) {
    x<-a[[num]];y<-a[[den]]
    if(sum(y)==0) return(c(estimate=NA_real_,lower=NA_real_,upper=NA_real_,n_boot_valid=0))
    est<-scale*sum(x)/sum(y)
    bs<-scale*colSums(matrix(x[boot],100))/colSums(matrix(y[boot],100))
    bs<-bs[is.finite(bs)]
    ci<-if(length(bs))quantile(bs,c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)
    c(estimate=est,lower=ci[1],upper=ci[2],n_boot_valid=length(bs))
  }
  specs<-list(sensitivity=c("sum_sensitivity","n_sensitivity"),
              conditional_mean_fdp=c("sum_fdp","n_any"),
              mean_fdp=c("sum_fdp","n_eval"),
              pooled_sensitivity=c("tp","n_target"),
              fp_per_100k=c("fp","n_non_target"),
              flag_per_100k=c("n_flag","n_snp"),any_flag=c("n_any","n_eval"),
              false_discovery_fraction=c("fp","n_flag"))
  rbindlist(lapply(names(specs),function(m) {
    v<-ratio(specs[[m]][1],specs[[m]][2],if(grepl("100k",m))1e5 else 1)
    data.table(metric=m,estimate=v[1],lower=v[2],upper=v[3],n_boot_valid=as.integer(v[4]),
      n_loci=nrow(a),n_evaluations=sum(a$n_eval),
      n_evaluations_with_targets=sum(a$n_sensitivity),
      n_evaluations_with_flags=sum(a$n_any),
      n_target=sum(a$n_target),tp=sum(a$tp),fp=sum(a$fp))
  }))
}

validate_summaries <- function(d,summaries,groups) {
  safe_mean<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)
  direct<-d[,.(sensitivity=safe_mean(sensitivity_eval),
               conditional_mean_fdp=safe_mean(fdp_conditional),mean_fdp=mean(fdp),
               n_flagged=sum(n_flag>0),n_evaluations=.N),by=groups]
  for(m in c("sensitivity","conditional_mean_fdp","mean_fdp")) {
    check<-merge(summaries[metric==m],direct[,c(groups,m,"n_flagged"),with=FALSE],by=groups)
    expected<-check[[m]]
    stopifnot(nrow(check)==nrow(direct),identical(is.na(check$estimate),is.na(expected)),
              all(abs(check$estimate-expected)<1e-12,na.rm=TRUE),
              all(check$n_evaluations_with_flags==check$n_flagged),
              all(check$lower>=0 & check$upper<=1 & check$lower<=check$upper,na.rm=TRUE))
  }
  # E[FDP with zero for no flags] = Pr(any flags) * E[FDP | any flags].
  active<-direct[n_flagged>0]
  stopifnot(all(abs(active$mean_fdp-active$conditional_mean_fdp*
                     active$n_flagged/active$n_evaluations)<1e-12))
  invisible(TRUE)
}

# Edge-case checks independent of the actual simulations.
toy<-data.table(tp=c(9L,1L,0L),fp=c(1L,0L,0L),n_target=c(10L,100L,10L),n_flag=c(10L,1L,0L))
add_evaluation_metrics(toy)
stopifnot(abs(mean(toy$sensitivity_eval)-(.9+.01)/3)<1e-12,
          abs(mean(toy$fdp_conditional,na.rm=TRUE)-.05)<1e-12,
          abs(mean(toy$fdp)-.1/3)<1e-12,
          is.na(toy$fdp_conditional[3]))
rm(toy)
