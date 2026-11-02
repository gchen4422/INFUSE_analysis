#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Summarize conditional and unconditional FDP from sparse LD perturbations.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
revision <- .infuse_work
input <- file.path(revision, "res/LD_mismatch_sensitivity/no_z_filter_ld_r0p1/LD_mismatch_sparse_no_z_filter_summary.RData")
out <- file.path(revision, "res/LD_mismatch_sensitivity/no_z_filter_ld_r0p1/conditional_fdp")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
e <- new.env()
load(input, envir=e)
tests <- c("separate_union", "joint", "joint_separate_s")
d <- as.data.table(e$metrics)[scenario_class == "sparse" & test %in% tests]
key <- c("scenario", "test", "locus", "causal_num", "h2_num")
stopifnot(nrow(d) == 4L*3L*100L*6L, !anyDuplicated(d[, ..key]),
          all(d$n_flagged == d$TP+d$FP), all(d$n_truth == d$TP+d$FN),
          all(is.finite(d$TP)), all(is.finite(d$FP)), all(d$n_truth > 0))
d[, `:=`(fdp_conditional= fifelse(n_flagged > 0, FP/pmax(1, n_flagged), NA_real_),
          fdp_unconditional=FP/pmax(1, n_flagged), sensitivity_eval=TP/n_truth,
          causal_count=c(2L,3L,5L)[causal_num], h2=c(1e-4,2e-4)[h2_num])]
stopifnot(all(abs(d$FDR - d$fdp_unconditional) < 1e-12))
loci <- sort(unique(d$locus))
stopifnot(length(loci)==100L)
set.seed(20260917L)
B <- 2000L
boot <- replicate(B, tabulate(sample.int(100L, 100L, replace=TRUE), nbins=100L))
summarize <- function(x) {
  a <- x[, .(sum_fdp=sum(fdp_unconditional), n_discovery=sum(n_flagged>0),
             n_eval=.N, sum_sensitivity=sum(sensitivity_eval)), by=locus]
  a <- a[match(loci, locus)]
  stopifnot(!anyNA(a))
  den <- drop(crossprod(a$n_discovery, boot))
  num <- drop(crossprod(a$sum_fdp, boot))
  valid <- den > 0
  ci <- if(any(valid)) quantile(num[valid]/den[valid], c(.025,.975), names=FALSE) else c(NA_real_, NA_real_)
  n_positive <- sum(x$n_flagged>0)
  list(n_evaluations=nrow(x), n_with_discoveries=n_positive,
       n_without_discoveries=sum(x$n_flagged==0), TP=sum(x$TP), FP=sum(x$FP),
       n_flagged=sum(x$n_flagged), conditional_mean_fdp=if(n_positive) mean(x$fdp_conditional, na.rm=TRUE) else NA_real_,
       conditional_fdp_low=ci[1], conditional_fdp_high=ci[2], valid_bootstrap_replicates=sum(valid),
       unconditional_mean_fdp=mean(x$fdp_unconditional),
       pooled_fdp=if(sum(x$n_flagged)>0) sum(x$FP)/sum(x$n_flagged) else NA_real_,
       mean_sensitivity=mean(x$sensitivity_eval), pooled_sensitivity=sum(x$TP)/sum(x$n_truth))
}
by_arch <- d[, summarize(.SD), by=.(scenario, anchor_scale, test, causal_count, h2)]
overall <- d[, summarize(.SD), by=.(scenario, anchor_scale, test)]
# Validate that pooling the raw counts reproduces the existing figure's source.
old <- as.data.table(e$sparse_metrics)
old[, `:=`(causal_count=c(2L,3L,5L)[causal_num], h2=c(1e-4,2e-4)[h2_num])]
check <- merge(by_arch, old[, .(scenario,test,causal_count,h2,old_TP=TP,old_FP=FP,old_FDR=FDR)],
               by=c("scenario","test","causal_count","h2"))
stopifnot(nrow(check)==72L, all(check$TP==check$old_TP), all(check$FP==check$old_FP))
finite <- is.finite(check$pooled_fdp) & is.finite(check$old_FDR)
stopifnot(all(abs(check$pooled_fdp[finite]-check$old_FDR[finite]) < 1e-12))
fwrite(by_arch, file.path(out,"conditional_fdp_by_architecture.tsv"), sep="\t", na="NA")
fwrite(overall, file.path(out,"conditional_fdp_overall.tsv"), sep="\t", na="NA")
fwrite(d[, c(key,"TP","FP","n_flagged","n_truth","fdp_conditional","fdp_unconditional","sensitivity_eval"), with=FALSE],
       file.path(out,"evaluation_metrics.tsv"), sep="\t", na="NA")
saveRDS(list(seed=20260917L, reps=B, loci=loci, weights=boot, source=input,
             source_md5=unname(tools::md5sum(input))), file.path(out,"provenance_bootstrap.rds"))
print(overall)
