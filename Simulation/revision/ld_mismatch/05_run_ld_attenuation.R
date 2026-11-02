#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Joint and separate LD diagnostics under bilateral or one-ancestry attenuation.
stopifnot(nzchar(Sys.getenv("SLURM_JOB_ID")))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
base <- .infuse_project
module <- file.path(.infuse_work, "LD_mismatch_sensitivity")
source(file.path(.infuse_repo, "Simulation/revision/ld_mismatch/ld_test_engine.R"))
locus <- as.integer(commandArgs(TRUE)[1]); stopifnot(locus %in% 1:100)
outdir <- file.path(module,"results/ld_attenuation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
outfile <- file.path(outdir,sprintf("locus_%03d.rds",locus))
if(file.exists(outfile)) stop("Refusing overwrite: ",outfile)
zpath <- function(c,h,rho) file.path(base,"More_annotes",paste0("causal_num_",c),"summary_data",
  paste0("flipped_0.01_rho_",rho),sprintf("CAUSAL_%d_LOCI_%d_h2_%d_flipped_0.01_rho_%d.txt",c,locus,h,rho))
ldpath <- function(c,h,k) file.path(Sys.getenv("INFUSE_SHARED_LD", file.path(dirname(base),"shared_50")),paste0("causal_num_",c),"summary_data",
  sprintf("CAUSAL_%d_LOCI_%d_h2_%d.LD%d",c,locus,h,k))
target_locus <- locus
saved_truth <- fread(file.path(module,"inputs/flip_truth.tsv"))[locus==target_locus]
ref <- fread(zpath(1,1,4)); p <- nrow(ref)
start <- proc.time()[3]; eigens <- list(); input_paths <- character()
for(k in 1:2) {
  path <- ldpath(1,1,k); R <- as.matrix(fread(path))
  stopifnot(nrow(R)==p,ncol(R)==p,max(abs(R-t(R)))<1e-7,max(abs(diag(R)-1))<1e-7)
  refhash <- unname(tools::md5sum(path))
  for(c in 1:3) for(h in 1:2) {
    path_ch <- ldpath(c,h,k)
    if(unname(tools::md5sum(path_ch))!=refhash) stopifnot(max(abs(as.matrix(fread(path_ch))-R))<1e-12)
    input_paths <- c(input_paths,path_ch)
  }
  e <- eigen(R,symmetric=TRUE); e$squared <- e$vectors^2; eigens[[k]] <- e
  rm(R,e); gc()
}
cat("Eigen complete",locus,"SNPs",p,"seconds",proc.time()[3]-start,"\n");flush.console()
metrics <- list(); flags <- list(); parameters <- list()
record <- function(ans,c,h,fam,lambda,truth,z) {
  info <- data.table(locus=locus,causal_index=c,h2_index=h,family=fam,level=lambda)
  if(!is.null(ans$s)) parameters[[length(parameters)+1L]] <<- cbind(info,as.data.table(as.list(ans$s)))
  for(mid in names(ans$methods)) for(sc in c("WB","BB","Either")) {
    fl <- ans$methods[[mid]]$flags[[sc]]
    mm <- as.data.table(count_metrics(fl,truth[[sc]],p))
    ii <- cbind(info,data.table(method=mid,scope=sc))
    metrics[[length(metrics)+1L]] <<- cbind(ii,mm)
    if(length(fl)) flags[[length(flags)+1L]] <<- cbind(ii,data.table(index=fl,RSID=z$RSID[fl],is_target=fl %in% truth[[sc]]))
  }
}
for(c in 1:3) for(h in 1:2) {
  z <- fread(zpath(c,h,4)); rho1 <- fread(zpath(c,h,1))
  input_paths <- c(input_paths,zpath(c,h,4),zpath(c,h,1))
  stopifnot(identical(z$RSID,ref$RSID),identical(z$RSID,rho1$RSID),all(z$N_1==300000),all(z$N_2==300000))
  truth <- list()
  for(k in 1:2) {
    col <- paste0("zscore_",k)
    idx <- which(abs(z[[col]]-sign(rho1[[col]])*abs(z[[col]]))>1e-10)
    saved <- saved_truth[causal_index==c & h2_index==h & ancestry==k]
    stopifnot(identical(idx,sort(saved$row_index)),length(idx)==ceiling(.01*p))
    truth[[c("WB","BB")[k]]] <- idx
  }
  truth$Either <- sort(union(truth$WB,truth$BB))
  zz <- list(adjust_z(z$zscore_1,z$N_1),adjust_z(z$zscore_2,z$N_2))
  for(fam in c("shrink_both","shrink_WB","shrink_BB")) {
    levels <- if(fam=="shrink_both") c(0,.01,.05,.10,.25,.50) else c(0,.01,.05,.10,.25,.50,.75)
    for(lambda in levels) {
      l1 <- if(fam=="shrink_BB") 0 else lambda
      l2 <- if(fam=="shrink_WB") 0 else lambda
      ans <- run_diagnostics(spectral_input(eigens[[1]],zz[[1]],l1),
        spectral_input(eigens[[2]],zz[[2]],l2),include_union=TRUE,include_separate=TRUE)
      record(ans,c,h,fam,lambda,truth,z)
    }
  }
  cat("Completed",locus,c,h,"seconds",proc.time()[3]-start,"\n");flush.console()
}
new <- rbindlist(metrics)
keys <- c("locus","causal_index","h2_index","family","level","method","scope")
stopifnot(!anyDuplicated(new[,..keys]))
complete <- new
stopifnot(nrow(complete)==1080L,!anyDuplicated(complete[,..keys]),
  all(complete$tp+complete$fp==complete$n_flag),all(complete$tp+complete$fn==complete$n_target))
out <- list(metrics=complete,new_flags=rbindlist(flags),parameters=rbindlist(parameters),
  seconds=proc.time()[3]-start,input_md5=tools::md5sum(unique(input_paths)),
  source_md5=tools::md5sum(file.path(.infuse_repo, c("Simulation/revision/ld_mismatch/05_run_ld_attenuation.R", "Simulation/revision/ld_mismatch/ld_test_engine.R"))),
  session=capture.output(sessionInfo()))
saveRDS(out,paste0(outfile,".partial"),compress=FALSE)
stopifnot(file.rename(paste0(outfile,".partial"),outfile))
