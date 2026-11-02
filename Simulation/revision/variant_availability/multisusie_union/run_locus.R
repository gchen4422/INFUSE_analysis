#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Fit the union benchmark using saved simulations and causal configurations.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 2L))
locus <- as.integer(args[[1L]])
pilot <- length(args) == 2L && args[[2L]] == "pilot"
stopifnot(locus %in% 1:100)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L]))
OUT <- file.path(.infuse_work, "INFUSE_variant_availability/simulation_native_assv/multisusie_union_test_20261020")
NATIVE <- dirname(OUT)
ROOT <- dirname(NATIVE)
suppressPackageStartupMessages({library(data.table); library(snpStats); library(reticulate)})
ANCESTRY <- c("EUR", "AFR")
BASE <- c(2L, 3L, 5L)
cell_path <- function(b, h, a, type="cells") file.path(NATIVE, "results_native_sim", type,
  paste0("base_", BASE[b]), paste0("h2_", h), paste0("assv_each_", a), sprintf("locus_%03d.rds", locus))
output_path <- function(b,h,a) file.path(OUT,"results",paste0("base_",BASE[b]),
  paste0("h2_",h),paste0("assv_each_",a),sprintf("locus_%03d.rds",locus))
atomic_rds <- function(x,p) {
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  q <- paste0(p,".partial_",Sys.getpid())
  saveRDS(x,q,compress=FALSE)
  if(!file.rename(q,p)) stop("Could not finalize ",p)
}
grid <- expand.grid(base_index=if(pilot)2L else 1:3, h2_index=if(pilot)1L else 1:2,
                    assv_each=0:2)
if(all(vapply(seq_len(nrow(grid)),function(i)file.exists(output_path(grid[i,1],grid[i,2],grid[i,3])),logical(1)))) {
  cat("Already complete; no files overwritten\n"); quit(status=0)
}
panel <- file.path(NATIVE,"panels_source_preserved",sprintf("locus_%03d",locus))
stopifnot(file.exists(file.path(panel,"COMPLETE")))
annot <- readRDS(file.path(panel,"annotation.rds"))
ids <- rownames(annot); rm(annot)
vm <- fread(file.path(panel,"variant_map.tsv"))
hit <- match(ids,vm$SNP); stopifnot(!anyNA(hit),!anyDuplicated(hit),!anyDuplicated(ids))
vm <- vm[hit]
available <- as.matrix(vm[,.(available_EUR,available_AFR)])
storage.mode(available) <- "logical"
dimnames(available) <- list(ids,ANCESTRY)
stopifnot(all(rowSums(available)>0L))
p <- length(ids)
common <- which(rowSums(available)==2L)
observed_idx <- setNames(lapply(ANCESTRY,function(a)which(available[,a])),ANCESTRY)
input_dir <- file.path(OUT,"inputs"); dir.create(input_dir,showWarnings=FALSE)
maf_path <- file.path(input_dir,sprintf("locus_%03d_maf.rds",locus))
if(file.exists(maf_path)) {
  cache <- readRDS(maf_path)
  stopifnot(identical(cache$snp,ids),identical(cache$available,available))
  maf <- cache$maf
} else {
  maf <- setNames(vector("list",2L),ANCESTRY)
  for(a in ANCESTRY) {
    pl <- read.plink(file.path(panel,"plink",paste0(a,".bed")))
    ii <- observed_idx[[a]]
    jj <- match(ids[ii],as.character(pl$map$snp.name))
    stopifnot(!anyNA(jj),!anyDuplicated(jj))
    observed_maf <- as.numeric(col.summary(pl$genotypes[,jj,drop=FALSE])$MAF)
    stopifnot(all(is.finite(observed_maf)),all(observed_maf>0 & observed_maf<=.5))
    # Zero outside availability is a censoring marker, not a biological MAF claim.
    maf[[a]] <- setNames(numeric(p),ids)
    maf[[a]][ii] <- observed_maf
    rm(pl); invisible(gc())
  }
  old_maf <- readRDS(file.path(NATIVE,"multisusie_inputs",sprintf("locus_%03d_shared_maf.rds",locus)))
  stopifnot(identical(old_maf$snp,ids[common]))
  for(a in ANCESTRY)stopifnot(max(abs(maf[[a]][common]-old_maf$maf[[a]]))<1e-12)
  atomic_rds(list(snp=ids,available=available,maf=maf,
                 missing_marker="MAF=0 only to censor unavailable entries; not inferred monomorphism"),maf_path)
}
R_union <- setNames(vector("list",2L),ANCESTRY)
for(a in ANCESTRY) {
  rr <- readRDS(file.path(panel,"ld",paste0(a,".rds")))
  ii <- observed_idx[[a]]; jj <- match(ids[ii],rownames(rr))
  stopifnot(!anyNA(jj),identical(rownames(rr),colnames(rr)))
  r <- matrix(0,p,p)
  r[ii,ii] <- rr[jj,jj,drop=FALSE]
  stopifnot(all(is.finite(r)),max(abs(diag(r)[ii]-1))<1e-6,
            all(r[!available[,a],,drop=FALSE]==0),all(r[,!available[,a],drop=FALSE]==0))
  R_union[[a]] <- r
  rm(rr,r);invisible(gc())
}
use_python(Sys.getenv("INFUSE_PYTHON_BIN", Sys.which("python")),required=TRUE)
ms <- import("MultiSuSiE",convert=FALSE)
np <- import("numpy",convert=FALSE)
random <- import("random",convert=FALSE)
hashlib <- import("hashlib",convert=FALSE)
builtins <- import_builtins(convert=FALSE)
pkg_path <- normalizePath(py_to_r(ms$`__file__`))
source_path <- file.path(dirname(pkg_path),"susiepy_ss.py")
source_sha <- py_to_r(hashlib$sha256(builtins$open(source_path,"rb")$read())$hexdigest())
stopifnot(identical(names(R_union),ANCESTRY),identical(names(maf),ANCESTRY))
R_py <- unname(lapply(R_union,function(r)np$array(r,dtype=np$float32)))
maf_py <- unname(lapply(maf,function(x)np$array(unname(x))))
rho_py <- np$array(matrix(c(1,.8,.8,1),2L,byrow=TRUE))
fit_ms <- function(z, r=R_py, ma=maf_py, n, seed) {
  random$seed(as.integer(seed));np$random$seed(as.integer(seed))
  ms$multisusie_rss(z_list=unname(lapply(z,np$array)),R_list=unname(r),rho=rho_py,
    population_sizes=np$array(as.integer(n)),L=as.integer(10L),
    scaled_prior_variance=.2,low_memory_mode=FALSE,min_abs_corr=.5,
    single_population_mac_thresh=as.integer(20L),maf_list=unname(ma),coverage=.95)
}
metric_tables <- function(fits,truth,b,h,a) {
  truth_flag <- truth$Signal!=0L
  origin <- as.character(truth$causal_origin)
  cs_rows <- roc_rows <- cal_rows <- causal_rows <- list()
  for(m in names(fits)) {
    f <- fits[[m]];pip <- as.numeric(f$pip)
    # Some saved MESuSiE PIPs equal 1 + 2.22e-16. Preserve these exact
    # values and threshold semantics; tolerate only machine-scale excursions.
    stopifnot(length(pip)==p,all(is.finite(pip)),all(pip>=-1e-12 & pip<=1+1e-12))
    selected <- sort(unique(as.integer(unlist(f$cs,use.names=FALSE))))
    stopifnot(all(selected>=1L & selected<=p))
    tp_cs <- sum(truth_flag[selected]);ncs <- length(f$cs)
    cs_rows[[m]] <- data.table(Method=m,Power=tp_cs/sum(truth_flag),
      FDR=if(length(selected))1-tp_cs/length(selected) else 0,Size=length(selected),
      n_cs=ncs,median_component_size=if(ncs)median(lengths(f$cs)) else NA_real_)
    roc_rows[[m]] <- rbindlist(lapply(seq(1,.01,-.01),function(cut) {
      keep <- pip>cut;tp <- sum(keep & truth_flag);fp <- sum(keep & !truth_flag)
      data.table(Method=m,Cutoff=cut,Power=tp/sum(truth_flag),
        FDR=if(tp+fp)fp/(tp+fp) else 0,TP=tp,FP=fp,n_causal=sum(truth_flag))
    }))
    # Matches existing calibration cut(): zero PIPs excluded from (0,0.1].
    cb <- cut(pip,breaks=seq(0,1,.1),labels=FALSE)
    cal_rows[[m]] <- data.table(Method=m,bin=cb,pip=pip,causal=truth_flag)[!is.na(bin),
      .(n=.N,sum_pip=sum(pip),n_causal=sum(causal)),by=.(Method,bin)]
    jj <- which(truth_flag)
    causal_rows[[m]] <- data.table(Method=m,SNP=ids[jj],Signal=truth$Signal[jj],
      causal_origin=origin[jj],PIP=pip[jj],in_cs=jj%in%selected,
      available_EUR=available[jj,"EUR"],available_AFR=available[jj,"AFR"])
  }
  add <- function(x) {d<-rbindlist(x);d[,`:=`(locus=locus,base_index=b,n_base=BASE[b],
    h2_index=h,h2_per_causal=c(1e-4,2e-4)[h],assv_each=a)];d}
  list(cs=add(cs_rows),roc=add(roc_rows),calibration=add(cal_rows),causal=add(causal_rows))
}
for(g in seq_len(nrow(grid))) {
  b <- grid$base_index[g];h <- grid$h2_index[g];a <- grid$assv_each[g]
  out <- output_path(b,h,a)
  if(file.exists(out))next
  input <- cell_path(b,h,a)
  cell <- readRDS(input)
  stopifnot(identical(as.character(cell$truth$SNP),ids),
    identical(unname(cell$availability),unname(available)),identical(cell$metadata$L,10L))
  z <- setNames(lapply(ANCESTRY,function(anc) {
    col <- c(EUR="zscore_1",AFR="zscore_2")[[anc]]
    zz <- as.numeric(cell$zfile[[col]])
    stopifnot(all(is.finite(zz[available[,anc]])))
    zz[!available[,anc]] <- 0
    zz
  }),ANCESTRY)
  n <- vapply(ANCESTRY,function(anc) {
    col <- c(EUR="N_1",AFR="N_2")[[anc]]
    nn <- cell$zfile[[col]][available[,anc]]
    stopifnot(all(is.finite(nn)),length(unique(nn))==1L)
    as.numeric(unique(nn))
  },numeric(1))
  stopifnot(all(n==300000))
  old <- readRDS(cell_path(b,h,a,"multisusie_intersect"))
  seed <- as.integer(20261020+locus*1000+b*100+h*10+a)
  # Reproduce an existing intersection fit before accepting the union pilot.
  if(pilot && a==0L) {
    rp <- lapply(R_union,function(r)np$array(r[common,common,drop=FALSE]))
    mp <- lapply(maf,function(v)np$array(unname(v[common])))
    same <- fit_ms(lapply(z,function(v)v[common]),rp,mp,n,seed)
    delta <- max(abs(as.numeric(py_to_r(same$pip))-old$pip[common]))
    stopifnot(is.finite(delta),delta<1e-6)
    atomic_rds(list(locus=locus,pip_max_absolute_difference=delta,
      package_source_sha256=source_sha,settings="same as existing intersection runner"),
      file.path(OUT,"validation",sprintf("intersection_parity_locus_%03d.rds",locus)))
    rm(same,rp,mp);invisible(gc())
  }
  cat(sprintf("FIT locus=%d base=%d h2=%d assv=%d union=%d common=%d\n",locus,BASE[b],h,a,p,length(common)))
  started <- proc.time()[[3L]]
  f <- fit_ms(z,n=n,seed=seed)
  runtime <- (proc.time()[[3L]]-started)/60
  pip <- as.numeric(py_to_r(f$pip))
  alpha <- py_to_r(f$alpha)
  stopifnot(length(pip)==p,all(is.finite(pip)),all(pip>=0 & pip<=1),
            all(is.finite(alpha)),max(abs(rowSums(alpha)-1))<1e-5)
  sets <- py_to_r(f$sets);inc <- as.logical(sets[[4L]])
  cs <- lapply(which(inc),function(j)as.integer(sets[[1L]][[j]])+1L)
  stopifnot(all(unlist(cs)>=1L & unlist(cs)<=p))
  union_fit <- list(pip=pip,cs=cs)
  fits <- setNames(cell$methods,c("INFUSE union","MESuSiE union","INFUSE intersect",
    "MESuSiE intersect","SuSiE","SuSiEx intersect","SuSiEx union"))
  stopifnot(identical(names(cell$methods),c("infuse_union","mesusie_union","infuse_intersect",
    "mesusie_intersect","susie","susiex_intersect","susiex_union")))
  fits[["INFUSE-H"]] <- readRDS(cell_path(b,h,a,"infuse_union_d187"))
  fits[["MultiSuSiE intersect"]] <- old
  fits[["MultiSuSiE union"]] <- union_fit
  elbo <- as.numeric(py_to_r(f$elbo))
  stopifnot(all(is.finite(elbo)))
  result <- list(metadata=list(locus=locus,base_index=b,n_base=BASE[b],h2_index=h,
      h2_per_causal=c(1e-4,2e-4)[h],assv_each=a,variant_set="union",n_union=p,
      n_common=length(common),n_EUR_only=sum(available[,1]&!available[,2]),
      n_AFR_only=sum(available[,2]&!available[,1]),seed=seed,
      population_sizes=n,L=10L,rho=.8,scaled_prior_variance=.2,coverage=.95,
      min_abs_corr=.5,single_population_mac_thresh=20L,package_path=pkg_path,
      package_source_sha256=source_sha,source_cell=input,
      padding="Z=0 and full LD row/column including diagonal=0 when unavailable; MAF=0 censoring marker"),
    pip=pip,cs=cs,alpha=alpha,coef=py_to_r(f$coef),coef_sd=py_to_r(f$coef_sd),
    converged=isTRUE(py_to_r(f$converged)),niter=as.integer(py_to_r(f$niter)),
    elbo=elbo,runtime_minutes=runtime,metrics=metric_tables(fits,cell$truth,b,h,a))
  atomic_rds(result,out)
  cat(sprintf("DONE locus=%d base=%d h2=%d assv=%d n_cs=%d converged=%s minutes=%.3f\n",
    locus,BASE[b],h,a,length(cs),result$converged,runtime))
  rm(cell,old,f,result,fits,union_fit,alpha);invisible(gc())
}
