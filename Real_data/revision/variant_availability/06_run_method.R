#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Run fine-mapping on QC-filtered observed inputs with matched annotations.
suppressPackageStartupMessages({library(data.table);library(INFUSE);library(susieR)})
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
args <- commandArgs(TRUE); task <- as.integer(args[1]); method <- args[2]
susie_max_iter <- if(length(args)>2L)as.integer(args[3]) else 100L
stopifnot(length(susie_max_iter)==1L,is.finite(susie_max_iter),susie_max_iter>=100L)
stopifnot(method %in% c("INFUSE_union","MESuSiE_union","SuSiE_posthoc"))
man <- fread(file.path(root,"inputs","frozen_manifests","region_manifest.tsv"))[task_id==task]
stopifnot(nrow(man)==1L,isTRUE(man$qc_pass))
input_dir <- file.path(root,"inputs","prepared",sprintf("%s_%03d",man$trait,man$locus))
stopifnot(file.exists(file.path(input_dir,"PASS")))
input_file <- file.path(input_dir,"input.rds");x <- readRDS(input_file)
outdir <- file.path(root,"results",method,sprintf("%s_%03d",man$trait,man$locus))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
attempt <- file.path(outdir,paste0("attempt_",Sys.getenv("SLURM_JOB_ID",unset=Sys.getpid())))
dir.create(attempt,showWarnings=FALSE)
if(file.exists(file.path(outdir,"SUCCESS.rds"))) stop("Successful result exists; refusing rerun")
warnings <- character();set.seed(20260910L+task)
t0 <- proc.time()[["elapsed"]]
tryCatch(withCallingHandlers({
  p <- length(x$union_ids)
  settings <- list(L=10L,coverage=.95,min_abs_corr=.5,seed=20260910L+task)
  if(method=="SuSiE_posthoc"){
    fits <- lapply(c("EUR","AFR"),function(anc){
      ss <- x$summary_stat_list[[anc]]
      susieR::susie_rss(z=ss$Z,R=x$R_mat_list[[anc]],n=median(ss$N),L=10L,
        coverage=.95,min_abs_corr=.5,max_iter=susie_max_iter)
    });names(fits)<-c("EUR","AFR")
    P <- matrix(0,p,2,dimnames=list(x$union_ids,c("EUR","AFR")))
    cs <- list()
    for(anc in names(fits)){
      ss <- x$summary_stat_list[[anc]];ii <- match(ss$SNP,x$union_ids)
      P[ii,anc] <- fits[[anc]]$pip
      for(j in seq_along(fits[[anc]]$sets$cs)) cs[[paste(anc,j,sep="_")]] <- ii[fits[[anc]]$sets$cs[[j]]]
    }
    tab <- data.table(SNP=x$union_ids,pip_either=pmax(P[,1],P[,2]),pip_shared=pmin(P[,1],P[,2]),
      pip_EUR=P[,1],pip_AFR=P[,2])
    tab[,pip_unique:=fifelse((pip_EUR>.5 & pip_AFR<.5)|(pip_AFR>.5 & pip_EUR<.5),pip_either,0)]
    converged <- all(vapply(fits,function(f)isTRUE(f$converged),logical(1)))
    diagnostics <- list(converged=converged,niter=lapply(fits,`[[`,"niter"),elbo=lapply(fits,`[[`,"elbo"))
    settings$n <- vapply(x$summary_stat_list,function(ss)median(ss$N),numeric(1))
    settings$max_iter <- susie_max_iter
    settings$unique_definition <- "fixed strict >0.5/<0.5 ancestry PIPs; Shared=min; Either=max"
    saveRDS(fits,file.path(attempt,"raw_fit.rds"))
  }else{
    prep <- INFUSE::prepare_infuse_union(x$summary_stat_list,x$R_mat_list,
      annot=if(method=="INFUSE_union")x$annot else NULL,union_ids=x$union_ids,
      allele_ref=x$allele_ref,ld_allele_coding="ancestry_effect",verbose=FALSE)
    settings <- c(settings,list(optim_method="optim",est_annot_prior="fixed",
      annot_method=if(method=="INFUSE_union")"mlk" else NULL,n_annotations=if(method=="INFUSE_union")186L else 0L,
      max_iter=100L,estimate_residual_variance=FALSE,ancestry_weight=c(3/7,3/7,1/7)))
    fit <- INFUSE::infuse_core_avail(R_mat_list=prep$R_mat_list,
      summary_stat_list=prep$summary_stat_list,L=10L,variant_available=prep$variant_available,
      annot=prep$annot,annot_method=settings$annot_method,est_annot_prior="fixed",optim_method="optim",
      ancestry_weight=c(3/7,3/7,1/7),estimate_residual_variance=FALSE,max_iter=100L,cor_threshold=.5,verbose=TRUE)
    cfg <- fit$pip_config
    stopifnot(all(c("EUR","AFR","EUR+AFR") %in% colnames(cfg)))
    tab <- data.table(SNP=x$union_ids,pip_either=as.numeric(fit$pip),pip_shared=cfg[,"EUR+AFR"],pip_EUR=cfg[,"EUR"],pip_AFR=cfg[,"AFR"])
    tab[,pip_unique:=pmax(pip_EUR,pip_AFR)]
    for(a in fit$alpha) stopifnot(abs(sum(a)-1)<1e-8,all(a[!fit$availability$config_mask]==0))
    stopifnot(all(tab$pip_EUR[!x$variant_available[,"EUR"]]==0),all(tab$pip_AFR[!x$variant_available[,"AFR"]]==0),
      all(tab$pip_shared[rowSums(x$variant_available)<2]==0))
    ev <- fit$ELBO[is.finite(fit$ELBO)]
    # The upstream object does not expose a reliable explicit convergence flag.
    # Preserve ELBO and covariance diagnostics rather than claiming convergence.
    diagnostics <- list(converged=NA,elbo=fit$ELBO,finite_elbo_length=length(ev),
      min_elbo_step=if(length(ev)>1)min(diff(ev)) else NA_real_,
      covariance_min_eigenvalue=vapply(fit$V,function(V)min(eigen(V,symmetric=TRUE,only.values=TRUE)$values),numeric(1)))
    cs <- fit$cs$cs
    saveRDS(fit,file.path(attempt,"raw_fit.rds"))
  }
  vals <- unlist(tab[,.(pip_either,pip_shared,pip_EUR,pip_AFR,pip_unique)],use.names=FALSE)
  stopifnot(all(is.finite(vals)),all(vals>=-1e-12 & vals<=1+1e-12),identical(tab$SNP,x$union_ids))
  ci <- unique(as.integer(unlist(cs,use.names=FALSE)));stopifnot(all(ci %in% seq_len(p)))
  tab[,in_cs:=seq_len(.N) %in% ci]
  tab[,`:=`(available_EUR=x$variant_available[,"EUR"],available_AFR=x$variant_available[,"AFR"]) ]
  elapsed <- proc.time()[["elapsed"]]-t0
  result <- list(method=method,trait=x$trait,locus=x$locus,table=tab,cs=cs,
    settings=settings,diagnostics=diagnostics,warnings=warnings,elapsed_seconds=elapsed,
    input_file=input_file,input_md5=unname(tools::md5sum(input_file)),
    session=capture.output(sessionInfo()),raw_fit=file.path(attempt,"raw_fit.rds"),
    status=if(identical(diagnostics$converged,FALSE))"NONCONVERGED" else "COMPLETED")
  saveRDS(result,file.path(attempt,"result.rds"))
  fwrite(tab,file.path(attempt,"pip.tsv"),sep="\t")
  writeLines(warnings,file.path(attempt,"warnings.txt"))
  if(result$status=="NONCONVERGED") stop("SuSiE did not converge; raw result preserved for investigation")
  saveRDS(list(result=file.path(attempt,"result.rds")),file.path(outdir,"SUCCESS.rds"))
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
error=function(e){writeLines(c(conditionMessage(e),warnings),file.path(attempt,"ERROR.txt"));stop(e)})
