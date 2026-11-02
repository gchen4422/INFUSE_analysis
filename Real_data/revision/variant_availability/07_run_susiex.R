#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages({library(data.table);library(snpStats)})
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
src <- Sys.getenv("INFUSE_AOU_INPUTS", "/path/to/aou_inputs")
task <- as.integer(commandArgs(TRUE)[1])
man <- fread(file.path(root,"inputs","frozen_manifests","region_manifest.tsv"))[task_id==task]
stopifnot(nrow(man)==1L,isTRUE(man$qc_pass))
input_dir <- file.path(root,"inputs","prepared",sprintf("%s_%03d",man$trait,man$locus))
stopifnot(file.exists(file.path(input_dir,"PASS")))
input_file <- file.path(input_dir,"input.rds"); x <- readRDS(input_file)
outdir <- file.path(root,"results","SuSiEx_union",sprintf("%s_%03d",man$trait,man$locus))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(outdir,"SUCCESS.rds"))) stop("Successful result exists; refusing rerun")
attempt <- file.path(outdir,paste0("attempt_",Sys.getenv("SLURM_JOB_ID",unset=Sys.getpid())))
dir.create(attempt,showWarnings=FALSE)
tryCatch({
  sst <- ref <- ld <- character()
  for(anc in c("EUR","AFR")){
    ss <- as.data.table(x$summary_stat_list[[anc]])
    # Fixed explicit column order: SNP CHR BP A1 A2 BETA SE P.
    s <- ss[,.(SNP,CHR,BP=POS,A1,A2,BETA=Beta,SE=Se,P=PVAL)]
    stopifnot(all(is.finite(s$P)),all(s$P>=0 & s$P<=1))
    sf <- file.path(attempt,paste0(anc,".sumstats.tsv"));fwrite(s,sf,sep="\t")
    sst <- c(sst,sf)
    ref <- c(ref,file.path(src,paste0("risk_loci_ld_",tolower(anc),"_",x$trait),paste0("loci_",x$locus)))
    ld <- c(ld,file.path(attempt,paste0(anc,"_reference_ld")))
  }
  n <- vapply(x$summary_stat_list,function(s)median(s$N),numeric(1))
  observed <- rbindlist(x$summary_stat_list)
  chr <- unique(observed$CHR);stopifnot(length(chr)==1L)
  cmdargs <- c(paste0("--sst_file=",paste(sst,collapse=",")),paste0("--ref_file=",paste(ref,collapse=",")),
    paste0("--ld_file=",paste(ld,collapse=",")),paste0("--n_gwas=",paste(n,collapse=",")),
    paste0("--out_dir=",attempt),"--out_name=SuSiEx_union",paste0("--chr=",chr),paste0("--bp=",min(observed$POS),",",max(observed$POS)),
    "--snp_col=1,1","--chr_col=2,2","--bp_col=3,3","--a1_col=4,4","--a2_col=5,5","--eff_col=6,6","--se_col=7,7","--pval_col=8,8",
    file.path(.infuse_repo, "Real_data/revision/variant_availability/plink_limited.sh"),"--keep-ambig=True","--mult-step=True","--maf=0.01",
    "--level=0.95","--min_purity=0.5","--pval_thresh=1e-05","--tol=0.0001","--n_sig=10","--max_iter=200","--threads=2")
  binary <- Sys.getenv("INFUSE_SUSIEX_BIN", unname(Sys.which("SuSiEx")))
  writeLines(c(binary,cmdargs),file.path(attempt,"command.txt"))
  t0<-proc.time()[["elapsed"]]
  exit <- system2(binary,args=shQuote(cmdargs),stdout=file.path(attempt,"stdout.log"),stderr=file.path(attempt,"stderr.log"))
  elapsed <- proc.time()[["elapsed"]]-t0
  stopifnot(exit==0L)
  snpf <- file.path(attempt,"SuSiEx_union.snp");csf <- file.path(attempt,"SuSiEx_union.cs")
  stopifnot(file.exists(snpf),file.exists(csf))
  csline <- readLines(csf,n=1L)
  if(any(grepl("FAIL",csline,fixed=TRUE))) stop("SuSiEx reported FAIL/nonconvergence")
  no_cs <- identical(trimws(csline),"NULL")
  raw <- fread(snpf,check.names=FALSE)
  stopifnot(all(c("SNP","BP") %in% names(raw)),!anyDuplicated(raw$SNP),all(raw$SNP %in% x$union_ids))
  pipcols <- grep("^PIP\\(CS[0-9]+\\)$",names(raw),value=TRUE)
  if(!no_cs && length(pipcols)==0L) stop("SuSiEx non-null output lacks expected PIP columns")
  if(length(pipcols)){
    pp <- as.matrix(raw[,..pipcols]);stopifnot(all(is.finite(pp)),all(pp>=0 & pp<=1))
    pip <- 1-apply(1-pp,1,prod)
  }else pip <- rep(0,nrow(raw))
  # Track effective ancestry-specific reference/SST intersections. These may be
  # smaller than the common intended universe due to SuSiEx's reference MAF filter.
  effective <- matrix(FALSE,length(x$union_ids),2,dimnames=list(x$union_ids,c("EUR","AFR")))
  ld_qc <- list()
  for(k in 1:2){
    anc <- c("EUR","AFR")[k]
    bim <- fread(paste0(ld[k],"_ref.bim"),header=FALSE)
    ss <- x$summary_stat_list[[anc]]
    ids <- ss$SNP[ss$SNP %in% bim$V2]
    effective[,k] <- x$union_ids %in% ids
    idx <- match(ids,bim$V2); sh <- match(ids,ss$SNP)
    same <- ss$A1[sh]==bim$V5[idx] & ss$A2[sh]==bim$V6[idx]
    swap <- ss$A1[sh]==bim$V6[idx] & ss$A2[sh]==bim$V5[idx]
    stopifnot(all(same|swap))
    con <- file(paste0(ld[k],".ld.bin"),"rb")
    v <- readBin(con,"numeric",n=nrow(bim)^2,size=4L);close(con)
    stopifnot(length(v)==nrow(bim)^2,all(is.finite(v)))
    R <- matrix(v,nrow(bim),nrow(bim))[idx,idx,drop=FALSE]
    sign <- ifelse(swap,-1,1);R <- R*outer(sign,sign)
    target <- x$R_mat_list[[anc]][ids,ids,drop=FALSE]
    delta <- max(abs(R-target))
    # Preserve the original pipeline's PLINK founder-only convention. The saved
    # R LD instead uses all reference individuals; compare against the actual
    # PLINK sample set independently, rather than relaxing the saved-LD gate.
    pl <- read.plink(paste0(ref[k],".bed"))
    fam <- fread(paste0(ref[k],".fam"),header=FALSE)
    founder <- fam$V3==0 & fam$V4==0
    source_bim <- fread(paste0(ref[k],".bim"),header=FALSE)
    gi <- match(ids,pl$map$snp.name)
    stopifnot(!anyNA(gi),identical(as.character(pl$map$snp.name),as.character(source_bim$V2)),sum(founder)>2)
    G <- as(pl$genotypes[founder,gi,drop=FALSE],"numeric")
    source_same <- ss$A1[sh]==source_bim$V5[gi] & ss$A2[sh]==source_bim$V6[gi]
    source_swap <- ss$A1[sh]==source_bim$V6[gi] & ss$A2[sh]==source_bim$V5[gi]
    stopifnot(all(source_same|source_swap))
    expected <- cor(G,use="pairwise.complete.obs")
    source_sign <- ifelse(source_swap,-1,1)
    expected <- expected*outer(source_sign,source_sign)
    founder_delta <- max(abs(R-expected))
    if(!is.finite(founder_delta)||founder_delta>1e-6)
      stop("SuSiEx LD differs from independently reconstructed founder LD for ",anc,": ",founder_delta)
    ld_qc[[anc]] <- data.table(ancestry=anc,n_intended=nrow(ss),n_effective=length(ids),
      n_reference_all=nrow(fam),n_reference_founders=sum(founder),
      founder_missing_genotype_calls=sum(is.na(G)),
      reference_vs_saved_LD_max_abs_delta=delta,reference_vs_founder_LD_max_abs_delta=founder_delta)
    rm(pl,G,expected,R,target,v);gc()
  }
  stopifnot(setequal(raw$SNP,x$union_ids[rowSums(effective)>0]))
  tab <- data.table(SNP=x$union_ids,pip_either=NA_real_,analyzed=rowSums(effective)>0,
    available_EUR=x$variant_available[,1],available_AFR=x$variant_available[,2],
    effective_EUR=effective[,1],effective_AFR=effective[,2],in_cs=FALSE)
  tab$pip_either[match(raw$SNP,tab$SNP)] <- pip
  cs <- list()
  if(!no_cs){
    cst <- fread(csf)
    stopifnot(all(c("CS_ID","SNP") %in% names(cst)),all(cst$SNP %in% raw$SNP))
    cs <- lapply(split(cst$SNP,cst$CS_ID),function(s)match(s,x$union_ids))
    tab$in_cs[unique(unlist(cs))] <- TRUE
  }
  fwrite(rbindlist(ld_qc),file.path(attempt,"reference_ld_audit.tsv"),sep="\t")
  fwrite(tab,file.path(attempt,"pip.tsv"),sep="\t")
  result <- list(method="SuSiEx_union",trait=x$trait,locus=x$locus,table=tab,cs=cs,
    status=if(no_cs)"COMPLETED_NO_CS" else "COMPLETED",elapsed_seconds=elapsed,
    command=cmdargs,binary=binary,binary_md5=unname(tools::md5sum(binary)),input_file=input_file,
    input_md5=unname(tools::md5sum(input_file)),diagnostics=list(reference_ld=rbindlist(ld_qc)),
    pip_definition="1-product(1-per-CS PIP); unanalysed SNPs remain NA, not failed-fit zero",
    reference_policy="Original PLINK founder-only LD and MAF; other methods use saved all-reference-sample LD")
  saveRDS(result,file.path(attempt,"result.rds"))
  saveRDS(list(result=file.path(attempt,"result.rds")),file.path(outdir,"SUCCESS.rds"))
},error=function(e){writeLines(conditionMessage(e),file.path(attempt,"ERROR.txt"));stop(e)})
