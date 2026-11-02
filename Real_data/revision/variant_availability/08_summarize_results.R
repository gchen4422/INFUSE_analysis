#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Summarize saved inputs and fit results; write outputs under summary/.
# Enrichment point estimates only: no independence-based CI or p-value is claimed.
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")), mustWork=TRUE)
args <- commandArgs(TRUE)
stopifnot(all(args %in% "--partial"))
partial <- "--partial" %in% args
outdir <- file.path(root,"summary")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
methods <- c("INFUSE_union","MESuSiE_union","SuSiEx_union","SuSiE_posthoc")
accepted_nonconverged <- function(r) {
  identical(r$status,"ACCEPTED_NONCONVERGED_ITER100") &&
    identical(r$method,"SuSiE_posthoc") && identical(r$trait,"BMI") &&
    identical(as.integer(r$locus),182L) && identical(r$diagnostics$converged,FALSE) &&
    identical(as.integer(r$settings$max_iter),100L) && isTRUE(r$acceptance$user_authorized)
}
allowed_status <- c("COMPLETED","COMPLETED_NO_CS","ACCEPTED_NONCONVERGED_ITER100")
thresholds <- round(seq(.5,.9,.05),2)
features <- c("non_synonymous_common","missense","synonymous","promotor","utr_comb","CRE",
              "heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")
manifest_file <- file.path(root,"inputs","frozen_manifests","region_manifest.tsv")
man <- fread(manifest_file)
stopifnot(nrow(man)==289L,!anyDuplicated(man[,.(trait,locus)]),
  !anyDuplicated(man$task_id),setequal(man$trait,c("BMI","DBP","SBP")),
  nrow(man[trait=="BMI"])==204L,nrow(man[trait=="DBP"])==39L,nrow(man[trait=="SBP"])==46L)
inside <- function(path) {
  p <- normalizePath(path,mustWork=TRUE)
  if(!startsWith(p,paste0(root,"/"))) stop("Result/input path escapes task directory: ",p)
  p
}
records <- list(); status <- list(); z <- 0L
for(i in seq_len(nrow(man))) for(method in methods) {
  z <- z+1L
  d <- file.path(root,"results",method,sprintf("%s_%03d",man$trait[i],man$locus[i]))
  success <- file.path(d,"SUCCESS.rds")
  row <- data.table(task_id=man$task_id[i],trait=man$trait[i],locus=man$locus[i],method=method,
    qc_pass=isTRUE(man$qc_pass[i]),success_file=success,result_file=NA_character_,
    status="MISSING",message=NA_character_)
  if(file.exists(success)) {
    tryCatch({
      f <- inside(readRDS(success)$result); r <- readRDS(f)
      stopifnot(identical(r$method,method),identical(as.character(r$trait),as.character(man$trait[i])),
                identical(as.integer(r$locus),as.integer(man$locus[i])),
                r$status %in% allowed_status,
                r$status!="ACCEPTED_NONCONVERGED_ITER100" || accepted_nonconverged(r))
      row[,`:=`(result_file=f,status=r$status)]
      records[[paste(man$trait[i],man$locus[i],method,sep="|")]] <- r
    },error=function(e) {row[,`:=`(status="INVALID_SUCCESS",message=conditionMessage(e))]})
  } else {
    errors <- list.files(d,pattern="^ERROR[.]txt$",recursive=TRUE,full.names=TRUE)
    if(length(errors)) row[,`:=`(status="FAILED_ATTEMPT_NO_SUCCESS",message=paste(errors,collapse=";"))]
  }
  status[[z]] <- row
}
status <- rbindlist(status)
fwrite(status,file.path(outdir,"completion_status.tsv"),sep="\t",na="NA")
fwrite(status[,.(regions=.N),by=.(trait,method,status)],file.path(outdir,"completion_counts.tsv"),sep="\t")
if(partial) {
  cat("Partial mode: status files only; no scientific summaries generated. Existing full summaries, if any, are not refreshed.\n")
  quit(save="no",status=0L)
}
if(!all(status$qc_pass) || !all(status$status %in% allowed_status))
  stop("Full summary requires eligible results for all 289 regions and 4 methods; see completion_status.tsv")
eval_file <- file.path(root,"inputs","evaluation_annotations","saved_evaluation_annotations.rds")
evaluation <- readRDS(eval_file)
eval <- copy(as.data.table(evaluation$data))
setnames(eval,c("PHENONAME","Region"),c("trait","locus"))
stopifnot(all(c("trait","locus","SNP","CHR","POS","ALT","REF",features) %in% names(eval)),
          !anyDuplicated(eval[,.(trait,locus,SNP)]))
for(nm in features) {
  v <- eval[[nm]]
  stopifnot(is.numeric(v)||is.logical(v),all(is.na(v)|v %in% c(0,1)))
}
pip_parts <- region_parts <- annotation_parts <- list()
validation <- list(); z <- 0L
for(i in seq_len(nrow(man))) {
  tr <- man$trait[i]; lo <- man$locus[i]
  input_file <- inside(file.path(root,"inputs","prepared",sprintf("%s_%03d",tr,lo),"input.rds"))
  x <- readRDS(input_file); input_md5 <- unname(tools::md5sum(input_file))
  stopifnot(identical(as.character(x$trait),as.character(tr)),as.integer(x$locus)==as.integer(lo),
    !anyDuplicated(x$union_ids),length(x$union_ids)==man$retained[i],
    identical(rownames(x$variant_available),x$union_ids),
    identical(colnames(x$variant_available),c("EUR","AFR")),all(rowSums(x$variant_available)>0))
  a <- eval[trait==tr & locus==lo]
  hit <- match(x$union_ids,a$SNP);stopifnot(!anyNA(hit));a <- a[hit]
  ref <- as.data.table(x$allele_ref)[match(x$union_ids,SNP)]
  stopifnot(identical(as.character(a$ALT),as.character(ref$A1)),
            identical(as.character(a$REF),as.character(ref$A2)))
  for(anc in c("EUR","AFR")) {
    ss <- x$summary_stat_list[[anc]]; jj <- match(ss$SNP,a$SNP)
    stopifnot(!anyNA(jj),identical(as.character(ss$CHR),as.character(a$CHR[jj])),
      identical(as.character(ss$POS),as.character(a$POS[jj])),
      identical(as.character(ss$A1),as.character(a$ALT[jj])),
      identical(as.character(ss$A2),as.character(a$REF[jj])))
  }
  a[,availability_pattern:=fcase(x$variant_available[,1]&x$variant_available[,2],"shared",
       x$variant_available[,1],"EUR-only",default="AFR-only")]
  annotation_parts[[i]] <- a[,c("trait","locus","SNP","CHR","POS","ALT","REF","availability_pattern",features),with=FALSE]
  for(method in methods) {
    z <- z+1L;r <- records[[paste(tr,lo,method,sep="|")]];tab <- copy(as.data.table(r$table))
    stopifnot(identical(inside(r$input_file),input_file),identical(r$input_md5,input_md5),
      identical(as.character(tab$SNP),x$union_ids),!anyDuplicated(tab$SNP),
      identical(as.logical(tab$available_EUR),as.logical(x$variant_available[,1])),
      identical(as.logical(tab$available_AFR),as.logical(x$variant_available[,2])),
      is.numeric(r$elapsed_seconds),length(r$elapsed_seconds)==1L,is.finite(r$elapsed_seconds),r$elapsed_seconds>=0)
    targets <- if(method=="SuSiEx_union") "pip_either" else c("pip_either","pip_shared","pip_unique","pip_EUR","pip_AFR")
    stopifnot(all(targets %in% names(tab)))
    for(nm in targets) {
      # Match the runner's floating-point tolerance; preserve raw values rather
      # than clipping or renormalizing any method's posterior probabilities.
      v <- tab[[nm]];stopifnot(is.numeric(v),all(is.na(v)|(is.finite(v)&v>=-1e-12&v<=1+1e-12)))
      if(method!="SuSiEx_union") stopifnot(all(is.finite(v)))
    }
    if(method=="SuSiEx_union") {
      stopifnot(!anyNA(tab$analyzed),identical(is.finite(tab$pip_either),as.logical(tab$analyzed)),
        identical(as.logical(tab$analyzed),as.logical(tab$effective_EUR|tab$effective_AFR)),
        !any(tab$effective_EUR & !tab$available_EUR),!any(tab$effective_AFR & !tab$available_AFR))
    } else {
      tab[,`:=`(analyzed=TRUE,effective_EUR=available_EUR,effective_AFR=available_AFR)]
      stopifnot(all(tab$pip_EUR[!tab$available_EUR]==0),all(tab$pip_AFR[!tab$available_AFR]==0),
                all(tab$pip_shared[!tab$available_EUR|!tab$available_AFR]==0))
    }
    stopifnot(!identical(r$diagnostics$converged,FALSE) || accepted_nonconverged(r),is.null(r$cs)||is.list(r$cs),all(vapply(r$cs,function(cs)
      is.numeric(cs)&&!anyNA(cs)&&all(cs==as.integer(cs))&&!anyDuplicated(cs)&&all(cs %in% seq_len(nrow(tab))),logical(1))))
    cs_ids <- unique(as.integer(unlist(r$cs,use.names=FALSE)))
    stopifnot(identical(as.logical(tab$in_cs),seq_len(nrow(tab)) %in% cs_ids),all(tab$analyzed[cs_ids]))
    for(nm in setdiff(c("pip_shared","pip_unique","pip_EUR","pip_AFR"),names(tab))) tab[,(nm):=NA_real_]
    tab[,`:=`(trait=tr,locus=lo,method=method,availability_pattern=a$availability_pattern)]
    pip_parts[[z]] <- tab[,.(trait,locus,SNP,method,availability_pattern,available_EUR,available_AFR,
      analyzed,effective_EUR,effective_AFR,pip_either,pip_shared,pip_unique,pip_EUR,pip_AFR,in_cs)]
    region_parts[[z]] <- data.table(trait=tr,locus=lo,method=method,n_intended=nrow(tab),
      n_effective=sum(tab$analyzed),n_unanalyzed=sum(!tab$analyzed),
      n_EUR_effective=sum(tab$effective_EUR),n_AFR_effective=sum(tab$effective_AFR),
      n_cs=length(r$cs),distinct_cs_snps=length(cs_ids),elapsed_seconds=r$elapsed_seconds,
      converged=if(is.null(r$diagnostics$converged))NA else r$diagnostics$converged,
      warning_count=length(r$warnings),status=r$status)
    validation[[z]] <- data.table(trait=tr,locus=lo,method=method,input_md5=input_md5,check="PASS")
  }
  rm(x);gc(verbose=FALSE)
}
annotation <- rbindlist(annotation_parts);pips <- rbindlist(pip_parts);region <- rbindlist(region_parts)
overlaps <- annotation[,.(n_loci=uniqueN(locus)),by=.(trait,SNP)][n_loci>1]
if(nrow(overlaps)) {
  fwrite(overlaps,file.path(outdir,"overlapping_trait_snp_REVIEW_REQUIRED.tsv"),sep="\t")
  stop("Variants occur in multiple loci within a trait; do not count duplicated discoveries. Review overlap policy first.")
}
stopifnot(!anyDuplicated(pips[,.(trait,locus,SNP,method)]))
common <- pips[,.(common_four_method=.N==length(methods)&&all(is.finite(pip_either))),by=.(trait,locus,SNP)]
pips <- merge(pips,common,by=c("trait","locus","SNP"),sort=FALSE)
joined <- merge(pips,annotation[,c("trait","locus","SNP",features),with=FALSE],by=c("trait","locus","SNP"),sort=FALSE)
stopifnot(nrow(joined)==nrow(pips))
target_cols <- c(Either="pip_either",Shared="pip_shared",Unique="pip_unique")
enrich <- counts <- list();ie <- ic <- 0L
for(tr in c("BMI","DBP","SBP")) for(method in methods) for(target in names(target_cols)) {
  if(method=="SuSiEx_union" && target!="Either") next
  current_method <- method
  d <- joined[trait==tr & method==current_method]
  score <- target_cols[[target]]
  for(threshold in thresholds) {
    dd <- copy(d);dd[,selected:=is.finite(get(score)) & get(score)>threshold]
    ic <- ic+1L
    counts[[ic]] <- dd[,.(n_intended=.N,n_finite_score=sum(is.finite(get(score))),
      n_selected=sum(selected),n_selected_in_cs=sum(selected&in_cs)),by=.(trait,locus,method,availability_pattern)][,
        `:=`(target=target,threshold=threshold)]
    for(universe in c("method_effective","common_four_method")) {
      keep <- is.finite(dd[[score]])
      if(universe=="common_four_method") keep <- keep & dd$common_four_method
      v <- dd[keep]; sel <- v$selected
      for(ann in features) {
        av <- v[[ann]];observed <- !is.na(av)
        ns <- sum(sel & observed);nb <- sum(!sel & observed)
        as <- sum(av[sel & observed]);ab <- sum(av[!sel & observed])
        fs <- if(ns>0)as/ns else NA_real_;fb <- if(nb>0)ab/nb else NA_real_
        reason <- if(ns==0)"no_selected_with_annotation" else if(nb==0)"no_background_with_annotation" else if(fb==0)"zero_background_annotation_frequency" else "defined"
        value <- if(reason=="defined")fs/fb else NA_real_
        ie <- ie+1L
        enrich[[ie]] <- data.table(trait=tr,method=method,target=target,threshold=threshold,universe=universe,
          annotation=ann,n_universe=nrow(v),n_selected_total=sum(sel),n_background_total=sum(!sel),
          n_selected_annotated=ns,n_background_annotated=nb,
          n_selected_missing=sum(sel & !observed),n_background_missing=sum(!sel & !observed),
          selected_annotation_count=as,background_annotation_count=ab,
          selected_frequency=fs,background_frequency=fb,enrichment=value,status=reason)
      }
    }
  }
}
enrichment <- rbindlist(enrich);threshold_counts <- rbindlist(counts)
overlap <- list();io <- 0L
for(tr in c("BMI","DBP","SBP")) for(threshold in thresholds) {
  d <- pips[trait==tr & common_four_method==TRUE]
  pairs <- combn(methods,2,simplify=FALSE)
  for(pair in pairs) {
    sa <- d[method==pair[1] & pip_either>threshold,SNP]
    sb <- d[method==pair[2] & pip_either>threshold,SNP]
    io <- io+1L;un <- union(sa,sb)
    overlap[[io]] <- data.table(trait=tr,threshold=threshold,method_1=pair[1],method_2=pair[2],
      universe="common_four_method",n_universe=uniqueN(d$SNP),selected_1=length(sa),selected_2=length(sb),
      overlap=length(intersect(sa,sb)),selected_union=length(un),jaccard=if(length(un))length(intersect(sa,sb))/length(un) else NA_real_)
  }
}
trait_metrics <- region[,.(regions=.N,regions_with_cs=sum(n_cs>0),median_locus_cs_snps=as.numeric(median(distinct_cs_snps)),
  median_locus_cs_snps_nonempty=if(any(distinct_cs_snps>0))as.numeric(median(distinct_cs_snps[distinct_cs_snps>0])) else NA_real_,
  total_cs=sum(n_cs),mean_elapsed_seconds=mean(elapsed_seconds),median_elapsed_seconds=median(elapsed_seconds),
  total_elapsed_seconds=sum(elapsed_seconds),convergence_unknown=sum(is.na(converged)),
  nonconverged=sum(!converged,na.rm=TRUE)),by=.(trait,method)]
feature_definitions <- data.table(annotation=features,definition=c(
  "Saved baselineLF non_synonymous_common; overlaps the INFUSE fitting annotation set",
  "Saved script-7 category combining LOFTEE high/low-confidence loss-of-function and VEP missense calls; distinct from baseline non_synonymous_common",
  "Saved synonymous category, excluding variants in the script-7 missense/LoF category",
  "Saved promoter category after hierarchical exclusion of coding and UTR categories",
  "Union of saved utr_3 and utr_5; missing if either constituent is missing",
  "Saved CRE category from H3K27ac and DHS overlap after hierarchical exclusions",
  "Saved heart independent eQTL indicator","Saved artery independent eQTL indicator",
  "Saved brain independent eQTL indicator"))
tables <- list(region_metrics=region,trait_metrics=trait_metrics,region_threshold_counts=threshold_counts,
  trait_threshold_counts=threshold_counts[,.(n_intended=sum(n_intended),n_finite_score=sum(n_finite_score),
    n_selected=sum(n_selected),n_selected_in_cs=sum(n_selected_in_cs)),by=.(trait,method,target,threshold,availability_pattern)],
  functional_enrichment=enrichment,method_overlap=rbindlist(overlap),input_validation=rbindlist(validation),
  evaluation_feature_definitions=feature_definitions)
for(nm in names(tables)) fwrite(tables[[nm]],file.path(outdir,paste0(nm,".tsv")),sep="\t",na="NA")
fwrite(pips,file.path(outdir,"variant_pips.tsv.gz"),sep="\t",na="NA")
fwrite(annotation,file.path(outdir,"evaluation_annotation_rows.tsv.gz"),sep="\t",na="NA")
saveRDS(list(tables=tables,pips=pips,evaluation_annotations=annotation,
  provenance=list(manifest=manifest_file,manifest_md5=unname(tools::md5sum(manifest_file)),
    evaluation_file=eval_file,evaluation_md5=unname(tools::md5sum(eval_file)),evaluation_original_source=evaluation$source,
    methods=methods,thresholds=thresholds,features=features,
    enrichment="Selected annotation frequency / nonselected annotation frequency; annotation-specific complete cases; undefined ratios NA; point estimates only",
    universe="Primary finite method-specific score universe; sensitivity common finite Either PIP across all four methods",
    unique="Stored scores; INFUSE/MESuSiE max exclusive-configuration PIPs; SuSiE strict 0.5 post hoc gate",
    timestamp=as.character(Sys.time()),session=capture.output(sessionInfo()))),file.path(outdir,"summary_data.rds"))
