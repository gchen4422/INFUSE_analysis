#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
status <- fread(file.path(root,"qc","annotation_match_status.tsv"))
manifest <- fread(file.path(root,"qc","locus_manifest.tsv"))
evaluation <- readRDS(file.path(root,"inputs","evaluation_annotations","saved_evaluation_annotations.rds"))
e <- as.data.table(evaluation$data)
setnames(e,c("PHENONAME","Region"),c("trait","locus"))
stopifnot(!anyDuplicated(status[,.(trait,locus,SNP)]),!anyDuplicated(e[,.(trait,locus,SNP)]))
status <- merge(status,e[,.(trait,locus,SNP,POS,ALT,REF)],by=c("trait","locus","SNP"),all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(status$POS),nrow(status)==sum(manifest$n_shared+manifest$n_EUR_only+manifest$n_AFR_only))
status[,availability_pattern:=fcase(availability=="afr+eur","shared",availability=="eur","EUR-only",availability=="afr","AFR-only",default=NA_character_)]
stopifnot(!anyNA(status$availability_pattern))
status[,included:=annotation_status=="matched"]
status[,exclusion_reason:=fcase(annotation_status=="absent_rsid","absent_from_baselineLF",annotation_status=="allele_mismatch","baselineLF_allele_mismatch",annotation_status=="ambiguous_annotation","ambiguous_baselineLF_annotation",default="none")]
setorder(status,trait,locus,CHR,POS,SNP)
out <- file.path(root,"inputs","frozen_manifests")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(status,file.path(out,"all_variant_decisions.tsv.gz"),sep="\t")
fwrite(status[included==FALSE],file.path(out,"excluded_variants.tsv.gz"),sep="\t")
fwrite(status[,.(n=.N),by=.(trait,availability_pattern,included,exclusion_reason)],file.path(out,"counts_by_trait_availability_reason.tsv"),sep="\t")
fwrite(status[,.(n=.N),by=.(trait,locus,availability_pattern,included,exclusion_reason)],file.path(out,"counts_by_locus_availability_reason.tsv"),sep="\t")
reg <- status[,.(before=.N,retained=sum(included),excluded=sum(!included),
  shared_before=sum(availability_pattern=="shared"),shared_retained=sum(included & availability_pattern=="shared"),
  EUR_only_before=sum(availability_pattern=="EUR-only"),EUR_only_retained=sum(included & availability_pattern=="EUR-only"),
  AFR_only_before=sum(availability_pattern=="AFR-only"),AFR_only_retained=sum(included & availability_pattern=="AFR-only")),by=.(trait,locus)]
manifest <- merge(manifest,reg,by=c("trait","locus"),sort=FALSE)
manifest[,audit_file:=vapply(seq_len(.N),function(i){
  candidates <- list.files(file.path(root,"qc_v2",sprintf("%s_%03d",trait[i],locus[i])),pattern="^PASS$",recursive=TRUE,full.names=TRUE)
  if(length(candidates)==1L && !file.exists(file.path(dirname(candidates),"ERROR.txt"))) file.path(dirname(candidates),"audited_observed.rds") else NA_character_
},character(1))]
manifest[,qc_pass:=!is.na(audit_file)]
setorder(manifest,task_id)
fwrite(manifest,file.path(out,"region_manifest.tsv"),sep="\t")
for(i in seq_len(nrow(manifest))){
  tr <- manifest$trait[i]; lo <- manifest$locus[i]
  d <- file.path(out,sprintf("%s_%03d",tr,lo));dir.create(d,showWarnings=FALSE)
  fwrite(status[trait==tr & locus==lo],file.path(d,"variant_decisions.tsv"),sep="\t")
}
fwrite(data.table(file=c(file.path(root,"qc","annotation_match_status.tsv"),evaluation$source),
  md5=unname(tools::md5sum(c(file.path(root,"qc","annotation_match_status.tsv"),evaluation$source)))),file.path(out,"source_checksums.tsv"),sep="\t")
print(status[,.(n=.N),by=.(trait,availability_pattern,included,exclusion_reason)])
if(!all(manifest$qc_pass)) cat("Pending QC regions cannot be prepared/fitted; regenerate region_manifest after all QC terminates.\n")
