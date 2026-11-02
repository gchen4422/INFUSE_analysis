#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages({library(data.table);library(INFUSE)})
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
task <- as.integer(commandArgs(TRUE)[1])
manifest <- fread(file.path(root,"inputs","frozen_manifests","region_manifest.tsv"))
m <- manifest[task_id==task]
stopifnot(nrow(m)==1L,isTRUE(m$qc_pass))
outdir <- file.path(root,"inputs","prepared",sprintf("%s_%03d",m$trait,m$locus))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
outfile <- file.path(outdir,"input.rds")
if(file.exists(outfile)) stop("Prepared input already exists; refusing overwrite: ",outfile)
x <- readRDS(m$audit_file)
decisions <- fread(file.path(root,"inputs","frozen_manifests",sprintf("%s_%03d",m$trait,m$locus),"variant_decisions.tsv"))
ids <- decisions[included==TRUE,SNP]
stopifnot(length(ids)>0L,!anyDuplicated(ids))
annotation <- fread(file.path(root,"inputs","baseline_186",paste0("chr",m$CHR,".tsv")))
features <- readLines(file.path(root,"inputs","baseline_186","annotation_names.txt"))
stopifnot(length(features)==186L,!"GTEx_eQTL_MaxCPP_common" %in% features,!anyDuplicated(annotation$SNP))
hit <- match(ids,annotation$SNP);stopifnot(!anyNA(hit))
A <- as.matrix(annotation[hit,..features]);rownames(A)<-ids
stopifnot(is.numeric(A),all(is.finite(A)))
for(anc in c("EUR","AFR")){
  ss <- as.data.table(x$summary_stat_list[[anc]])
  setnames(ss,c("BETA","SE","ALT","REF"),c("Beta","Se","A1","A2"))
  ss <- ss[match(ids[ids %in% ss$SNP],SNP)]
  stopifnot(nrow(ss)>0L,!anyNA(ss$SNP))
  x$summary_stat_list[[anc]] <- as.data.frame(ss)
  x$R_mat_list[[anc]] <- x$R_mat_list[[anc]][ss$SNP,ss$SNP,drop=FALSE]
}
ref <- as.data.frame(decisions[match(ids,SNP),.(SNP,A1=ALT,A2=REF)])
prep <- INFUSE::prepare_infuse_union(x$summary_stat_list,x$R_mat_list,
  annot=A,union_ids=ids,allele_ref=ref,ld_allele_coding="ancestry_effect",verbose=TRUE)
stopifnot(identical(prep$union_ids,ids),all(rowSums(prep$variant_available)>0),
  identical(prep$variant_available[,"EUR"],setNames(ids %in% x$summary_stat_list$EUR$SNP,ids)),
  identical(prep$variant_available[,"AFR"],setNames(ids %in% x$summary_stat_list$AFR$SNP,ids)),
  identical(as.numeric(prep$annot),as.numeric(A)))
# Store observed matrices only to avoid duplicating padded union matrices on disk.
saveRDS(list(trait=m$trait,locus=m$locus,task_id=task,union_ids=ids,
  summary_stat_list=x$summary_stat_list,R_mat_list=x$R_mat_list,
  variant_available=prep$variant_available,variant_map=prep$variant_map,
  annot=A,allele_ref=ref,decisions=decisions,source_audit=m$audit_file,
  annotation_excluded="GTEx_eQTL_MaxCPP_common",fitting_annotation_count=186L,
  input_md5=unname(tools::md5sum(m$audit_file))),outfile,compress=FALSE)
fwrite(prep$variant_map,file.path(outdir,"variant_map.tsv"),sep="\t")
writeLines("PASS",file.path(outdir,"PASS"))
