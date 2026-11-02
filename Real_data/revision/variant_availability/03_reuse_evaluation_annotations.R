#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Reuse saved downstream annotations; this is NOT the 186-feature fitting matrix.
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")), mustWork=TRUE)
src <- Sys.getenv("INFUSE_AOU_INPUTS", "/path/to/aou_inputs")
f <- file.path(src,"summary_res","res_all_mf_updated.RData")
env <- new.env(parent=baseenv()); load(f,envir=env)
dat <- as.data.table(env$res_all)
original_names <- names(dat)
stopifnot(!anyDuplicated(trimws(original_names)))
setnames(dat,trimws(original_names))
features <- c("missense","synonymous","utr_5","utr_3","promotor","CRE",
              "GERP4","GERP","brain_ind_eQTL","heart_ind_eQTL","artery_ind_eQTL",
              "non_synonymous_common")
identity <- c("PHENONAME","Region","SNP","CHR","POS","CHR_POS","ALT","REF")
stopifnot(all(c(identity,features) %in% names(dat)))
dat <- copy(dat[,c(identity,features),with=FALSE])
stopifnot(!anyDuplicated(dat[,.(PHENONAME,Region,SNP)]))
manifest <- fread(file.path(root,"qc","locus_manifest.tsv"))
qc <- list()
for(i in seq_len(nrow(manifest))){
  trait <- manifest$trait[i]; locus <- manifest$locus[i]
  ss <- rbindlist(lapply(c("eur","afr"),function(k)
    fread(file.path(src,paste0("results_",trait),"summary_data",paste0(k,"_loci_",locus,"_summ")))
  ))
  ss <- unique(ss[,.(SNP,CHR,POS,CHR_POS,ALT,REF)])
  stopifnot(!anyDuplicated(ss$SNP))
  a <- dat[PHENONAME==trait & Region==locus]
  stopifnot(setequal(a$SNP,ss$SNP))
  a <- a[match(ss$SNP,SNP)]
  for(nm in c("SNP","CHR","POS","CHR_POS","ALT","REF"))
    stopifnot(identical(as.character(a[[nm]]),as.character(ss[[nm]])))
  qc[[i]] <- data.table(trait=trait,locus=locus,n=nrow(a),
    missing_annotation_cells=sum(is.na(a[,..features])))
}
dat[,utr_comb:=fifelse(is.na(utr_3)|is.na(utr_5),NA_integer_,as.integer(utr_3+utr_5>0))]
dir.create(file.path(root,"inputs","evaluation_annotations"),showWarnings=FALSE,recursive=TRUE)
saveRDS(list(data=dat,source=f,source_md5=unname(tools::md5sum(f)),
  original_column_names=original_names,reused_columns=features,
  derived=list(utr_comb="as.integer(utr_3 + utr_5 > 0); preserve missing"),
  purpose="downstream enrichment only; not INFUSE baseline186 fitting annotations"),
  file.path(root,"inputs","evaluation_annotations","saved_evaluation_annotations.rds"))
fwrite(rbindlist(qc),file.path(root,"qc","saved_evaluation_annotation_alignment.tsv"),sep="\t")
print(rbindlist(qc)[,.(regions=.N,variant_region_records=sum(n),missing_annotation_cells=sum(missing_annotation_cells)),by=trait])
