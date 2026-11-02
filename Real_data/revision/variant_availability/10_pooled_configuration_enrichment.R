#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Pool BMI/DBP/SBP records BEFORE ranking; retain the approved background definition.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
out<-file.path(root,"summary/matched_size/pooled_BMI_BP")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
features<-c("non_synonymous_common","synonymous","promotor","utr_comb",
 "Conserved_LindbladToh_common","H3K4me3_Trynka_common")
methods<-c("MESuSiE_union","INFUSE_union")
targets<-c(Shared="pip_shared",Unique="pip_unique")
pf<-file.path(root,"summary/variant_pips.tsv.gz")
af<-file.path(root,"summary/evaluation_annotation_rows.tsv.gz")
p<-fread(pf,select=c("trait","locus","SNP","method",unname(targets)))[method %in% methods]
a<-fread(af,select=c("trait","locus","SNP",features))
stopifnot(!anyDuplicated(p[,.(trait,locus,SNP,method)]),!anyDuplicated(a[,.(trait,locus,SNP)]),
 setequal(p$trait,c("BMI","DBP","SBP")),nrow(p)==2*nrow(a))
d<-merge(p,a,by=c("trait","locus","SNP"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(d)==nrow(p))
for(f in features)stopifnot(all(is.na(d[[f]])|d[[f]] %in% c(0,1)))
metrics<-sizes<-selected<-list()
for(tar in names(targets))for(cut in c(.5,.9)){
 score<-targets[[tar]]
 inf<-d[method=="INFUSE_union"]
 stopifnot(all(is.finite(inf[[score]])))
 N<-sum(inf[[score]]>cut)
 for(m in methods){
  x<-copy(d[method==m & is.finite(get(score))])
  x[,ranking_score:=get(score)]
  setorderv(x,c("ranking_score","trait","locus","SNP"),c(-1,1,1,1))
  stopifnot(nrow(x)>=N)
  sel<-seq_len(nrow(x))<=N
  if(m=="INFUSE_union")stopifnot(identical(sel,x$ranking_score>cut))
  sizes[[length(sizes)+1L]]<-data.table(group="BMI_BP_pooled",target=tar,method=m,INFUSE_threshold=cut,N=N,
   n_universe=nrow(x),n_background=nrow(x)-N,selected_BMI=sum(sel & x$trait=="BMI"),
   selected_DBP=sum(sel & x$trait=="DBP"),selected_SBP=sum(sel & x$trait=="SBP"))
  selected[[length(selected)+1L]]<-x[sel,.(trait,locus,SNP,method,ranking_score)][,
   `:=`(target=tar,INFUSE_threshold=cut,rank=seq_len(.N),N=N)]
  for(f in features){
   ok<-!is.na(x[[f]]);ns<-sum(sel & ok);nb<-sum(!sel & ok)
   vs<-sum(x[[f]][sel & ok]);vb<-sum(x[[f]][!sel & ok])
   sf<-if(ns>0)vs/ns else NA_real_;bf<-if(nb>0)vb/nb else NA_real_
   status<-if(ns==0)"no_selected_with_annotation" else if(nb==0)"no_background_with_annotation" else if(bf==0)"zero_background_annotation_frequency" else "defined"
   metrics[[length(metrics)+1L]]<-data.table(group="BMI_BP_pooled",target=tar,method=m,INFUSE_threshold=cut,N=N,annotation=f,
    n_selected_annotated=ns,n_background_annotated=nb,selected_annotation_count=vs,background_annotation_count=vb,
    selected_frequency=sf,background_frequency=bf,enrichment=if(status=="defined")sf/bf else NA_real_,status=status)
  }
 }
}
z<-rbindlist(metrics);sz<-rbindlist(sizes)
stopifnot(nrow(z)==48L,nrow(sz)==8L)
fwrite(z,file.path(out,"functional_enrichment.tsv"),sep="\t",na="NA")
fwrite(sz,file.path(out,"selection_sizes.tsv"),sep="\t")
fwrite(rbindlist(selected),file.path(out,"selected_records.tsv"),sep="\t")
fwrite(data.table(file=c(pf,af),md5=unname(tools::md5sum(c(pf,af)))),file.path(out,"source_checksums.tsv"),sep="\t")
