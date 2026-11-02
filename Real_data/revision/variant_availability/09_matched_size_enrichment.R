#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
out<-file.path(root,"summary/matched_size");dir.create(out,recursive=TRUE,showWarnings=FALSE)
methods<-c("INFUSE_union","MESuSiE_union","SuSiEx_union","SuSiE_posthoc")
thresholds<-round(seq(.5,.9,.05),2)
features<-c("non_synonymous_common","missense","synonymous","promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL","Conserved_LindbladToh_common","H3K4me3_Trynka_common")
source(file.path(.infuse_repo, "Real_data/revision/variant_availability/annotate_enrichment_panels.R"),local=new.env())
p<-fread(file.path(root,"summary/variant_pips.tsv.gz"),select=c("trait","locus","SNP","method","pip_either","pip_shared","pip_unique"))
a<-fread(file.path(root,"summary/evaluation_annotation_rows.tsv.gz"))
stopifnot(!anyDuplicated(p[,.(trait,locus,SNP,method)]),!anyDuplicated(a[,.(trait,locus,SNP)]),setequal(p$method,methods))
for(f in features)stopifnot(all(is.na(a[[f]])|a[[f]] %in% c(0,1)))
# BP is pooled before ranking; cross-trait SNP records are not deduplicated.
p[,group:=ifelse(trait=="BMI","BMI","Blood_pressure")]
d<-merge(p,a[,c("trait","locus","SNP",features),with=FALSE],by=c("trait","locus","SNP"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(d)==nrow(p))
rank_parts<-selected_parts<-metrics<-sizes<-list()
targets<-c(Either="pip_either",Shared="pip_shared",Unique="pip_unique")
for(g in c("BMI","Blood_pressure"))for(tar in names(targets)){
  score<-targets[[tar]]
  inf<-d[group==g & method=="INFUSE_union"]
  stopifnot(all(is.finite(inf[[score]])))
  NN<-vapply(thresholds,function(t)sum(inf[[score]]>t),integer(1))
  for(m in if(tar=="Either")methods else setdiff(methods,"SuSiEx_union")){
    x<-copy(d[group==g & method==m & is.finite(get(score))])
    x[,ranking_score:=get(score)]
    setorderv(x,c("ranking_score","trait","locus","SNP"),c(-1,1,1,1))
    x[,rank:=seq_len(.N)];stopifnot(nrow(x)>=max(NN))
    rank_parts[[length(rank_parts)+1L]]<-x[,.(group,method,trait,locus,SNP,ranking_score,rank)][,target:=tar]
    for(i in seq_along(thresholds)){
      N<-NN[i];t<-thresholds[i];sel<-x$rank<=N;stopifnot(sum(sel)==N)
      if(m=="INFUSE_union")stopifnot(identical(sel,x$ranking_score>t))
      sizes[[length(sizes)+1L]]<-data.table(group=g,target=tar,method=m,INFUSE_threshold=t,N=N,
        n_universe=nrow(x),n_background=nrow(x)-N,
        n_unrankable=sum(d$group==g & d$method==m & !is.finite(d[[score]])),
        selected_BMI=sum(sel & x$trait=="BMI"),selected_DBP=sum(sel & x$trait=="DBP"),selected_SBP=sum(sel & x$trait=="SBP"))
      selected_parts[[length(selected_parts)+1L]]<-x[sel,.(group,method,trait,locus,SNP,ranking_score,rank)][,`:=`(target=tar,INFUSE_threshold=t,N=N)]
      for(f in features){
        av<-x[[f]];ok<-!is.na(av);ns<-sum(sel & ok);nb<-sum(!sel & ok)
        vs<-sum(av[sel & ok]);vb<-sum(av[!sel & ok])
        fs<-if(ns>0)vs/ns else NA_real_;fb<-if(nb>0)vb/nb else NA_real_
        status<-if(ns==0)"no_selected_with_annotation" else if(nb==0)"no_background_with_annotation" else if(fb==0)"zero_background_annotation_frequency" else "defined"
        metrics[[length(metrics)+1L]]<-data.table(group=g,target=tar,method=m,INFUSE_threshold=t,N=N,annotation=f,
          n_universe=nrow(x),n_selected=N,n_background=nrow(x)-N,n_selected_annotated=ns,n_background_annotated=nb,
          selected_annotation_count=vs,background_annotation_count=vb,n_selected_missing=N-ns,n_background_missing=nrow(x)-N-nb,
          selected_frequency=fs,background_frequency=fb,enrichment=if(status=="defined")fs/fb else NA_real_,status=status)
      }
    }
  }
}
sz<-rbindlist(sizes);z<-rbindlist(metrics)
stopifnot(all(sz[,.(n_sizes=uniqueN(N)),by=.(group,target,INFUSE_threshold)]$n_sizes==1L),nrow(z)==2*10*9*length(features))
fwrite(sz,file.path(out,"selection_sizes.tsv"),sep="\t")
fwrite(z,file.path(out,"functional_enrichment.tsv"),sep="\t",na="NA")
fwrite(rbindlist(selected_parts),file.path(out,"selected_records.tsv.gz"),sep="\t")
fwrite(rbindlist(rank_parts),file.path(out,"ranked_universe.tsv.gz"),sep="\t")
