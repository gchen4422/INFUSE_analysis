#!/usr/bin/env Rscript
# Figure labels use the manuscript name; stored result keys remain compatible.
availability_method_labels <- function(x) {
  sub("^MESuSiE(?=([ ~]|$))", "INFUSE\u2212", as.character(x), perl = TRUE)
}

.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1L]))
OUT <- file.path(.infuse_work, "INFUSE_variant_availability/simulation_native_assv/multisusie_union_test_20261020")
NATIVE <- dirname(OUT)
ROOT <- dirname(NATIVE)
REVISION <- dirname(ROOT)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
# Use base R means like the original dplyr summaries. GForce's alternate
# accumulation can turn an exact reporting tie (e.g. 0.335) into 0.33499999999999996.
options(datatable.optimize=1L)
base_counts <- c(2L,3L,5L)
grid <- CJ(locus=1:100,base_index=1:3,h2_index=1:2,assv_each=0:2)
grid[,path:=file.path(OUT,"results",paste0("base_",base_counts[base_index]),
  paste0("h2_",h2_index),paste0("assv_each_",assv_each),sprintf("locus_%03d.rds",locus))]
grid[,exists:=file.exists(path)]
dir.create(file.path(OUT,"tables"),showWarnings=FALSE)
fwrite(grid,file.path(OUT,"tables","completion_manifest.tsv"),sep="\t")
if(!all(grid$exists))stop("Incomplete grid: ",sum(grid$exists),"/1800. No partial figures made.")
cache_path<-file.path(OUT,"metric_cache.rds")
if(file.exists(cache_path) && file.info(cache_path)$mtime>=max(file.info(grid$path)$mtime)) {
  cached<-readRDS(cache_path)
  metrics<-cached$metrics;diagnostics<-cached$diagnostics
  rm(cached)
  cat("Using cached cell metrics\n")
} else {
metrics <- setNames(lapply(c("cs","roc","calibration","causal"),function(x)vector("list",nrow(grid))),
                    c("cs","roc","calibration","causal"))
diag <- vector("list",nrow(grid))
for(i in seq_len(nrow(grid))) {
  f <- readRDS(grid$path[i]);m<-f$metadata
  stopifnot(m$locus==grid$locus[i],m$base_index==grid$base_index[i],m$h2_index==grid$h2_index[i],
    m$assv_each==grid$assv_each[i],m$variant_set=="union",all(is.finite(f$pip)))
  for(k in names(metrics))metrics[[k]][[i]]<-f$metrics[[k]]
  diag[[i]]<-data.table(locus=m$locus,base_index=m$base_index,h2_index=m$h2_index,
    assv_each=m$assv_each,n_union=m$n_union,n_common=m$n_common,converged=f$converged,
    niter=f$niter,runtime_minutes=f$runtime_minutes,source_sha256=m$package_source_sha256,
    max_alpha_sum_error=max(abs(rowSums(f$alpha)-1)),last_elbo=tail(f$elbo,1))
  if(i%%100L==0L)cat("Loaded ",i,"/1800\n",sep="")
}
metrics<-lapply(metrics,rbindlist)
cs<-metrics$cs;rl<-metrics$roc;cal<-metrics$calibration;causal<-metrics$causal
diagnostics<-rbindlist(diag)
saveRDS(list(metrics=metrics,diagnostics=diagnostics),cache_path,compress=FALSE)
}
cs<-metrics$cs;rl<-metrics$roc;cal<-metrics$calibration;causal<-metrics$causal
stopifnot(nrow(cs)==18000L,all(cs[, .N,by=.(assv_each,Method)]$N==600L),
  uniqueN(diagnostics$source_sha256)==1L,all(diagnostics$max_alpha_sum_error<1e-5))
methods<-c("SuSiE","SuSiEx intersect","SuSiEx union","MultiSuSiE intersect",
  "MultiSuSiE union","MESuSiE intersect","MESuSiE union","INFUSE intersect","INFUSE union","INFUSE-H")
# The legacy script adds each locus's ROC matrix sequentially and then divides
# by 100. Preserve that floating-point summation order: near-equidistant FDR
# points otherwise select a different cutoff under the legacy exact tie rule.
legacy_mean<-function(x)Reduce(`+`,as.list(x))/length(x)
roc<-rl[,.(Power=legacy_mean(Power),FDR=legacy_mean(FDR),pooled_FDR=sum(FP)/max(1,sum(TP+FP)),
            TP=sum(TP),FP=sum(FP),n_causal=sum(n_causal),n_loci=.N),
         by=.(assv_each,Method,h2_index,base_index,n_base,Cutoff)]
stopifnot(all(roc$n_loci==100L))
setorder(roc,assv_each,Method,h2_index,base_index,-Cutoff)
fdr<-rbindlist(lapply(c(.01,.05,.1,.5),function(target) {
  z<-roc[,{
    jj<-which(abs(FDR-target)==min(abs(FDR-target)))
    .SD[jj[which.max(Power[jj])]]
  },by=.(assv_each,Method,h2_index,base_index,n_base)]
  z[,`:=`(Target_FDR=target,Achieved_FDR=FDR,Power_unrounded=Power,Power=round(Power,2))]
  z
}))
cs_summary<-cs[,.(n_cells=.N,median_CS_union_size=median(Size),mean_CS_union_size=mean(Size),
    CS_power=mean(Power),CS_FDR=mean(FDR),mean_n_cs=mean(n_cs)),by=.(assv_each,Method)]
fd5<-fdr[Target_FDR==.05,.(either_power_5pct=mean(Power),
  either_power_5pct_unrounded=mean(Power_unrounded),mean_achieved_FDR_5pct=mean(Achieved_FDR),
  mean_pooled_FDR_at_selected_cutoff=mean(pooled_FDR)),by=.(assv_each,Method)]
summary<-merge(cs_summary,fd5,by=c("assv_each","Method"))
cal_summary<-cal[,.(n=sum(n),sum_pip=sum(sum_pip),n_causal=sum(n_causal)),
                 by=.(assv_each,Method,h2_index,bin)]
cal_summary[,`:=`(mean_pip=sum_pip/n,observed=n_causal/n)]
cal_summary[,se:=sqrt(observed*(1-observed)/n)]
cal_summary[,`:=`(CI_low=pmax(0,observed-1.96*se),CI_high=pmin(1,observed+1.96*se))]
origin<-causal[,.(n=.N,mean_PIP=mean(PIP),median_PIP=median(PIP),
  recovered_PIP05=sum(PIP>.5),power_PIP05=mean(PIP>.5),CS_power=mean(in_cs)),
  by=.(assv_each,Method,causal_origin)]
origin_threshold<-causal[,rbindlist(lapply(c(.1,.5,.9),function(t)
  data.table(PIP_threshold=t,n=.N,n_recovered=sum(PIP>t),power=mean(PIP>t)))),
  by=.(assv_each,Method,causal_origin)]
fixed_pip<-rl[abs(Cutoff-.5)<1e-10,.(Power=mean(Power),FDR=mean(FDR),
  pooled_FDR=sum(FP)/max(1,sum(TP+FP))),by=.(assv_each,Method)]

# Check comparator metrics against the saved summaries.
validation<-list()
for(a in 0:2) {
  ref<-new.env(parent=emptyenv())
  load(file.path(REVISION,"res","INFUSE_variant_availability",paste0("native_ASSV_each_",a),
    paste0("variant_availability_native_ASSV_each_",a,"_spca.RData")),envir=ref)
  rcs<-as.data.table(ref$all_Set_data_dataframe)
  rcs[,Method:=trimws(as.character(Method))]
  rcs[Method=="MultiSuSiE",Method:="MultiSuSiE intersect"]
  refcs<-rcs[,.(reference_median=median(Size),reference_CS_power=mean(Power)),by=Method]
  rf<-as.data.table(ref$FDR_Power_either)
  rf[,Method:=trimws(as.character(Method))]
  rf[Method=="MultiSuSiE",Method:="MultiSuSiE intersect"]
  rf<-rf[abs(Target_FDR-.05)<1e-10,.(reference_either_power=mean(Power)),by=Method]
  chk<-merge(merge(summary[assv_each==a],refcs,by="Method"),rf,by="Method")
  chk[,`:=`(CS_size_error=abs(median_CS_union_size-reference_median),
    CS_power_error=abs(CS_power-reference_CS_power),
    either_power_error=abs(either_power_5pct-reference_either_power))]
  validation[[a+1L]]<-chk[,.(assv_each,Method,CS_size_error,CS_power_error,either_power_error)]
  rm(ref,rcs,rf,chk);invisible(gc())
}
validation<-rbindlist(validation)
fwrite(validation,file.path(OUT,"tables","frozen_comparator_validation.tsv"),sep="\t")
stopifnot(nrow(validation)==27L,all(validation$CS_size_error<1e-10),
          all(validation$CS_power_error<1e-10),all(validation$either_power_error<1e-10))
tables<-list(summary=summary,CS_per_cell=cs,ROC_per_setting=roc,
  nominal_FDR_per_setting=fdr,calibration=cal_summary,causal_SNP_PIPs=causal,
  causal_origin_recovery=origin,causal_origin_threshold_recovery=origin_threshold,
  PIP05_power_FDR=fixed_pip,fit_diagnostics=diagnostics)
for(k in names(tables))fwrite(tables[[k]],file.path(OUT,"tables",paste0(k,".tsv")),sep="\t")
saveRDS(list(tables=tables,validation=validation,definitions=list(
  CS_size="union of retained component CSs per locus",
  FDR="mean per-locus FDP with zero when no discoveries; nearest-target selection matches original",
  calibration="(0,.1],...,(.9,1] bins; exact zeros excluded, matching original",
  native_recovery="overall PIP truth-stratum recovery; not configuration posterior")),
  file.path(OUT,"summary.rds"),compress=FALSE)

cols<-c("SuSiE"="#777777","SuSiEx intersect"="#80B1D3","SuSiEx union"="#0072B2",
  "MultiSuSiE intersect"="#E69F00","MultiSuSiE union"="#D55E00",
  "MESuSiE intersect"="#A6D8C9","MESuSiE union"="#009E73",
  "INFUSE intersect"="#C8B3D8","INFUSE union"="#CC79A7","INFUSE-H"="#222222")
labels<-setNames(gsub(" ","\n",availability_method_labels(methods)),methods)
for(k in names(tables))if("Method"%in%names(tables[[k]]))set(tables[[k]],j="Method",value=factor(tables[[k]]$Method,levels=methods))
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
  axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",
  plot.tag=element_text(face="bold"),strip.background=element_rect(fill="grey95")))
fill<-function()scale_fill_manual(values=cols,drop=FALSE, labels = availability_method_labels)
colour<-function()scale_colour_manual(values=cols,drop=FALSE, labels = availability_method_labels)
xx<-function()scale_x_discrete(labels=labels,drop=FALSE)
export<-function(plot,stem,width=14,height=9) {
  ggsave(paste0(stem,".pdf"),plot=plot,width=width,height=height,device=cairo_pdf,limitsize=FALSE)
  ggsave(paste0(stem,".png"),plot=plot,width=width,height=height,dpi=300,limitsize=FALSE,bg="white")
}
for(a in 0:2) {
  figdir<-file.path(OUT,"figures",paste0("native_ASSV_each_",a));dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  dcs<-tables$CS_per_cell[assv_each==a]
  ds<-tables$summary[assv_each==a]
  dr<-tables$ROC_per_setting[assv_each==a]
  df<-tables$nominal_FDR_per_setting[assv_each==a]
  cr<-dr[,.(Power=mean(Power),FDR=mean(FDR)),by=.(Method,Cutoff)]
  psize<-ggplot(dcs,aes(Method,log2(Size+1),fill=Method))+geom_boxplot(outlier.size=.3)+
    fill()+xx()+labs(x=NULL,y="log2(Credible-set union size + 1)")+theme(legend.position="none")
  ppow<-ggplot(ds,aes(Method,CS_power,fill=Method))+geom_col()+
    geom_text(aes(label=sprintf("%.3f",CS_power)),vjust=-.3,size=3)+
    fill()+xx()+coord_cartesian(ylim=c(0,1))+labs(x=NULL,y="Credible-set power")+theme(legend.position="none")
  proc<-ggplot(cr,aes(Power,1-FDR,colour=Method,group=Method))+
    geom_path(linewidth=.8)+colour()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(x="Detection power",y="1 - FDR")+theme(axis.text.x=element_text(angle=0),legend.text=element_text(size=8))
  pfdr<-ggplot(ds,aes(Method,either_power_5pct,fill=Method))+geom_col()+
    geom_text(aes(label=sprintf("%.2f",either_power_5pct)),vjust=-.3,size=3)+fill()+xx()+
    coord_cartesian(ylim=c(0,1))+labs(x=NULL,y="Detection power at nominal 5% FDR")+theme(legend.position="none")
  overall<-(psize+ppow)/(proc+pfdr)+plot_annotation(
    title=paste0("Additional ancestry-unique causal variants per ancestry: ",a),
    caption="600 simulation evaluations per method.",tag_levels="a")
  export(overall,file.path(figdir,"summarize"),16,11)
  export(psize+ppow+plot_annotation(tag_levels="a"),file.path(figdir,"Set_size_power"),16,6)
  dr[,h2:=factor(h2_index,labels=c("h2 = 1e-4","h2 = 2e-4"))]
  pr<-ggplot(dr,aes(Power,1-FDR,colour=Method))+geom_path(linewidth=.7)+
    facet_grid(h2~n_base,labeller=label_both)+colour()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(x="Detection power",y="1 - FDR")+theme(axis.text.x=element_text(angle=0))
  fd<-df[Target_FDR%in%c(.01,.05,.1),.(Power=mean(Power)),by=.(Method,Target_FDR)]
  pb<-ggplot(fd,aes(Method,Power,fill=Method))+geom_col()+facet_wrap(~Target_FDR,labeller=label_both)+
    fill()+xx()+coord_cartesian(ylim=c(0,1))+labs(x=NULL,y="Detection power at nominal FDR")+theme(legend.position="none")
  export(pr/pb+plot_layout(heights=c(2,1.2))+plot_annotation(tag_levels="a"),
    file.path(figdir,"PIP_ROC_FDR_Power_Either"),17,12)
  dc<-tables$calibration[assv_each==a]
  dc[,h2:=factor(h2_index,labels=c("h2 = 1e-4","h2 = 2e-4"))]
  pc<-ggplot(dc,aes(mean_pip,observed,colour=Method))+geom_abline(slope=1,intercept=0,colour="grey60",linetype=2)+
    geom_errorbar(aes(ymin=CI_low,ymax=CI_high),width=.02,linewidth=.3)+geom_line()+geom_point(size=1)+
    facet_grid(h2~Method)+colour()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(x="Mean predicted PIP",y="Observed causal proportion",caption="PIP=0 excluded, matching original calibration bins.")+
    theme(legend.position="none",axis.text.x=element_text(angle=0),strip.text.x=element_text(size=7))
  export(pc,file.path(figdir,"pip_either_calibration"),20,6)
  if(a>0L) {
    dor<-tables$causal_origin_recovery[assv_each==a & grepl("^native_ASSV",causal_origin)]
    po<-ggplot(dor,aes(Method,power_PIP05,fill=Method))+geom_col()+
      facet_wrap(~causal_origin)+fill()+xx()+coord_cartesian(ylim=c(0,1))+
      labs(x=NULL,y="Recovery at overall PIP > 0.5",caption="Truth-stratum recovery, not ancestry-configuration classification.")+theme(legend.position="none")
    export(po,file.path(figdir,"native_unique_recovery"),16,6)
  }
  dd<-tables$fit_diagnostics[assv_each==a]
  pd<-ggplot(dd,aes(factor(n_union),runtime_minutes,colour=converged))+
    geom_point(size=.6)+scale_colour_manual(values=c(`TRUE`="#0072B2",`FALSE`="#D55E00"))+
    labs(x="Locus SNP count",y="Runtime (minutes)",colour="Converged")+
    scale_x_discrete(breaks=function(x)x[seq(1,length(x),length.out=min(10,length(x)))])
  export(pd,file.path(figdir,"fit_diagnostics"),12,6)
}
allsummary<-tables$summary
combined<-ggplot(allsummary,aes(Method,either_power_5pct,fill=Method))+geom_col()+
  facet_wrap(~assv_each,labeller=label_both)+fill()+xx()+coord_cartesian(ylim=c(0,1))+
  labs(x=NULL,y="Detection power at nominal 5% FDR")+theme(legend.position="none")
export(combined,file.path(OUT,"figures","comparison_across_scenarios"),18,6)
print(summary[Method%in%c("INFUSE union","INFUSE-H","MESuSiE union","SuSiEx union","MultiSuSiE intersect","MultiSuSiE union")])
