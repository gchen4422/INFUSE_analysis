#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
input<-file.path(root,"summary/matched_size/functional_enrichment.tsv");d<-fread(input)
# Plot the selected methods and annotation categories; retain full summaries separately.
d<-d[method!="SuSiE_posthoc" & annotation!="missense"]
d<-d[annotation!="CRE"]
fwrite(d,file.path(root,"summary/enrichment_plot_data.tsv"),sep="\t",na="NA")
features<-c("non_synonymous_common","synonymous","promotor","utr_comb","Conserved_LindbladToh_common","H3K4me3_Trynka_common","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")
fl<-c("Non-synonymous","Synonymous","Promoter","UTR","Conserved","H3K4me3","CRE","Heart eQTL","Artery eQTL","Brain eQTL")
methods<-c("SuSiEx_union","MESuSiE_union","INFUSE_union")
labs<-c("SuSiEx union","INFUSE\u2212 union","INFUSE union")
pal<-setNames(c("#CC79A7","#009E73","#0072B2"),labs)
# Either panels follow the original Figure 4 pastel palette. Preserve the
# independently requested shared/ancestry-specific palette above.
either_pal<-setNames(c("#CC79A7","#84c3b7","#71b8ed"),labs)
d[,annotation:=factor(annotation,levels=features,labels=fl)];d[,method:=factor(method,levels=methods,labels=labs)]
stopifnot(!anyNA(d$annotation),!anyNA(d$method))
out<-file.path(root,"figures/enrichment");dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Adapt reference BMI/BP/combined custom_theme: black axes, no grid/border,
# bold strip/axis titles, method colors, and bottom method/configuration legends.
reference_theme<-theme_bw(base_size=9,base_family="sans")+theme(
 legend.position="bottom",legend.title=element_text(size=9,face="bold"),legend.text=element_text(size=9),
 strip.background=element_blank(),strip.text=element_text(size=9,face="bold"),
 axis.text=element_text(size=7.5,color="black"),axis.title=element_text(size=9,face="bold"),
 plot.title=element_text(size=11,hjust=.5),plot.subtitle=element_text(size=8,hjust=.5),
 plot.caption=element_text(size=7,hjust=0),plot.tag=element_text(size=12,face="bold"),
 panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(color="black",linewidth=.3),
 panel.spacing=grid::unit(.8,"lines"),plot.margin=margin(6,8,6,8))
theme_set(reference_theme)
saved<-list();warn<-character();plots<-list()
export<-function(p,name,width=8.8,height=5.5){for(ext in c("pdf","png")){
 f<-file.path(out,paste0(name,".",ext))
 withCallingHandlers(ggsave(f,p,width=width,height=height,units="in",dpi=500,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png"),
   warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")})
 saved[[length(saved)+1L]]<<-data.table(file=f,source_table=input)
}}
for(g in c("BMI","Blood_pressure")){
 x<-d[group==g & target=="Either"];ttl<-if(g=="BMI")"BMI" else "Blood pressure (DBP + SBP)"
 x[,annotation:=droplevels(annotation)]
 breaks<-c(.5,.6,.7,.8,.9)
 sizes<-unique(x[,.(INFUSE_threshold,N)])
 axis_labels<-paste0(sprintf("%.2f",breaks),"\n(n=",sizes$N[match(breaks,sizes$INFUSE_threshold)],")")
 # Keep the plot canvas free of calculation notes.
 caption<-NULL
 p<-ggplot(x,aes(INFUSE_threshold,enrichment,color=method,shape=method,linetype=method,group=method))+
 geom_hline(yintercept=1,linetype="dotted",color="grey60",linewidth=.3)+geom_line(linewidth=.55,na.rm=TRUE)+geom_point(size=1.4,na.rm=TRUE)+
 facet_wrap(~annotation,ncol=4,scales="free_y",drop=FALSE)+scale_color_manual(values=either_pal,name="Method")+scale_x_continuous(breaks=breaks,labels=axis_labels)+
 labs(title=paste(ttl,"- Matched-size enrichment"),x="PIP threshold (matched N)",y="Enrichment (selected / background)",caption=NULL,shape="Method",linetype="Method")
 export(p,paste0(g,"_matched_size_either_thresholds"))
 b<-x[abs(INFUSE_threshold-.5)<1e-8];N<-unique(b$N);stopifnot(length(N)==1L)
 p<-ggplot(b,aes(method,enrichment,fill=method))+geom_col(width=.7,na.rm=TRUE,colour=NA)+geom_hline(yintercept=1,linetype="dotted",color="grey60",linewidth=.3)+
 facet_wrap(~annotation,ncol=4,scales="free_y",drop=FALSE)+scale_fill_manual(values=either_pal,name="Method")+scale_y_continuous(expand=expansion(mult=c(0,.1)))+
 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
 labs(title=paste(ttl,"- Matched-size enrichment"),subtitle=NULL,x=NULL,y="Enrichment (selected / background)",caption=NULL)
 export(p,paste0(g,"_matched_size_either_cutoff_0.50"));plots[[g]]<-p
 x<-copy(d[group==g & target %in% c("Shared","Unique")]);x[,configuration:=factor(target,levels=c("Shared","Unique"),labels=c("Shared","Ancestry-specific"))]
 x[,annotation:=droplevels(annotation)]
 cfg_sizes<-unique(x[,.(target,INFUSE_threshold,N)])
 n_shared<-cfg_sizes[target=="Shared"]$N[match(breaks,cfg_sizes[target=="Shared"]$INFUSE_threshold)]
 n_unique<-cfg_sizes[target=="Unique"]$N[match(breaks,cfg_sizes[target=="Unique"]$INFUSE_threshold)]
 cfg_labels<-paste0(sprintf("%.2f",breaks),"\nS=",n_shared,"\nA=",n_unique)
 p<-ggplot(x,aes(INFUSE_threshold,enrichment,color=method,linetype=configuration,group=interaction(method,configuration)))+
 geom_line(linewidth=.55,na.rm=TRUE)+facet_wrap(~annotation,ncol=4,scales="free_y",drop=FALSE)+scale_color_manual(values=pal,drop=TRUE,name="Method")+
 scale_linetype_manual(values=c("Shared"="solid","Ancestry-specific"="dashed"),name="Variant")+
 scale_x_continuous(breaks=breaks,labels=cfg_labels)+
 labs(title=paste(ttl,"- Shared and ancestry-specific enrichment"),x="PIP threshold (matched N)",y="Enrichment (selected / background)",caption=NULL)
 export(p,paste0(g,"_matched_size_shared_ancestry_specific"))
}
# Figure5 pools all three traits before ranking: six annotations, two thresholds.
pooled_path<-file.path(root,"summary/matched_size/pooled_BMI_BP/functional_enrichment.tsv")
pooled_sources<-fread(file.path(dirname(pooled_path),"source_checksums.tsv"))
stopifnot(identical(pooled_sources$md5,unname(tools::md5sum(pooled_sources$file))))
pooled<-fread(pooled_path)
stopifnot(nrow(pooled)==48L,setequal(pooled$group,"BMI_BP_pooled"))
fwrite(pooled,file.path(root,"summary/pooled_configuration_enrichment_plot_data.tsv"),sep="\t",na="NA")
pooled[,annotation:=factor(annotation,levels=features[1:6],labels=fl[1:6])]
pooled[,method:=factor(method,levels=c("MESuSiE_union","INFUSE_union"),labels=c("INFUSE\u2212","INFUSE"))]
pooled[,configuration:=factor(target,levels=c("Shared","Unique"),labels=c("Shared","Ancestry-specific"))]
stopifnot(!anyNA(pooled$annotation),!anyNA(pooled$method))
config_bars<-lapply(c(.5,.9),function(cut){
 b<-pooled[abs(INFUSE_threshold-cut)<1e-8]
 ggplot(b,aes(x=factor(INFUSE_threshold),y=enrichment,fill=method,linetype=configuration,
   group=interaction(method,configuration)))+
  geom_col(position=position_dodge(width=.85),width=.78,color="black",linewidth=.6,na.rm=TRUE)+
  facet_wrap(~annotation,ncol=3,scales="free_y")+
  scale_fill_manual(values=c("INFUSE\u2212"="#009E73","INFUSE"="#0072B2"),name="Method")+
  scale_linetype_manual(values=c("Shared"="solid","Ancestry-specific"="dashed"),name="Variant")+
  scale_y_continuous(expand=expansion(mult=c(0,.08)))+
  labs(title=paste0("Shared and ancestry-specific enrichment at PIP = ",sprintf("%.1f",cut)),
    x=NULL,y="Enrichment (selected / background)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text=element_text(size=11,face="bold"),
    plot.title=element_text(size=11,hjust=.5))
})
combined<-wrap_plots(config_bars,ncol=1,guides="collect")+plot_annotation(tag_levels="a")&theme(legend.position="bottom")
export(combined,"BMI_BP_matched_size_shared_ancestry_specific_cutoffs_0.50_0.90",width=8.8,height=8.8)
for(ext in c("png","pdf")){
 f<-file.path(root,"figures",paste0("Figure_5_shared_ancestry_specific_BMI_BP_union.",ext))
 stopifnot(file.copy(file.path(out,paste0("BMI_BP_matched_size_shared_ancestry_specific_cutoffs_0.50_0.90.",ext)),f,overwrite=TRUE))
 saved[[length(saved)+1L]]<-data.table(file=f,source_table=pooled_path)
}
# Keep superseded single-group bars recoverable, outside the current figure set.
archive<-file.path(root,"archive/configuration_bars_before_pooled_20260910")
for(g in c("BMI","Blood_pressure"))for(ext in c("png","pdf")){
 f<-file.path(out,paste0(g,"_matched_size_shared_ancestry_specific_cutoffs_0.50_0.90.",ext))
 if(file.exists(f)){
  dir.create(archive,recursive=TRUE,showWarnings=FALSE)
  dest<-file.path(archive,basename(f))
  if(file.exists(dest))stop("Archive target exists; refusing overwrite: ",dest)
  stopifnot(file.rename(f,dest))
 }
}
for(i in seq_along(saved))if(grepl("BMI_BP_matched_size_shared",saved[[i]]$file))saved[[i]]$source_table<-pooled_path
fwrite(rbindlist(saved),file.path(out,"figure_manifest.tsv"),sep="\t")
saveRDS(plots,file.path(root,"summary/matched_size/enrichment_plot_objects.rds"))
cat("Matched-size plots complete:",length(saved),"files;",length(warn),"captured warnings\n")
