#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Combine credible-set and enrichment panels for BMI and blood pressure.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
out<-file.path(root,"figures/real_data_summary")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sources<-file.path(root,c("summary/cs_panels/panel_b_nonempty_locus_cs_sizes.tsv",
 "summary/cs_panels/panel_c_median_max_abs_Z.tsv",
 "summary/matched_size/functional_enrichment.tsv",
 "summary/matched_size/enrichment_plot_objects.rds"))
before<-tools::md5sum(sources)
cs<-fread(sources[1]);zz<-fread(sources[2]);ee<-readRDS(sources[4])
methods<-c("SuSiEx_union","MESuSiE_union","INFUSE_union","SuSiE_posthoc")
short<-c("SuSiEx","INFUSE\u2212","INFUSE","SuSiE")
pal<-setNames(c("#D95F3B","#84c3b7","#71b8ed","#E69F00"),short)
stopifnot(nrow(cs)==1058L,nrow(zz)==12L,!anyDuplicated(zz[,.(trait,method)]),
 all(cs$distinct_cs_snps>0),max(abs(cs$log2_distinct_cs_snps_plus1-log2(cs$distinct_cs_snps+1)))<1e-12)
style<-theme_classic(base_size=9,base_family="sans")+theme(
 axis.text=element_text(size=7,color="black"),axis.title=element_text(size=9,face="bold"),
 strip.text=element_text(size=9,face="bold"),strip.background=element_blank(),
 panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(color="black",linewidth=.3),
 plot.title=element_text(size=12,hjust=.5),plot.tag=element_text(size=11),
 legend.position="none",plot.margin=margin(5,7,5,7),panel.spacing=grid::unit(.7,"lines"))
theme_set(style)
exports<-list();display_cs<-list();display_z<-list();rows<-list()
save_plot<-function(p,path,height){
 for(ext in c("pdf","png")){
  f<-paste0(path,".",ext)
  ggsave(f,p,width=12,height=height,units="in",dpi=500,bg="white",
    device=if(ext=="pdf")grDevices::cairo_pdf else "png")
  exports[[length(exports)+1L]]<<-data.table(file=f,width_inches=12,height_inches=height)
 }
}
for(g in c("BMI","Blood_pressure")){
 tr<-if(g=="BMI") "BMI" else c("DBP","SBP")
 b<-copy(cs[trait %in% tr]);z<-copy(zz[trait %in% tr])
 stopifnot(setequal(unique(b$trait),tr),setequal(unique(z$trait),tr),
   nrow(z)==length(tr)*4L,setequal(z$method,methods))
 # Omit post-hoc SuSiE from the BMI and blood-pressure credible-set and Z-score panels.
 # SuSiEx, all underlying fits, and the source summaries are unchanged.
 b<-b[method!="SuSiE_posthoc"];z<-z[method!="SuSiE_posthoc"]
 stopifnot(nrow(z)==length(tr)*3L,setequal(z$method,methods[1:3]),
   !any(b$method=="SuSiE_posthoc"))
 # BP CS/Z show DBP and SBP separately; only BP enrichment pools before ranking.
 display_cs[[g]]<-copy(b);display_z[[g]]<-copy(z)
 b[,Method:=factor(method,levels=methods,labels=short)]
 z[,Method:=factor(method,levels=methods,labels=short)]
 b[,trait:=factor(trait,levels=tr)];z[,trait:=factor(trait,levels=tr)]
 # Stagger adjacent BP labels slightly; retain their unmodified numerical values.
 z[,label_y:=median_max_abs_observed_Z+if(g=="BMI") 0 else .22*(as.integer(Method) %% 2L==0L)]
 tags<-if(g=="BMI") c("a","b","c") else c("d","e","f")
 x_cs<-if(g=="BMI") aes(x=Method,y=log2_distinct_cs_snps_plus1,fill=Method) else
   aes(x=trait,y=log2_distinct_cs_snps_plus1,fill=Method)
 x_z<-if(g=="BMI") aes(x=Method,y=median_max_abs_observed_Z,fill=Method) else
   aes(x=trait,y=median_max_abs_observed_Z,fill=Method)
 p_cs<-ggplot(b,x_cs)+geom_boxplot(width=.75,position=position_dodge2(width=.85,preserve="single"),
   outlier.size=.45,outlier.alpha=.5,linewidth=.3,color="darkgray")+
   scale_fill_manual(values=pal,drop=FALSE)+
   labs(tag=tags[1],x=if(g=="BMI") "Method" else "Trait",y=expression(bold(log[2]("Set size + 1"))))
 p_z<-ggplot(z,x_z)+geom_col(width=.75,position=position_dodge(width=.8))+
   geom_text(aes(y=label_y,label=sprintf("%.2f",median_max_abs_observed_Z)),
     position=position_dodge(width=.8),vjust=-.4,size=2.1)+
   scale_fill_manual(values=pal,drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(0,.12)))+
   labs(tag=tags[2],x=if(g=="BMI") "Method" else "Trait",y="Median maximum |Z|")
 # Reference annotations with trait-relevant eQTL panels; preserve stacked layout.
 omit<-if(g=="BMI") c("Heart eQTL","Artery eQTL") else "Brain eQTL"
 ee[[g]]$data<-copy(ee[[g]]$data)[!as.character(annotation) %in% omit]
 ee[[g]]$data[,annotation:=droplevels(annotation)]
 stopifnot(uniqueN(ee[[g]]$data$annotation)==if(g=="BMI") 7L else 8L)
 e<-ee[[g]]+style+labs(tag=tags[3],title=NULL,subtitle=NULL,caption=NULL,x="Method")+
   scale_x_discrete(labels=c("INFUSE union"="INFUSE","INFUSE\u2212 union"="INFUSE\u2212","SuSiEx union"="SuSiEx"))+
   scale_fill_manual(values=setNames(unname(pal[1:3]),c("SuSiEx union","INFUSE\u2212 union","INFUSE union")))+
   guides(fill="none")+theme(axis.text.x=element_text(size=7,color="black"),
     axis.ticks.x=element_line(linewidth=.3),legend.position="none")
 left<-wrap_plots(p_cs,p_z,ncol=1,heights=c(1,1))
 left<-wrap_elements(full=patchworkGrob(left))
 p<-wrap_plots(left,e,nrow=1,widths=c(1,3))+plot_annotation(
   title=paste("Multi-ancestry fine-mapping on All of Us",if(g=="BMI") "BMI" else "BP","data"),
   theme=theme(plot.title=element_text(size=12,hjust=.5),plot.margin=margin(5,5,5,5)))
 save_plot(p,file.path(out,paste0("real_data_summary_",g)),5)
 rows[[g]]<-wrap_elements(full=patchworkGrob(p))
}
combined<-wrap_plots(rows,ncol=1,heights=c(1,1))
save_plot(combined,file.path(root,"figures/Figure_4_set_size_enrichment_BMI_BP_union"),10)
stopifnot(identical(before,tools::md5sum(sources)))
fwrite(rbindlist(exports),file.path(out,"figure_manifest.tsv"),sep="\t")
# Exact plotting records, not new estimates; source numerical tables are read-only.
fwrite(rbindlist(display_cs,idcol="figure_group"),file.path(out,"displayed_cs_size_records.tsv"),sep="\t")
fwrite(rbindlist(display_z,idcol="figure_group"),file.path(out,"displayed_z_summary.tsv"),sep="\t")
fwrite(data.table(file=names(before),md5=unname(before)),file.path(out,"source_checksums.tsv"),sep="\t")
