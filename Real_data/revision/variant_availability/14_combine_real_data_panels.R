#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
x<-readRDS(file.path(root,"figures/cs_summary/cs_panel_objects.rds"))
ee<-readRDS(file.path(root,"summary/matched_size/enrichment_plot_objects.rds"))
out<-file.path(root,"figures/real_data_summary");dir.create(out,recursive=TRUE,showWarnings=FALSE)
style<-theme_bw(base_size=10,base_family="sans")+theme(
 axis.text=element_text(size=8,color="black"),axis.title=element_text(size=10,face="bold"),
 strip.text=element_text(size=10,face="bold"),strip.background=element_blank(),
 panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(color="black",linewidth=.3),
 plot.title=element_text(size=11,hjust=.5),plot.tag=element_text(size=12,face="bold"),
 legend.position="bottom",legend.title=element_text(size=9,face="bold"),legend.text=element_text(size=9))
# Restyle saved plots only; avoid rereading all genotype/LD inputs for typography.
for(nm in names(x$plots)){
 x$plots[[nm]]<-x$plots[[nm]]+style
 for(ext in c("pdf","png"))ggsave(file.path(root,"figures/cs_summary",paste0("panel_",nm,".",ext)),x$plots[[nm]],
   width=if(nm=="a")6.2 else 5.5,height=3.8,units="in",dpi=500,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png")
}
abc<-wrap_plots(x$plots,nrow=1,guides="collect")+plot_annotation(caption=NULL)
for(ext in c("pdf","png"))ggsave(file.path(root,"figures/cs_summary",paste0("cs_panels_abc.",ext)),abc,
 width=14,height=4.4,units="in",dpi=500,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png")
for(g in names(ee)){
 a<-x$plots$a+labs(title=NULL,x="Regions with credible sets")+
   guides(fill="none",pattern="none")+
   theme(legend.position="none",axis.title.x=element_text(size=8),axis.text.y=element_text(size=7))
 b<-x$plots$b+labs(title=NULL,fill="Method",y=expression(log[2]("Set size + 1")))+
   guides(fill="none")+theme(legend.position="none")
 c<-x$plots$c+labs(title=NULL,
   y="Median maximum |Z|")+guides(fill="none")+theme(legend.position="none",
   plot.title=element_text(size=10,hjust=.5))
 # Narrow reference-style left column: keep adjacent numerical labels distinct.
 c$layers[[2]]<-geom_text(aes(label=ifelse(is.finite(median_max_abs_observed_Z),
   sprintf("%.1f",median_max_abs_observed_Z),"")),
   position=position_dodge(width=.8),vjust=-.35,size=1.8,na.rm=TRUE)
 d<-ee[[g]]+labs(tag="d",title=NULL,subtitle=NULL,caption=NULL)+style+
   theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text=element_text(size=9))
 # Reference 9_plot_enrichment_either_bmi/bp: stacked CS/Z on the left,
 # annotation facets on the right. Retain the additionally requested region-count
 # panel at the top of the left column.
 left<-wrap_plots(a,b,c,ncol=1,heights=c(1,1,1),guides="keep")
 left<-wrap_elements(full=patchworkGrob(left))
 trait_label<-if(g=="BMI") "BMI" else "blood pressure"
 p<-wrap_plots(left,d,nrow=1,widths=c(1.4,3))+plot_annotation(
   title=paste0("All of Us: credible-set summaries and ",trait_label," enrichment"),
   caption=NULL,theme=theme(plot.title=element_text(size=12,hjust=.5),plot.margin=margin(8,8,8,8)))
 for(ext in c("pdf","png"))ggsave(file.path(out,paste0("real_data_summary_",g,".",ext)),p,
   width=12,height=7.6,units="in",dpi=500,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png")
}
# Latest requested final composition: trait-specific BMI and BP blocks on one page.
source(file.path(.infuse_repo, "Real_data/revision/variant_availability/15_final_bmi_bp_figure.R"),local=new.env())
