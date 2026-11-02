#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
stopifnot(nzchar(Sys.getenv("SLURM_JOB_ID")))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
setDTthreads(1L)
revision <- .infuse_work
module <- file.path(revision,"LD_mismatch_sensitivity")
outdir <- file.path(module,"results/ld_attenuation/summary")
figdir <- file.path(revision,"figures/LD_mismatch_sensitivity/ld_attenuation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
files <- file.path(module,"results/ld_attenuation",sprintf("locus_%03d.rds",1:100))
stopifnot(all(file.exists(files)))
objects <- lapply(files,readRDS)
d <- rbindlist(lapply(objects,`[[`,"metrics"))
keys <- c("locus","causal_index","h2_index","family","level","method","scope")
stopifnot(nrow(d)==108000L,!anyDuplicated(d[,..keys]),
  all(d$tp+d$fp==d$n_flag),all(d$tp+d$fn==d$n_target))
coverage <- d[,.(n_evaluations=.N,n_loci=uniqueN(locus)),by=.(family,level,method,scope)]
stopifnot(all(coverage$n_evaluations==600L),all(coverage$n_loci==100L))
fwrite(coverage,file.path(outdir,"coverage.tsv"),sep="\t")
fwrite(rbindlist(lapply(objects,`[[`,"parameters")),file.path(outdir,"parameters.tsv"),sep="\t")
fwrite(data.table(locus=1:100,seconds=vapply(objects,function(x)as.numeric(x$seconds),0)),file.path(outdir,"runtime.tsv"),sep="\t")
source(file.path(.infuse_repo, "Simulation/revision/ld_mismatch/summary_metrics.R"))
set.seed(20260917L)
boot <- replicate(2000L,sample.int(100L,100L,replace=TRUE))
saveRDS(boot,file.path(outdir,"bootstrap_locus_indices.rds"))
stopifnot(identical(dim(boot),c(100L,2000L)))
add_evaluation_metrics(d)
groups <- c("family","level","method","scope")
overall <- d[,cluster_summary(.SD),by=groups]
arch <- d[,cluster_summary(.SD),by=c(groups,"causal_index","h2_index")]
validate_summaries(d,overall,groups)
validate_summaries(d,arch,c(groups,"causal_index","h2_index"))
fwrite(d,file.path(outdir,"evaluation_metrics.tsv"),sep="\t")
fwrite(overall,file.path(outdir,"overall_estimates.tsv"),sep="\t")
fwrite(arch,file.path(outdir,"architecture_estimates.tsv"),sep="\t")
# Paired differences for all three method pairs, resampling common loci.
pairs <- combn(c("joint_shared","joint_separate","separate_union"),2,simplify=FALSE)
contrasts <- d[scope=="Either",{
  cell <- .SD
  rbindlist(lapply(pairs,function(pair) {
    aa <- lapply(pair,function(mid) {
      x <- cell[method==mid,.(ss=sum(sensitivity_eval),ns=.N,sf=sum(fdp),nf=sum(n_flag>0)),by=locus]
      setorder(x,locus); stopifnot(identical(x$locus,1:100));x
    })
    rbindlist(lapply(c("sensitivity","conditional_mean_fdp"),function(metric) {
      num <- if(metric=="sensitivity")"ss" else "sf"
      den <- if(metric=="sensitivity")"ns" else "nf"
      point <- vapply(aa,function(x)sum(x[[num]])/sum(x[[den]]),0)
      bs <- lapply(aa,function(x)colSums(matrix(x[[num]][boot],100))/colSums(matrix(x[[den]][boot],100)))
      dif <- bs[[1]]-bs[[2]]; dif <- dif[is.finite(dif)]
      ci <- if(length(dif))quantile(dif,c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)
      data.table(method_1=pair[1],method_2=pair[2],metric=metric,estimate=point[1]-point[2],lower=ci[1],upper=ci[2],n_boot_valid=length(dif))
    }))
  }))
},by=.(family,level)]
fwrite(contrasts,file.path(outdir,"paired_method_differences.tsv"),sep="\t")
overview <- dcast(overall[scope=="Either" & metric %in% c("sensitivity","conditional_mean_fdp","fp_per_100k")],
  family+level+method~metric,value.var="estimate")
fwrite(overview,file.path(outdir,"attenuation_overall_comparison.tsv"),sep="\t")
source(file.path(.infuse_repo, "Simulation/revision/common/plot_theme.R"))
spec <- infuse_plot_spec()
labels <- c(joint_shared="Joint test (shared s)",joint_separate="Joint test (ancestry-specific s)",separate_union="Separate tests (union)")
colors <- c(joint_shared="#0072B2",joint_separate="#009E73",separate_union="#D55E00")
shapes <- c(joint_shared=16,joint_separate=15,separate_union=17)
lines <- c(joint_shared="solid",joint_separate="dashed",separate_union="dashed")
wt <- theme_infuse_plot(spec)+theme(plot.title=element_text(size=11,face="plain",hjust=0),
  legend.title=element_blank(),legend.text=element_text(size=7),legend.key.width=unit(4,"mm"),legend.spacing.x=unit(1.5,"mm"))
plot_panel <- function(data,metric_id,ylabel,tag,architecture=FALSE,sparse=FALSE) {
  x <- copy(data[metric==metric_id]); x[,Method:=factor(method,levels=names(labels))]
  if(sparse) x[,X:=factor(scenario,levels=c("sparse_negative_100","sparse_negative_050","sparse_attenuation_050","sparse_attenuation_025"),
    labels=c("Full\nreversal","Partial\nreversal","Moderate\natten.","Mild\natten."))] else x[,X:=factor(level,levels=sort(unique(level)))]
  if(architecture) {
    x[,Causal:=factor(causal_index,levels=1:3,labels=c("2 causal SNPs","3 causal SNPs","5 causal SNPs"))]
    x[,Heritability:=factor(h2_index,levels=1:2,labels=c("h² = 10⁻⁴","h² = 2 × 10⁻⁴"))]
  } else x[,Scenario:=factor(family,levels=c("shrink_WB","shrink_BB"),labels=c("Shrink WB LD only","Shrink BB LD only"))]
  pd <- position_dodge(.23)
  p <- ggplot(x,aes(X,estimate,color=Method,shape=Method,linetype=Method,group=Method))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.10,linewidth=.25,position=pd)+
    geom_line(linewidth=.35,position=pd)+geom_point(size=1.35,position=pd)+
    scale_color_manual(values=colors,labels=labels)+scale_shape_manual(values=shapes,labels=labels)+
    scale_linetype_manual(values=lines,labels=labels)+labs(x=if(sparse)"LD perturbation" else "LD shrinkage λ",y=ylabel,title=tag)+wt
  if(architecture)p <- p+facet_grid(Heritability~Causal) else p <- p+facet_grid(~Scenario)
  if(metric_id=="sensitivity")p <- p+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25)) else
    p <- p+scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.08)))
  built <- ggplot_build(p)$data
  stopifnot(nrow(built[[3]])==nrow(x),isTRUE(all.equal(sort(built[[3]]$y),sort(x$estimate))),
    isTRUE(all.equal(sort(built[[1]]$ymin),sort(x$lower))),isTRUE(all.equal(sort(built[[1]]$ymax),sort(x$upper))))
  p
}
save_figure <- function(data,stem,architecture=TRUE,sparse=FALSE) {
  a <- plot_panel(data,"sensitivity",if(architecture)"Sensitivity" else "Overall sensitivity","a",architecture,sparse)
  b <- plot_panel(data,if(sparse)"fp_per_100k" else "conditional_mean_fdp",
    if(sparse)"False-positive calls\nper 100,000 SNPs" else "False discovery rate","b",architecture,sparse)
  p <- (a/b)+plot_layout(guides="collect") & theme(legend.position="bottom")
  for(ext in c("png","pdf"))ggsave(file.path(figdir,paste0(stem,".",ext)),p,width=170,
    height=if(architecture)155 else 130,units="mm",dpi=600,bg="white",device=if(ext=="pdf")cairo_pdf else "png")
  fwrite(data,file.path(outdir,paste0(stem,"_plot_data.tsv")),sep="\t")
}
plot_metrics <- c("sensitivity","conditional_mean_fdp")
save_figure(arch[family=="shrink_both" & scope=="Either" & metric %in% plot_metrics],"Bilateral_LD_attenuation_three_methods")
save_figure(overall[family %in% c("shrink_WB","shrink_BB") & level %in% c(0,.25,.5,.75) & scope=="Either" & metric %in% plot_metrics],
  "Asymmetric_LD_attenuation_three_methods_overall",FALSE)
for(fam in c("shrink_WB","shrink_BB"))save_figure(
  arch[family==fam & level %in% c(0,.25,.5,.75) & scope=="Either" & metric %in% plot_metrics],
  paste0("Asymmetric_LD_attenuation_",sub("shrink_","",fam),"_three_methods_by_architecture"))
print(overview[family=="shrink_both" | level %in% c(0,.25,.5,.75)])
cat("Complete benchmark written to",outdir,"\nFigures:",figdir,"\n")
