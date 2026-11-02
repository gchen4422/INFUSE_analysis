#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Panels a-c corresponding to the source real-data summary, with explicit estimands.
# Execute on Slurm: observed input matrices are loaded sequentially, never together.
# Outputs are confined to this task; no CI/p-values or causal-accuracy claims.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1)
root <- normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
args <- commandArgs(TRUE)
stopifnot(length(args)==0L || identical(args,"--replot"))
replot <- identical(args,"--replot")
methods <- c("SuSiEx_union","MESuSiE_union","INFUSE_union","SuSiE_posthoc")
method_labels <- c("SuSiEx union","INFUSE\u2212 union","INFUSE union","SuSiE post hoc")
method_colors <- setNames(c("#CC79A7","#009E73","#0072B2","#E69F00"),method_labels)
traits <- c("BMI","DBP","SBP")
region_file <- file.path(root,"summary","region_metrics.tsv")
pip_file <- file.path(root,"summary","variant_pips.tsv.gz")
if(!replot) {
region <- fread(region_file)
pips <- fread(pip_file,select=c("trait","locus","SNP","method","in_cs","analyzed",
                               "available_EUR","available_AFR"))
stopifnot(nrow(region)==289L*length(methods),!anyDuplicated(region[,.(trait,locus,method)]),
  setequal(region$method,methods),setequal(region$trait,traits),
  !anyDuplicated(pips[,.(trait,locus,SNP,method)]),!anyNA(pips$in_cs),
  !any(pips$in_cs & !pips$analyzed),!anyNA(region$distinct_cs_snps),
  all(region$distinct_cs_snps>=0),all(region$n_cs>=0))
accepted <- region$status=="ACCEPTED_NONCONVERGED_ITER100"
stopifnot(all(!accepted | (region$trait=="BMI" & region$locus==182L & region$method=="SuSiE_posthoc")))
if("converged" %in% names(region)) stopifnot(all(is.na(region$converged)|region$converged|accepted))
recount <- pips[,.(recomputed_cs_snps=sum(in_cs)),by=.(trait,locus,method)]
region <- merge(region,recount,by=c("trait","locus","method"),all=TRUE,sort=FALSE)
stopifnot(nrow(region)==1156L,!anyNA(region$recomputed_cs_snps),
          all(region$distinct_cs_snps==region$recomputed_cs_snps),
          all((region$n_cs>0)==(region$distinct_cs_snps>0)))
region[,log2_distinct_cs_snps_plus1:=log2(distinct_cs_snps+1)]
panel_a <- region[,.(regions_total=.N,regions_with_cs=sum(distinct_cs_snps>0),
                     total_credible_sets=sum(n_cs)),by=.(trait,method)]
stopifnot(all(panel_a$regions_total==c(BMI=204L,DBP=39L,SBP=46L)[panel_a$trait]))
panel_b <- region[distinct_cs_snps>0,.(trait,locus,method,distinct_cs_snps,
  log2_distinct_cs_snps_plus1,status)]
panel_b_stats <- panel_b[,.(regions_with_cs=.N,median_distinct_cs_snps=as.numeric(median(distinct_cs_snps)),
  q25_distinct_cs_snps=as.numeric(quantile(distinct_cs_snps,.25)),
  q75_distinct_cs_snps=as.numeric(quantile(distinct_cs_snps,.75))),by=.(trait,method)]
cs <- pips[in_cs==TRUE]
regions <- unique(region[,.(trait,locus)])
setorder(regions,trait,locus)
z_parts <- list();input_checks <- list();z_index <- 0L
for(i in seq_len(nrow(regions))) {
  tr <- regions$trait[i];lo <- regions$locus[i]
  input <- file.path(root,"inputs","prepared",sprintf("%s_%03d",tr,lo),"input.rds")
  x <- readRDS(input)
  stopifnot(identical(as.character(x$trait),as.character(tr)),as.integer(x$locus)==as.integer(lo),
    !anyDuplicated(x$union_ids),all(c("EUR","AFR") %in% names(x$summary_stat_list)))
  # NA denotes unavailable; padded zero Z is never used in this statistic.
  zz <- matrix(NA_real_,length(x$union_ids),2L,dimnames=list(x$union_ids,c("EUR","AFR")))
  for(anc in c("EUR","AFR")) {
    ss <- x$summary_stat_list[[anc]];hit <- match(ss$SNP,x$union_ids)
    stopifnot(!anyDuplicated(ss$SNP),!anyNA(hit),all(is.finite(ss$Z)))
    zz[hit,anc] <- ss$Z
  }
  stopifnot(all(rowSums(is.finite(zz))>=1L),identical(is.finite(zz),x$variant_available))
  zmax <- pmax(abs(zz[,"EUR"]),abs(zz[,"AFR"]),na.rm=TRUE)
  stopifnot(all(is.finite(zmax)))
  current <- copy(cs[trait==tr & locus==lo]);hit <- match(current$SNP,x$union_ids)
  stopifnot(!anyNA(hit),identical(as.logical(current$available_EUR),as.logical(is.finite(zz[hit,"EUR"]))),
    identical(as.logical(current$available_AFR),as.logical(is.finite(zz[hit,"AFR"]))))
  if(nrow(current)) {
    z_index <- z_index+1L
    current[,`:=`(Z_EUR_observed=as.numeric(zz[hit,"EUR"]),Z_AFR_observed=as.numeric(zz[hit,"AFR"]),
                  max_abs_observed_Z=as.numeric(zmax[hit]))]
    z_parts[[z_index]] <- current[,.(trait,locus,method,SNP,available_EUR,available_AFR,
                                    Z_EUR_observed,Z_AFR_observed,max_abs_observed_Z)]
  }
  input_checks[[i]] <- data.table(trait=tr,locus=lo,union_variants=length(x$union_ids),
    observed_EUR=sum(is.finite(zz[,"EUR"])),observed_AFR=sum(is.finite(zz[,"AFR"])),
    n_cs_variant_method_records=nrow(current),source_input=input,input_md5=unname(tools::md5sum(input)))
  rm(x,zz);gc(verbose=FALSE)
}
if(length(z_parts)) {
  cs_z <- rbindlist(z_parts)
  stopifnot(!anyDuplicated(cs_z[,.(trait,locus,method,SNP)]),nrow(cs_z)==nrow(cs))
  panel_c <- cs_z[,.(n_cs_variant_records=.N,n_regions_with_cs=uniqueN(locus),
    median_max_abs_observed_Z=median(max_abs_observed_Z)),by=.(trait,method)]
} else {
  cs_z <- data.table(trait=character(),locus=integer(),method=character(),SNP=character(),
    available_EUR=logical(),available_AFR=logical(),Z_EUR_observed=numeric(),Z_AFR_observed=numeric(),max_abs_observed_Z=numeric())
  panel_c <- data.table(trait=character(),method=character(),n_cs_variant_records=integer(),
    n_regions_with_cs=integer(),median_max_abs_observed_Z=numeric())
}
panel_c <- merge(panel_a[,.(trait,method)],panel_c,by=c("trait","method"),all.x=TRUE,sort=FALSE)
panel_c[is.na(n_cs_variant_records),`:=`(n_cs_variant_records=0L,n_regions_with_cs=0L)]
stopifnot(panel_c[,sum(n_cs_variant_records)]==nrow(cs))
table_dir <- file.path(root,"summary","cs_panels")
figure_dir <- file.path(root,"figures","cs_summary")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
tables <- list(panel_a_regions_with_cs=panel_a,panel_b_nonempty_locus_cs_sizes=panel_b,
  panel_b_nonempty_locus_cs_size_statistics=panel_b_stats,panel_c_median_max_abs_Z=panel_c,
  input_checks=rbindlist(input_checks))
for(nm in names(tables)) fwrite(tables[[nm]],file.path(table_dir,paste0(nm,".tsv")),sep="\t",na="NA")
fwrite(cs_z,file.path(table_dir,"cs_member_observed_Z.tsv.gz"),sep="\t",na="NA")
nonconvergence_note <- if(any(accepted)) {
  "The SuSiE result for BMI locus 182 was retained despite nonconvergence after 100 iterations."
} else "These summary tables contain no retained nonconverged results."
} else {
  # Reuse saved numerical summaries without rewriting them.
  table_dir <- file.path(root,"summary","cs_panels")
  figure_dir <- file.path(root,"figures","cs_summary")
  dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
  table_names <- c("panel_a_regions_with_cs","panel_b_nonempty_locus_cs_sizes",
    "panel_b_nonempty_locus_cs_size_statistics","panel_c_median_max_abs_Z","input_checks")
  tables <- setNames(lapply(table_names,function(nm) fread(file.path(table_dir,paste0(nm,".tsv")))),table_names)
  panel_a <- tables$panel_a_regions_with_cs
  panel_b <- tables$panel_b_nonempty_locus_cs_sizes
  panel_c <- tables$panel_c_median_max_abs_Z
  nonconvergence_note <- if(any(fread(region_file,select="status")$status=="ACCEPTED_NONCONVERGED_ITER100")) {
    "The SuSiE result for BMI locus 182 was retained despite nonconvergence after 100 iterations."
  } else "These summary tables contain no retained nonconverged results."
  stopifnot(nrow(panel_a)==12L,nrow(panel_c)==12L,nrow(panel_b)==1058L,
    all(panel_b$distinct_cs_snps>0),
    max(abs(panel_b$log2_distinct_cs_snps_plus1-log2(panel_b$distinct_cs_snps+1)))<1e-12)
}
style <- theme_classic(base_size=9,base_family="sans")+theme(
  axis.text=element_text(size=8,color="black"),axis.title=element_text(size=9),
  legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8),
  plot.title=element_text(size=10,face="bold"),plot.tag=element_text(size=12,face="bold"),
  plot.caption=element_text(size=7,hjust=0),plot.margin=margin(6,7,6,7))
format_plot_data <- function(d) {
  d <- copy(d);d[,method:=factor(method,levels=methods,labels=method_labels)]
  d[,trait:=factor(trait,levels=traits)];d
}
aa <- format_plot_data(panel_a);bb <- format_plot_data(panel_b);cc <- format_plot_data(panel_c)
aa[,trait:=factor(trait,levels=rev(traits))]
aa[,segment_label:=ifelse(regions_with_cs>0,paste0(trait,"\n",regions_with_cs),"")]
aa[,method:=factor(method,levels=rev(method_labels))]
# Method colors stay fixed across all panels. Trait labels provide a redundant
# encoding in panel a even when ggpattern is not available.
use_pattern <- requireNamespace("ggpattern",quietly=TRUE)
p_a <- ggplot(aa,aes(x=regions_with_cs,y=method,fill=method,group=trait))
if(use_pattern) {
  p_a <- p_a+ggpattern::geom_col_pattern(aes(pattern=trait),position="stack",width=.68,
    color="white",linewidth=.25,pattern_colour="white",pattern_fill="white",
    pattern_alpha=.45,pattern_density=.12,pattern_spacing=.04,pattern_key_scale_factor=.5)+
    ggpattern::scale_pattern_manual(values=c(BMI="none",DBP="stripe",SBP="circle"),name="Trait")
} else {
  p_a <- p_a+geom_col(position="stack",width=.68,color="white",linewidth=.25)
}
p_a <- p_a+geom_text(aes(label=segment_label),position=position_stack(vjust=.5),
  size=2.5,color="black",show.legend=FALSE)+scale_fill_manual(values=method_colors,guide="none")+
  scale_x_continuous(expand=expansion(mult=c(0,.03)))+
  labs(title="Regions with credible sets",tag="a",x="Number of regions with at least one credible set",y=NULL)+style
p_b <- ggplot(bb,aes(x=trait,y=log2_distinct_cs_snps_plus1,fill=method))+
  geom_boxplot(width=.75,position=position_dodge2(width=.85,preserve="single"),
    outlier.size=.65,outlier.alpha=.45,linewidth=.3)+
  scale_fill_manual(values=method_colors,drop=FALSE)+
  labs(title="Credible-set size",tag="b",x=NULL,y=expression(log[2]("Distinct credible-set SNPs + 1")))+style
p_c <- ggplot(cc,aes(x=trait,y=median_max_abs_observed_Z,fill=method))+
  geom_col(width=.75,position=position_dodge(width=.8),na.rm=TRUE)+
  geom_text(aes(label=ifelse(is.finite(median_max_abs_observed_Z),sprintf("%.1f",median_max_abs_observed_Z),"")),
    position=position_dodge(width=.8),vjust=-.35,size=2.5,na.rm=TRUE)+
  scale_fill_manual(values=method_colors,drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(0,.12)))+
  labs(title="Association strength of credible-set SNPs",tag="c",x=NULL,
       y="Median maximum absolute observed Z-score")+style
plots <- list(a=p_a,b=p_b,c=p_c)
exports <- list();z <- 0L
for(nm in names(plots)) for(ext in c("pdf","png")) {
  file <- file.path(figure_dir,paste0("panel_",nm,".",ext));z <- z+1L
  ggsave(file,plots[[nm]],width=if(nm=="a")6.2 else 5.5,height=3.8,units="in",dpi=600,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png")
  exports[[z]] <- data.table(panel=nm,file=file)
}
if(requireNamespace("patchwork",quietly=TRUE)) {
  combined <- patchwork::wrap_plots(plots,nrow=1,widths=c(1.15,1,1),guides="collect")+
    patchwork::plot_annotation(caption=paste("BMI: 204 regions; DBP: 39; SBP: 46.",nonconvergence_note))
  for(ext in c("pdf","png")) {
    file <- file.path(figure_dir,paste0("cs_panels_abc.",ext));z <- z+1L
    ggsave(file,combined,width=14,height=4.4,units="in",dpi=600,bg="white",device=if(ext=="pdf")grDevices::cairo_pdf else "png")
    exports[[z]] <- data.table(panel="abc",file=file)
  }
}
saveRDS(list(plots=plots,tables=tables,method_colors=method_colors,
  pattern_enabled=use_pattern,provenance=list(region_metrics=region_file,variant_pips=pip_file,
    summary_md5=tools::md5sum(c(region_file,pip_file)),generated=as.character(Sys.time()),
    session=capture.output(sessionInfo()))),file.path(figure_dir,"cs_panel_objects.rds"))
fwrite(rbindlist(exports),file.path(figure_dir,"figure_manifest.tsv"),sep="\t")
cat("Credible-set panels a-c and numerical source tables saved.\n")
