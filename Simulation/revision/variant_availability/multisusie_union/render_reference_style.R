#!/usr/bin/env Rscript
# Figure labels use the manuscript name; stored result keys remain compatible.
availability_method_labels <- function(x) {
  sub("^MESuSiE(?=([ ~]|$))", "INFUSE\u2212", as.character(x), perl = TRUE)
}

.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Summarize MultiSuSiE union enrichment and plot comparisons with saved method results.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
setDTthreads(1L); options(datatable.optimize=1L)
script <- normalizePath(sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
OUT <- file.path(.infuse_work, "INFUSE_variant_availability/simulation_native_assv/multisusie_union_test_20261020"); NATIVE <- dirname(OUT); ROOT <- dirname(NATIVE); REV <- dirname(ROOT)
a <- as.integer(commandArgs(TRUE)[1]); stopifnot(a %in% 0:2)
stopifnot(startsWith(normalizePath(tempdir()),OUT))
scenario <- paste0("native_ASSV_each_",a)
figdir <- file.path(OUT,"figures",scenario)
tabdir <- file.path(OUT,"reference_style_data",scenario)
dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
oldfig <- file.path(REV,"figures","INFUSE_variant_availability",scenario)
old_paths <- list.files(oldfig,recursive=TRUE,full.names=TRUE)
old_paths <- old_paths[!file.info(old_paths)$isdir]
before <- tools::md5sum(old_paths)
methods <- c("SuSiE","SuSiEx intersect","SuSiEx union","MultiSuSiE intersect",
             "MultiSuSiE union","MESuSiE intersect","MESuSiE union",
             "INFUSE intersect","INFUSE union","INFUSE-H")
palette <- c("SuSiE"="#fbe79e","SuSiEx intersect"="#f2b56e","SuSiEx union"="#d95f3b",
  "MultiSuSiE intersect"="#CC79A7","MultiSuSiE union"="#8C564B",
  "MESuSiE intersect"="#a8ddd3","MESuSiE union"="#3f8c7f",
  "INFUSE intersect"="#a5c8e1","INFUSE union"="#2b6cb0","INFUSE-H"="#6a3d9a")
rename_method <- function(d) {
  d <- as.data.table(d)
  d[,Method:=trimws(as.character(Method))]
  d[Method=="MultiSuSiE",Method:="MultiSuSiE intersect"]
  d
}
env <- new.env(); load(file.path(oldfig,"summarize_low_dim_assv.Rdata"),envir=env)
cs <- rename_method(env$p_size_box$data)
cp <- rename_method(env$p_power_bar$data)
roc <- rename_method(env$p_roc_either$data)
fd <- rename_method(env$p_FDR_Power$data)
fd[,Cat:=as.character(Cat)]
frozen_fd <- copy(fd)
union_cs <- fread(file.path(OUT,"tables","CS_per_cell.tsv"))[assv_each==a & Method=="MultiSuSiE union"]
cs <- rbindlist(list(cs,union_cs),fill=TRUE)
cp <- rbindlist(list(cp,union_cs[,.(Method="MultiSuSiE union",Power_name=mean(Power))]),fill=TRUE)
union_roc <- fread(file.path(OUT,"tables","ROC_per_setting.tsv"))[assv_each==a & Method=="MultiSuSiE union"]
roc <- rbindlist(list(roc,union_roc[,.(Power=mean(Power),FDR=mean(FDR)),by=.(Method,Cutoff)]),fill=TRUE)
union_fd <- fread(file.path(OUT,"tables","nominal_FDR_per_setting.tsv"))[
  assv_each==a & Method=="MultiSuSiE union" & Target_FDR==.05]
fd <- rbindlist(list(fd,union_fd[,.(FDR=.05,Method="MultiSuSiE union",Power=round(mean(Power),2),Cat="Either")]),fill=TRUE)

annotation_columns <- c("non_synonymous_common","Conserved_LindbladToh_common","DHS_Trynka_common",
                        "H3K27ac_Hnisz_common","H3K4me3_Trynka_common")
annotation_order <- c("Non-synonymous","Conserved","DHS","H3K27ac","H3K4me3")
old_detail <- fread(file.path(oldfig,"functional_enrichment_estimates_PIP_gt_0.5.tsv"))
old_detail[method=="MultiSuSiE",method:="MultiSuSiE intersect"]
stopifnot(nrow(old_detail)==270L)
cutoffs <- seq(1,.01,-.01)
# Mirror the reference's exact sequential accumulation and closest-FDR rule.
# For each=0, ancestry strata intentionally retain its legacy /100 normalization.
# Saved n_eligible and effective_cutoff make that legacy convention explicit.
roc_one <- function(pip,signal,target) {
  truth <- signal==target; n <- sum(truth)
  stopifnot(n>0L)
  t(vapply(cutoffs,function(cut) {
    selected <- pip>cut; tp <- sum(selected & truth); ns <- sum(selected)
    c(Cutoff=cut,Power=tp/n,FDR=if(ns) (ns-tp)/ns else 0)
  },numeric(3)))
}
detail_rows <- checks <- strata_rows <- list()
for(b in c(2L,3L,5L)) for(h in 1:2) {
  total_n <- causal_n <- 0L
  total_ann <- causal_ann <- numeric(5)
  selected_n <- c(union=0L,intersect=0L)
  selected_ann <- list(union=numeric(5),intersect=numeric(5))
  accum <- lapply(c("union","intersect"),function(z) lapply(1:3,function(t) matrix(0,100,3)))
  names(accum) <- c("union","intersect"); eligible <- integer(3)
  for(locus in 1:100) {
    rel <- file.path(paste0("base_",b),paste0("h2_",h),paste0("assv_each_",a),sprintf("locus_%03d.rds",locus))
    fit <- readRDS(file.path(OUT,"results",rel))
    old <- readRDS(file.path(NATIVE,"results_native_sim","multisusie_intersect",rel))
    A <- readRDS(file.path(NATIVE,"panels_source_preserved",sprintf("locus_%03d",locus),"annotation.rds"))
    ids <- rownames(A); A <- A[,annotation_columns,drop=FALSE]
    stopifnot(is.numeric(A),all(is.finite(A)),all(A %in% c(0,1)),!anyDuplicated(ids),
              length(fit$pip)==nrow(A),length(old$pip)==nrow(A),fit$metadata$locus==locus)
    ca <- as.data.table(fit$metrics$causal)[Method=="MultiSuSiE union"]
    ii <- match(ca$SNP,ids); stopifnot(!anyNA(ii),!anyDuplicated(ii))
    signal <- integer(length(ids)); signal[ii] <- ca$Signal
    stopifnot(sum(signal!=0)==b+2*a,all(signal[ii] %in% 1:3))
    total_n <- total_n+nrow(A); causal_n <- causal_n+length(ii)
    total_ann <- total_ann+colSums(A); causal_ann <- causal_ann+colSums(A[ii,,drop=FALSE])
    for(m in c("union","intersect")) {
      pip <- if(m=="union") fit$pip else old$pip
      stopifnot(all(is.finite(pip)),all(pip>=-1e-12 & pip<=1+1e-12))
      selected <- pmax(0,pmin(1,pip))>.5
      selected_n[m] <- selected_n[m]+sum(selected)
      selected_ann[[m]] <- selected_ann[[m]]+colSums(A[selected,,drop=FALSE])
      for(target in 1:3) if(any(signal==target))
        accum[[m]][[target]] <- accum[[m]][[target]]+roc_one(pip,signal,target)
    }
    eligible <- eligible+vapply(1:3,function(t) as.integer(any(signal==t)),integer(1))
  }
  true_freq <- causal_ann/causal_n
  background_freq <- (total_ann-causal_ann)/(total_n-causal_n)
  true_enrich <- true_freq/background_freq
  stopifnot(all(is.finite(true_enrich)))
  for(m in c("union","intersect")) {
    ns <- selected_n[m]
    sf <- selected_ann[[m]]/ns
    uf <- (total_ann-selected_ann[[m]])/(total_n-ns)
    estimated <- ifelse(uf>0,sf/uf,NA_real_)
    d <- data.table(assv_each=a,base_causal_count=b,h2_index=h,h2_per_causal=c(1e-4,2e-4)[h],
      method=paste("MultiSuSiE",m),annotation=annotation_order,annotation_column=annotation_columns,pip_cutoff=.5,
      n_variants=total_n,n_causal=causal_n,n_noncausal=total_n-causal_n,
      n_selected=unname(ns),n_unselected=total_n-ns,
      causal_annotation_frequency=true_freq,noncausal_annotation_frequency=background_freq,
      selected_annotation_frequency=sf,unselected_annotation_frequency=uf,
      true_enrichment=true_enrich,estimated_enrichment=estimated,squared_error=(estimated-true_enrich)^2)
    if(m=="union") detail_rows[[length(detail_rows)+1L]] <- d
    else {
      ref <- old_detail[method=="MultiSuSiE intersect" & base_causal_count==b & h2_index==h]
      ref <- ref[match(annotation_order,ref$annotation)]
      numeric_cols <- setdiff(names(d),c("method","annotation","annotation_column"))
      for(k in numeric_cols) {
        error <- max(abs(d[[k]]-ref[[k]]))
        checks[[length(checks)+1L]] <- data.table(base=b,h2=h,field=k,max_error=error)
        stopifnot(is.finite(error),error<1e-8)
      }
    }
    for(target in c(3L,1L,2L)) {
      denom <- if(target==3L) eligible[target] else 100L
      stopifnot(denom>0)
      rr <- as.data.table(accum[[m]][[target]]/denom)
      setnames(rr,c("Cutoff","Power","FDR"))
      jj <- which(abs(rr$FDR-.05)==min(abs(rr$FDR-.05)))
      pick <- jj[which.max(rr$Power[jj])]
      strata_rows[[length(strata_rows)+1L]] <- data.table(Method=paste("MultiSuSiE",m),base=b,h2=h,
        target=target,Cat=if(target==3L) "Shared" else "Ancestry",n_eligible=eligible[target],
        normalization_denominator=denom,Target_FDR=.05,Achieved_FDR=rr$FDR[pick],
        effective_cutoff=rr$Cutoff[pick],Power_unrounded=rr$Power[pick],Power=round(rr$Power[pick],2))
    }
  }
  cat("Scenario ",a,": enrichment and strata finished base=",b," h2=",h,"\n",sep="")
}
strata <- rbindlist(strata_rows)
extra_fd <- strata[Method=="MultiSuSiE union",.(FDR=.05,Power=round(mean(Power),2)),by=.(Method,Cat)]
fd <- rbindlist(list(fd,extra_fd),fill=TRUE)
ref_strata <- strata[Method=="MultiSuSiE intersect",.(recomputed=round(mean(Power),2)),by=.(Method,Cat)]
ref_check <- merge(ref_strata,frozen_fd[,.(Method,Cat,Power)],by=c("Method","Cat"))
stopifnot(nrow(ref_check)==2L,max(abs(ref_check$Power-ref_check$recomputed))<1e-12)
fwrite(ref_check,file.path(tabdir,"legacy_stratum_validation.tsv"),sep="\t")
fwrite(strata,file.path(tabdir,"MultiSuSiE_stratum_operating_points.tsv"),sep="\t")
md <- rbindlist(list(old_detail,rbindlist(detail_rows)),use.names=TRUE)
stopifnot(nrow(md)==300L,all(is.finite(md$squared_error)))
md[,Method:=factor(method,levels=methods)]
md[,annotation:=factor(annotation,levels=annotation_order)]
for(d in list(cs,cp,roc,fd)) d[,Method:=factor(Method,levels=methods)]
fd[,Cat:=factor(Cat,levels=c("Either","Shared","Ancestry"))]
stopifnot(nrow(cs)==6000L,nrow(cp)==10L,nrow(fd)==30L,!anyNA(md$Method))
fwrite(md,file.path(tabdir,"functional_enrichment_estimates_PIP_gt_0.5.tsv"),sep="\t")
fwrite(rbindlist(checks),file.path(tabdir,"enrichment_validation.tsv"),sep="\t")
fwrite(cs,file.path(tabdir,"panel_a_CS_sizes.tsv"),sep="\t")
fwrite(cp,file.path(tabdir,"panel_b_CS_power.tsv"),sep="\t")
fwrite(roc,file.path(tabdir,"panel_c_either_ROC.tsv"),sep="\t")
fwrite(fd,file.path(tabdir,"panel_d_FDR_power.tsv"),sep="\t")
fwrite(data.table(Method=methods,colour=palette[methods]),file.path(tabdir,"method_palette.tsv"),sep="\t")

family <- "Liberation Sans"
theme_word <- function() theme_classic(base_size=7.5,base_family=family)+theme(
 axis.text=element_text(size=6.5,colour="black"),axis.title=element_text(size=8,face="bold"),
 axis.title.x=element_text(margin=margin(t=3)),axis.title.y=element_text(margin=margin(r=3)),
 axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),
 axis.ticks.length=grid::unit(.7,"mm"),strip.background=element_blank(),
 strip.text=element_text(size=7.5,face="bold",margin=margin(3,1,3,1)),
 plot.tag=element_text(size=10,face="plain"),panel.spacing=grid::unit(2,"mm"),
 plot.margin=margin(2,2,2,2),legend.position="none")
fill_scale <- function() scale_fill_manual(values=palette,breaks=methods,drop=FALSE, labels = availability_method_labels)
no_x <- function() theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
med <- cs[,.(median_cs_size=median(Size),label_y=max(log2(Size+1))+.35),by=Method]
med[,label:=sprintf(if(any(median_cs_size%%1!=0)) "%.1f" else "%.0f",median_cs_size)]
fwrite(med,file.path(tabdir,"panel_a_median_labels.tsv"),sep="\t")
p_a <- ggplot(cs,aes(Method,log2(Size+1),fill=Method))+
 geom_boxplot(width=.7,linewidth=.25,outlier.size=.35,colour="darkgray")+
 geom_text(data=med,aes(x=Method,y=label_y,label=label),inherit.aes=FALSE,
   size=6.5/(72.27/25.4),family=family,vjust=0)+
 fill_scale()+scale_y_continuous(expand=expansion(mult=c(.05,.15)))+
 labs(x=NULL,y=expression(bold(log[2]("Set size + 1"))),tag="a")+theme_word()+no_x()
p_b <- ggplot(cp,aes(Method,Power_name,fill=Method))+geom_col(width=.75)+
 geom_text(aes(label=sprintf("%.2f",Power_name)),size=6/(72.27/25.4),vjust=-.4,family=family)+
 fill_scale()+scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(0,.02)))+
 labs(x=NULL,y="Power",tag="b")+theme_word()+no_x()
p_c <- ggplot(roc,aes(Power,1-FDR,colour=Method))+geom_line(linewidth=.3)+
 geom_point(data=roc[abs(Cutoff-.5)<1e-10],shape=21,size=1.5,fill=NA)+
 scale_colour_manual(values=palette,breaks=methods, labels = availability_method_labels)+
 scale_x_continuous(breaks=c(0,.25,.5,.75,1))+scale_y_continuous(breaks=c(0,.25,.5,.75,1))+
 coord_cartesian(xlim=c(0,1),ylim=c(0,1))+labs(x="Recall",y="Precision",tag="c")+theme_word()
p_d <- ggplot(fd,aes(Cat,Power,fill=Method))+geom_col(position=position_dodge(width=.85),width=.8)+
 fill_scale()+scale_x_discrete(labels=c(Either="Either",Shared="Shared",Ancestry="Ancestry-\nspecific"))+
 scale_y_continuous(limits=c(0,.6),breaks=c(0,.2,.4,.6),expand=expansion(mult=c(0,.02)))+
 labs(x="Nominal FDR = 0.05",y="Power",tag="d")+theme_word()
p_e <- ggplot(md,aes(Method,squared_error,fill=Method))+
 geom_hline(yintercept=0,linetype="dashed",linewidth=.25)+
 geom_boxplot(width=.64,linewidth=.25,outlier.size=.45,colour="#333333")+
 facet_wrap(~annotation,scales="free_y",ncol=3)+fill_scale()+
 scale_y_continuous(expand=expansion(mult=c(.03,.1)))+
 labs(x=NULL,y="Squared enrichment error",tag="e")+theme_word()+no_x()
# Keep the reference panel proportions; interpretation notes go below the legend.
combined <- ((p_a|p_c)/(p_b|p_d)/p_e)+plot_layout(heights=c(1,1.1,1.9))
legend_grob <- function() {
  children <- list(); gp <- grid::gpar(fontfamily=family,fontsize=6.6)
  tw <- function(s) grid::convertWidth(grid::grobWidth(grid::textGrob(s,gp=gp)),"mm",valueOnly=TRUE)
  widths <- vapply(availability_method_labels(methods),tw,numeric(1))
  cw <- vapply(1:5,function(j) max(widths[c(j,j+5)])+4,numeric(1))
  total <- sum(cw)+2*4
  cx <- function(x) grid::unit(.5,"npc")+grid::unit(x-total/2,"mm")
  for(row in 1:2) {
    x <- 0; yy <- if(row==1) 12.0 else 7.4
    for(j in 1:5) {
      i <- j+5*(row-1)
      children[[length(children)+1]] <- grid::rectGrob(x=cx(x+1),y=grid::unit(yy,"mm"),
        width=grid::unit(2.4,"mm"),height=grid::unit(2.4,"mm"),gp=grid::gpar(fill=palette[[methods[i]]],col=NA))
      children[[length(children)+1]] <- grid::textGrob(availability_method_labels(methods[i]),x=cx(x+3),y=grid::unit(yy,"mm"),just="left",gp=gp)
      x <- x+cw[j]+2
    }
  }
  children[[length(children)+1]] <- grid::circleGrob(x=grid::unit(.5,"npc")-grid::unit(12,"mm"),
    y=grid::unit(2.7,"mm"),r=grid::unit(.75,"mm"),gp=grid::gpar(fill="white",col="black",lwd=.8))
  children[[length(children)+1]] <- grid::textGrob("PIP = 0.5",x=grid::unit(.5,"npc")-grid::unit(9,"mm"),
    y=grid::unit(2.7,"mm"),just="left",gp=gp)
  stopifnot(total<168)
  grid::gTree(children=do.call(grid::gList,children))
}
makeContent.reference_style <- function(x) {
  plot <- patchworkGrob(x$source_plot)
  lg <- legend_grob()
  notes <- paste0("SuSiEx/MultiSuSiE shared/ancestry-specific bars use overall-PIP stratum recovery, not configuration posteriors.",
    if(a==0L) "\nAncestry-specific panel retains reference normalization (zero-target loci contribute zero; see accompanying notes)." else "")
  ng <- grid::textGrob(notes,x=grid::unit(3,"mm"),y=grid::unit(.6,"npc"),just=c("left","center"),
    gp=grid::gpar(fontfamily=family,fontsize=5.5,col="#444444",lineheight=1.1))
  lay <- grid::grid.layout(3,1,heights=grid::unit.c(grid::unit(1,"null"),grid::unit(16,"mm"),grid::unit(6,"mm")))
  grid::setChildren(x,grid::gList(grid::gTree(children=grid::gList(
    grid::grobTree(plot,vp=grid::viewport(layout.pos.row=1)),
    grid::grobTree(lg,vp=grid::viewport(layout.pos.row=2)),
    grid::grobTree(ng,vp=grid::viewport(layout.pos.row=3))),vp=grid::viewport(layout=lay))))
}
g <- grid::gTree(source_plot=combined,cl="reference_style")
stem <- file.path(figdir,"summarize_with_enrichment_supplement")
for(ext in c("pdf","png")) ggsave(paste0(stem,".",ext),g,width=170,height=170,
  units="mm",dpi=600,bg="white",device=if(ext=="pdf") cairo_pdf else "png",limitsize=FALSE)
# Data and numeric checks are separate from the figure; no original file is written.
stopifnot(identical(before,tools::md5sum(old_paths)))
