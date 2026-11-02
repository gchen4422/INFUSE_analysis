.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Result summary for one native causal-AS-SV level.
# Minimal edit of simulation_v2/5_result_summarizing_assv.R. The plotting,
# pooled-FDR, calibration, credible-set, and output structure are retained; only
# the native RDS adapter, AS-SV-level paths, and native truth semantics differ.
#
# Usage: Rscript 10_result_summarizing_native_assv.R <assv_each:0|1|2> [make_figures]

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1|length(args)>2){
  stop("Usage: Rscript 10_result_summarizing_native_assv.R <assv_each:0|1|2> [make_figures]")
}
assv_each<-as.integer(args[[1]])
if(is.na(assv_each)||!assv_each%in%0:2) stop("assv_each must be 0, 1, or 2")
assv_scenario<-paste0("native_ASSV_each_",assv_each)
make_figures<-if(length(args)==2) as.logical(args[[2]]) else TRUE

module_dir<-paste0(file.path(.infuse_work, "INFUSE_variant_availability/"), "/")
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/common/plot_helpers.R"))
simulation_dir<-module_dir
native_dir<-paste0(module_dir,"simulation_native_assv/")
revision_dir<-dirname(sub("/$","",module_dir))
plot_dir<-paste0(revision_dir,"/figures/INFUSE_variant_availability/",assv_scenario,"/")
dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
res_dir<-paste0(revision_dir,"/res/INFUSE_variant_availability/",assv_scenario,"/")
dir.create(res_dir,recursive=TRUE,showWarnings=FALSE)
res_out<-paste0("variant_availability_",assv_scenario,"_spca.RData")

base_counts<-c(2L,3L,5L)
native_result_file<-function(causal_num,h2,locus){
  file.path(native_dir,"results_native_sim","cells",
            paste0("base_",base_counts[causal_num]),paste0("h2_",h2),
            paste0("assv_each_",assv_each),sprintf("locus_%03d.rds",locus))
}
native_multisusie_file<-function(causal_num,h2,locus){
  file.path(native_dir,"results_native_sim","multisusie_intersect",
            paste0("base_",base_counts[causal_num]),paste0("h2_",h2),
            paste0("assv_each_",assv_each),sprintf("locus_%03d.rds",locus))
}
native_infuse_d187_file<-function(causal_num,h2,locus){
  file.path(native_dir,"results_native_sim","infuse_union_d187",
            paste0("base_",base_counts[causal_num]),paste0("h2_",h2),
            paste0("assv_each_",assv_each),sprintf("locus_%03d.rds",locus))
}

method_levels<-c("SuSiE","SuSiEx intersect","SuSiEx union","MultiSuSiE",
                 "MESuSiE intersect","MESuSiE union","INFUSE intersect","INFUSE union",
                 "INFUSE-H")
method_colors<-c("SuSiE"="#fbe79e","SuSiEx intersect"="#f2b56e",
                 "SuSiEx union"="#d95f3b","MultiSuSiE"="#CC79A7",
                 "MESuSiE intersect"="#a8ddd3","MESuSiE union"="#3f8c7f",
                 "INFUSE intersect"="#a5c8e1","INFUSE union"="#2b6cb0",
                 "INFUSE-H"="#6a3d9a")

# Adapt one native cell to the exact object names used by the V2 summarizer.
# SuSiEx has no SNP-level configuration posterior. Following the original INFUSE
# benchmark treatment of a non-configuration method (PAINTOR), its native overall
# PIP is reused as the shared/EUR/AFR STRATUM RECOVERY SCORE. These columns must
# not be described as SuSiEx configuration posteriors.
adapt_native_cell<-function(x,multisusie,infuse_d187){
  if(length(multisusie$pip)!=nrow(x$truth) || any(!is.finite(multisusie$pip)) ||
     any(multisusie$pip<0 | multisusie$pip>1))
    stop("Malformed MultiSuSiE PIP vector for locus ",x$metadata$locus)
  intersection<-x$availability[,"EUR"] & x$availability[,"AFR"]
  if(any(multisusie$pip[!intersection]!=0) ||
     any(!unlist(multisusie$cs,use.names=FALSE)%in%which(intersection)))
    stop("MultiSuSiE result is not restricted to the intersection at locus ",
         x$metadata$locus)
  if(length(infuse_d187$pip)!=nrow(x$truth) ||
     any(!is.finite(infuse_d187$pip)) ||
     any(infuse_d187$pip<0 | infuse_d187$pip>1) ||
     !identical(infuse_d187$metadata$n_annotations,187L))
    stop("Malformed INFUSE-H result for locus ",x$metadata$locus)
  legacy<-function(r,susiex=FALSE){
    pc<-if(susiex) cbind(r$pip,r$pip,r$pip) else
      cbind(r$pip_EUR,r$pip_AFR,r$pip_shared)
    list(pip=as.numeric(r$pip),pip_config=pc,cs=list(cs=r$cs))
  }
  list(
    zfile=data.frame(SNP=x$truth$SNP,Signal=x$truth$Signal),
    causal_origin=x$truth$causal_origin,
    variant_available=x$availability,
    MESuSiE_res=legacy(x$methods$mesusie_intersect),
    MESuSiE_union_res=legacy(x$methods$mesusie_union),
    INFUSE_res=legacy(x$methods$infuse_intersect),
    INFUSE_union_res=legacy(x$methods$infuse_union),
    INFUSE_H_res=legacy(infuse_d187),
    SuSiEx_intersect_res=legacy(x$methods$susiex_intersect,TRUE),
    SuSiEx_union_res=legacy(x$methods$susiex_union,TRUE),
    MultiSuSiE_res=list(pip=as.numeric(multisusie$pip),
                        pip_config=cbind(multisusie$pip,multisusie$pip,multisusie$pip),
                        cs=list(cs=multisusie$cs)),
    susie_EU=list(pip=x$methods$susie$pip_EUR,sets=list(cs=x$methods$susie$cs)),
    susie_BB=list(pip=x$methods$susie$pip_AFR,sets=list(cs=list()))
  )
}

if(identical(Sys.getenv("NATIVE_ADAPTER_TEST"),"1")){
  f<-native_result_file(1L,1L,1L)
  if(!file.exists(f)) stop("Adapter-test cell is missing: ",f)
  mf<-native_multisusie_file(1L,1L,1L)
  if(!file.exists(mf)) stop("Adapter-test MultiSuSiE cell is missing: ",mf)
  hf<-native_infuse_d187_file(1L,1L,1L)
  if(!file.exists(hf)) stop("Adapter-test INFUSE-H cell is missing: ",hf)
  raw<-readRDS(f); ms<-readRDS(mf); ih<-readRDS(hf); a<-adapt_native_cell(raw,ms,ih)
  stopifnot(identical(a$zfile$SNP,raw$truth$SNP),
            identical(a$causal_origin,raw$truth$causal_origin),
            all(a$SuSiEx_union_res$pip_config[,1]==a$SuSiEx_union_res$pip),
            all(a$SuSiEx_union_res$pip_config[,2]==a$SuSiEx_union_res$pip),
            all(a$SuSiEx_union_res$pip_config[,3]==a$SuSiEx_union_res$pip),
            all(a$MultiSuSiE_res$pip_config[,1]==a$MultiSuSiE_res$pip),
            all(a$MultiSuSiE_res$pip[!raw$availability[,"EUR"]|
                                     !raw$availability[,"AFR"]]==0),
            all(a$INFUSE_union_res$pip_config[,1]==raw$methods$infuse_union$pip_EUR),
            all(a$INFUSE_H_res$pip_config[,1]==ih$pip_EUR),
            all(raw$truth$Signal[raw$truth$causal_origin=="native_ASSV_EUR"]==1L),
            all(raw$truth$Signal[raw$truth$causal_origin=="native_ASSV_AFR"]==2L))
  quit(status=0L)
}


library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)

num_BB = 300000
num_EU = 300000
################################################################
#
#        Result Summarizing
#
#
################################################################
ancestry_all_ROC_data_list<-list()
either_all_ROC_data_list<-list()
shared_all_ROC_data_list<-list()
all_Set_data_list<-list()
pip_res_all_list<-list()
time_all<-c()
data_all<-c()
empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
MESuSiE_ROC_Either_All<-SuSiE_ROC_Either_All<-INFUSE_ROC_Either_All<-INFUSE_union_ROC_Either_All<-INFUSE_H_ROC_Either_All<-MESuSiE_union_ROC_Either_All<-SuSiEx_intersect_ROC_Either_All<-SuSiEx_union_ROC_Either_All<-MultiSuSiE_ROC_Either_All<-empty_mat



# Refuse to summarise a partially finished simulation: figures built from a
# subset of loci silently under-represent the harder settings.  Override with
# ALLOW_INCOMPLETE=1 only for a deliberate mid-run check.
n_expected<-3*2*100
n_present<-sum(vapply(1:3,function(cn) sum(vapply(1:2,function(hh)
  sum(file.exists(vapply(1:100,function(ll) native_result_file(cn,hh,ll),character(1)))),
  numeric(1))),numeric(1)))
cat("result files for",assv_scenario,":",n_present,"/",n_expected,"\n")
n_multisusie_present<-sum(vapply(1:3,function(cn) sum(vapply(1:2,function(hh)
  sum(file.exists(vapply(1:100,function(ll) native_multisusie_file(cn,hh,ll),character(1)))),
  numeric(1))),numeric(1)))
cat("MultiSuSiE result files for",assv_scenario,":",n_multisusie_present,"/",n_expected,"\n")
n_infuse_h_present<-sum(vapply(1:3,function(cn) sum(vapply(1:2,function(hh)
  sum(file.exists(vapply(1:100,function(ll) native_infuse_d187_file(cn,hh,ll),character(1)))),
  numeric(1))),numeric(1)))
cat("INFUSE-H result files for",assv_scenario,":",n_infuse_h_present,"/",n_expected,"\n")
if(n_present<n_expected || n_multisusie_present<n_expected || n_infuse_h_present<n_expected){
  if(!nzchar(Sys.getenv("ALLOW_INCOMPLETE"))){
    stop("Simulation or MultiSuSiE incomplete for ",assv_scenario," (simulation ",
         n_present,"/",n_expected,"; MultiSuSiE ",n_multisusie_present,"/",n_expected,
         "; INFUSE-H ",n_infuse_h_present,"/",n_expected,
         "). Set ALLOW_INCOMPLETE=1 to summarise anyway.")
  }
  warning("Summarising INCOMPLETE results: simulation ",n_present,"/",n_expected,
          "; MultiSuSiE ",n_multisusie_present,"/",n_expected,
          "; INFUSE-H ",n_infuse_h_present,"/",n_expected)
}

compare_table = data.frame(causal_num = NA,h2=NA,LOCI_num = NA,Category = NA)

for(causal_num in 1:3){
  for(h2 in 1:2){
    set = which(file.exists(vapply(1:100,function(ll)
      native_result_file(causal_num,h2,ll),character(1))) &
      file.exists(vapply(1:100,function(ll)
        native_multisusie_file(causal_num,h2,ll),character(1))) &
      file.exists(vapply(1:100,function(ll)
        native_infuse_d187_file(causal_num,h2,ll),character(1))))
    method_objects = c("mesusie_intersect","mesusie_union","infuse_intersect",
                       "infuse_union","susiex_intersect","susiex_union","susie")
    complete_set = set[sapply(set,function(L){
      x<-readRDS(native_result_file(causal_num,h2,L))
      all(sapply(method_objects,function(o) !is.null(x$methods[[o]]) &&
                   length(x$methods[[o]]$pip)>0))
    })]
    set = Reduce(intersect, list(set,complete_set))
    empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
    MESuSiE_ROC<-SuSiE_ROC<-INFUSE_ROC<-INFUSE_union_ROC<-INFUSE_H_ROC<-MESuSiE_union_ROC<-SuSiEx_intersect_ROC<-SuSiEx_union_ROC<-MultiSuSiE_ROC<-empty_mat
    MESuSiE_ROC_Either<-SuSiE_ROC_Either<-INFUSE_ROC_Either<-INFUSE_union_ROC_Either<-INFUSE_H_ROC_Either<-MESuSiE_union_ROC_Either<-SuSiEx_intersect_ROC_Either<-SuSiEx_union_ROC_Either<-MultiSuSiE_ROC_Either<-empty_mat
    MESuSiE_ROC_EU<-MESuSiE_ROC_BB<-SuSiE_ROC_EU<-SuSiE_ROC_BB<-INFUSE_ROC_EU<-INFUSE_ROC_BB<-empty_mat
    INFUSE_union_ROC_EU<-INFUSE_H_ROC_EU<-MESuSiE_union_ROC_EU<-INFUSE_union_ROC_BB<-INFUSE_H_ROC_BB<-MESuSiE_union_ROC_BB<-SuSiEx_intersect_ROC_EU<-SuSiEx_intersect_ROC_BB<-SuSiEx_union_ROC_EU<-SuSiEx_union_ROC_BB<-MultiSuSiE_ROC_EU<-MultiSuSiE_ROC_BB<-empty_mat
    num_shared_region<-num_EUR_region<-num_AFR_region<-0
    set_res_all<-c()
    pip_res_all<-c()

    for(LOCI_num in set){
      a<-adapt_native_cell(readRDS(native_result_file(causal_num,h2,LOCI_num)),
                           readRDS(native_multisusie_file(causal_num,h2,LOCI_num)),
                           readRDS(native_infuse_d187_file(causal_num,h2,LOCI_num)))
      list2env(a,envir=environment())

      # Native AS-SV effects are generated only in the ancestry where the SNP
      # exists. Therefore every causal effect is observable and the full union is
      # the valid truth universe for either and configuration-specific endpoints.
      # There is deliberately no V2 complete-case (`idf`) restriction.
      # `obs` is the complete native union, not a missingness filter.
      obs<-seq_len(nrow(zfile))
      effect_visible<-rep(TRUE,nrow(zfile))
      observable_truth<-as.integer(zfile$Signal!=0)
      restrict_set<-function(x) intersect(unique(as.integer(unlist(x))),obs)

      #Set Size/Power

      MESuSiE_SET = compute_set_susie(restrict_set(MESuSiE_res$cs$cs),observable_truth)
      SuSiE_SET = compute_set_susie(restrict_set(c(susie_EU$sets$cs,susie_BB$sets$cs)),observable_truth)
      INFUSE_SET = compute_set_susie(restrict_set(INFUSE_res$cs$cs),observable_truth)
      INFUSE_union_SET = compute_set_susie(restrict_set(INFUSE_union_res$cs$cs),observable_truth)
      INFUSE_H_SET = compute_set_susie(restrict_set(INFUSE_H_res$cs$cs),observable_truth)
      MESuSiE_union_SET = compute_set_susie(restrict_set(MESuSiE_union_res$cs$cs),observable_truth)
      SuSiEx_intersect_SET = compute_set_susie(restrict_set(SuSiEx_intersect_res$cs$cs),observable_truth)
      SuSiEx_union_SET = compute_set_susie(restrict_set(SuSiEx_union_res$cs$cs),observable_truth)
      MultiSuSiE_SET = compute_set_susie(restrict_set(MultiSuSiE_res$cs$cs),observable_truth)
      #Case 1: Great or equal power, lower or equal FDR
      #Case 2: Great or equal power, greater FDR
      #Case 3: Lower power,  lower or equal FDR
      #Case 4: Lower power,  greater FDR
      if (INFUSE_union_SET[1] >= INFUSE_SET[1] & INFUSE_union_SET[2] <= INFUSE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 1)
      } else if (INFUSE_union_SET[1] >= INFUSE_SET[1] & INFUSE_union_SET[2] > INFUSE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 2)
      } else if (INFUSE_union_SET[1] < INFUSE_SET[1] & INFUSE_union_SET[2] <= INFUSE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 3)
      } else {
        locus_summ = c(causal_num, h2, LOCI_num, 4)
      }
      compare_table = rbind(compare_table,locus_summ)

      set_res_all<-rbind(set_res_all,rbind(c(MESuSiE_SET,"MESuSiE intersect"),c( SuSiE_SET,"SuSiE"),c(INFUSE_SET,"INFUSE intersect"),c(INFUSE_union_SET,"INFUSE union"),c(INFUSE_H_SET,"INFUSE-H"),c(MESuSiE_union_SET,"MESuSiE union"),c(SuSiEx_intersect_SET,"SuSiEx intersect"),c(SuSiEx_union_SET,"SuSiEx union"),c(MultiSuSiE_SET,"MultiSuSiE")))

      #ROC for shared/either/ancestry-specific
      if(sum(zfile$Signal[obs]==3)>0){
        num_shared_region<-num_shared_region+1
        MESuSiE_ROC = MESuSiE_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        SuSiE_ROC = SuSiE_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = pmin(susie_EU$pip,susie_BB$pip)[obs]),ifelse(zfile$Signal[obs]==3,1,0))
        INFUSE_ROC = INFUSE_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        INFUSE_union_ROC = INFUSE_union_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_union_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        INFUSE_H_ROC = INFUSE_H_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_H_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        MESuSiE_union_ROC = MESuSiE_union_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_union_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        SuSiEx_intersect_ROC = SuSiEx_intersect_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_intersect_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        SuSiEx_union_ROC = SuSiEx_union_ROC + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_union_res$pip_config[obs,3]),ifelse(zfile$Signal[obs]==3,1,0))
        MultiSuSiE_ROC = MultiSuSiE_ROC + compute_ROC(data.frame("SNP"=zfile$SNP[obs],"PIP"=MultiSuSiE_res$pip[obs]),ifelse(zfile$Signal[obs]==3,1,0))
      }

      MESuSiE_ROC_Either = MESuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_res$pip[obs]),observable_truth[obs])
      SuSiE_ROC_Either = SuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = pmax(susie_EU$pip,susie_BB$pip)[obs]),observable_truth[obs])
      INFUSE_ROC_Either = INFUSE_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_res$pip[obs]),observable_truth[obs])
      INFUSE_union_ROC_Either = INFUSE_union_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_union_res$pip[obs]),observable_truth[obs])
      INFUSE_H_ROC_Either = INFUSE_H_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_H_res$pip[obs]),observable_truth[obs])
      MESuSiE_union_ROC_Either = MESuSiE_union_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_union_res$pip[obs]),observable_truth[obs])
      SuSiEx_intersect_ROC_Either = SuSiEx_intersect_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_intersect_res$pip[obs]),observable_truth[obs])
      SuSiEx_union_ROC_Either = SuSiEx_union_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_union_res$pip[obs]),observable_truth[obs])
      MultiSuSiE_ROC_Either = MultiSuSiE_ROC_Either + compute_ROC(data.frame("SNP"=zfile$SNP[obs],"PIP"=MultiSuSiE_res$pip[obs]),observable_truth[obs])

      MESuSiE_ROC_Either_All = MESuSiE_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_res$pip[obs]),observable_truth[obs])
      SuSiE_ROC_Either_All = SuSiE_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = pmax(susie_EU$pip,susie_BB$pip)[obs]),observable_truth[obs])
      INFUSE_ROC_Either_All = INFUSE_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_res$pip[obs]),observable_truth[obs])
      INFUSE_union_ROC_Either_All = INFUSE_union_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_union_res$pip[obs]),observable_truth[obs])
      INFUSE_H_ROC_Either_All = INFUSE_H_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_H_res$pip[obs]),observable_truth[obs])
      MESuSiE_union_ROC_Either_All = MESuSiE_union_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_union_res$pip[obs]),observable_truth[obs])
      SuSiEx_intersect_ROC_Either_All = SuSiEx_intersect_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_intersect_res$pip[obs]),observable_truth[obs])
      SuSiEx_union_ROC_Either_All = SuSiEx_union_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_union_res$pip[obs]),observable_truth[obs])
      MultiSuSiE_ROC_Either_All = MultiSuSiE_ROC_Either_All + compute_ROC(data.frame("SNP"=zfile$SNP[obs],"PIP"=MultiSuSiE_res$pip[obs]),observable_truth[obs])

      if(sum(zfile$Signal[obs]==1)>0){
        num_EUR_region<-num_EUR_region+1
        MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC_uni_susie(data.frame("SNP" =zfile$SNP[obs] ,"PIP_1" =susie_EU$pip[obs] ,"PIP_2" = susie_BB$pip[obs]),ifelse(zfile$Signal[obs]==1,1,0))
        INFUSE_ROC_EU = INFUSE_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        INFUSE_union_ROC_EU = INFUSE_union_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_union_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        INFUSE_H_ROC_EU = INFUSE_H_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_H_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        MESuSiE_union_ROC_EU = MESuSiE_union_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_union_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        SuSiEx_intersect_ROC_EU = SuSiEx_intersect_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_intersect_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        SuSiEx_union_ROC_EU = SuSiEx_union_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_union_res$pip_config[obs,1]),ifelse(zfile$Signal[obs]==1,1,0))
        MultiSuSiE_ROC_EU = MultiSuSiE_ROC_EU + compute_ROC(data.frame("SNP"=zfile$SNP[obs],"PIP"=MultiSuSiE_res$pip[obs]),ifelse(zfile$Signal[obs]==1,1,0))
      }

      if(sum(zfile$Signal[obs]==2)>0){
        num_AFR_region<-num_AFR_region+1
        MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC_uni_susie(data.frame("SNP" =zfile$SNP[obs] ,"PIP_2" =susie_EU$pip[obs] ,"PIP_1" = susie_BB$pip[obs]),ifelse(zfile$Signal[obs]==2,1,0))
        INFUSE_ROC_BB = INFUSE_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        INFUSE_union_ROC_BB = INFUSE_union_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_union_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        INFUSE_H_ROC_BB = INFUSE_H_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = INFUSE_H_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        MESuSiE_union_ROC_BB = MESuSiE_union_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = MESuSiE_union_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        SuSiEx_intersect_ROC_BB = SuSiEx_intersect_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_intersect_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        SuSiEx_union_ROC_BB = SuSiEx_union_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP[obs] ,"PIP" = SuSiEx_union_res$pip_config[obs,2]),ifelse(zfile$Signal[obs]==2,1,0))
        MultiSuSiE_ROC_BB = MultiSuSiE_ROC_BB + compute_ROC(data.frame("SNP"=zfile$SNP[obs],"PIP"=MultiSuSiE_res$pip[obs]),ifelse(zfile$Signal[obs]==2,1,0))
      }

      ##Sumstat
      locus_sumstat<-data.frame("RSID" = zfile$SNP,"Signal" = zfile$Signal,
                                "causal_origin" = causal_origin,
                                "effect_visible" = effect_visible,
                                "stratum" = causal_origin)%>%mutate(
                                                 MESuSiE_Either = MESuSiE_res$pip , MESuSiE_Shared = MESuSiE_res$pip_config[,3],MESuSiE_WB = MESuSiE_res$pip_config[,1],MESuSiE_BB = MESuSiE_res$pip_config[,2],
                                                 INFUSE_Either = INFUSE_res$pip , INFUSE_Shared = INFUSE_res$pip_config[,3],INFUSE_WB = INFUSE_res$pip_config[,1],INFUSE_BB = INFUSE_res$pip_config[,2],
                                                 INFUSE_union_Either = INFUSE_union_res$pip , INFUSE_union_Shared = INFUSE_union_res$pip_config[,3],INFUSE_union_WB = INFUSE_union_res$pip_config[,1],INFUSE_union_BB = INFUSE_union_res$pip_config[,2],
                                                 INFUSE_H_Either = INFUSE_H_res$pip , INFUSE_H_Shared = INFUSE_H_res$pip_config[,3],INFUSE_H_WB = INFUSE_H_res$pip_config[,1],INFUSE_H_BB = INFUSE_H_res$pip_config[,2],
                                                 MESuSiE_union_Either = MESuSiE_union_res$pip , MESuSiE_union_Shared = MESuSiE_union_res$pip_config[,3],MESuSiE_union_WB = MESuSiE_union_res$pip_config[,1],MESuSiE_union_BB = MESuSiE_union_res$pip_config[,2],
                                                 SuSiEx_intersect_Either = SuSiEx_intersect_res$pip , SuSiEx_intersect_Shared = SuSiEx_intersect_res$pip_config[,3],SuSiEx_intersect_WB = SuSiEx_intersect_res$pip_config[,1],SuSiEx_intersect_BB = SuSiEx_intersect_res$pip_config[,2],
                                                 SuSiEx_union_Either = SuSiEx_union_res$pip , SuSiEx_union_Shared = SuSiEx_union_res$pip_config[,3],SuSiEx_union_WB = SuSiEx_union_res$pip_config[,1],SuSiEx_union_BB = SuSiEx_union_res$pip_config[,2],
                                                 MultiSuSiE_Either = MultiSuSiE_res$pip, MultiSuSiE_Shared = MultiSuSiE_res$pip, MultiSuSiE_WB = MultiSuSiE_res$pip, MultiSuSiE_BB = MultiSuSiE_res$pip,
                                                 SuSiE_Either = pmax(susie_EU$pip,susie_BB$pip),SuSiE_Shared = pmin(susie_EU$pip,susie_BB$pip),SuSiE_WB = susie_EU$pip,SuSiE_BB =susie_BB$pip,N_WB = num_EU,N_BB = num_BB,h2 = h2,causal_num  = causal_num,locus = LOCI_num )

      data_all<-rbind(data_all,locus_sumstat)
      cat(LOCI_num)
    }

    #######################################
    #
    #   credible set combined
    #
    #######################################
    set_res_all_data = data.frame("Power" = as.numeric(set_res_all[,1]),"FDR" = as.numeric(set_res_all[,2]),"Size"=as.numeric(set_res_all[,4]),"Method"=set_res_all[,5])
    set_res_all_data$h2 = h2
    set_res_all_data$causal_num = causal_num
    set_res_summary<-set_res_all_data %>% group_by(Method)%>%summarise_at(vars(Power,Size),list(name = mean))

    set_res_summary$h2 = h2
    set_res_summary$causal_num = causal_num
    all_Set_data_list[[paste0("causal_",causal_num,"_h2_",h2)]]<-set_res_all_data

    #######################################
    #
    #   ROC combined
    #
    #######################################
    # Either signal part
    roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either,
                          SuSiE_ROC_Either = SuSiE_ROC_Either,
                          INFUSE_ROC_Either = INFUSE_ROC_Either,
                          INFUSE_union_ROC_Either = INFUSE_union_ROC_Either,
                          INFUSE_H_ROC_Either = INFUSE_H_ROC_Either,
                          MESuSiE_union_ROC_Either = MESuSiE_union_ROC_Either,
                          SuSiEx_intersect_ROC_Either = SuSiEx_intersect_ROC_Either,
                          SuSiEx_union_ROC_Either = SuSiEx_union_ROC_Either,
                          MultiSuSiE_ROC_Either = MultiSuSiE_ROC_Either)
    method_names <- c("MESuSiE intersect", "SuSiE", "INFUSE intersect","INFUSE union","INFUSE-H","MESuSiE union","SuSiEx intersect","SuSiEx union","MultiSuSiE")
    either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, length(set), h2, causal_num)

    # Shared signal part
    roc_data_list <- list(MESuSiE_ROC = MESuSiE_ROC, SuSiE_ROC = SuSiE_ROC, INFUSE_ROC = INFUSE_ROC, INFUSE_union_ROC = INFUSE_union_ROC, INFUSE_H_ROC = INFUSE_H_ROC, MESuSiE_union_ROC = MESuSiE_union_ROC, SuSiEx_intersect_ROC = SuSiEx_intersect_ROC, SuSiEx_union_ROC = SuSiEx_union_ROC, MultiSuSiE_ROC = MultiSuSiE_ROC)
    method_names <- c("MESuSiE intersect", "SuSiE", "INFUSE intersect","INFUSE union","INFUSE-H","MESuSiE union","SuSiEx intersect","SuSiEx union","MultiSuSiE")
    shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, num_shared_region, h2, causal_num)

    # Ancestry specific part
    roc_data_list <- list(MESuSiE_ROC_EU = MESuSiE_ROC_EU,MESuSiE_ROC_BB = MESuSiE_ROC_BB,
                          SuSiE_ROC_EU = SuSiE_ROC_EU, SuSiE_ROC_BB = SuSiE_ROC_BB,
                          INFUSE_ROC_EU = INFUSE_ROC_EU,INFUSE_ROC_BB = INFUSE_ROC_BB,
                          INFUSE_union_ROC_EU = INFUSE_union_ROC_EU,INFUSE_union_ROC_BB = INFUSE_union_ROC_BB,
                          INFUSE_H_ROC_EU = INFUSE_H_ROC_EU,INFUSE_H_ROC_BB = INFUSE_H_ROC_BB,
                          MESuSiE_union_ROC_EU = MESuSiE_union_ROC_EU,MESuSiE_union_ROC_BB = MESuSiE_union_ROC_BB,
                          SuSiEx_intersect_ROC_EU = SuSiEx_intersect_ROC_EU, SuSiEx_intersect_ROC_BB = SuSiEx_intersect_ROC_BB,
                          SuSiEx_union_ROC_EU = SuSiEx_union_ROC_EU, SuSiEx_union_ROC_BB = SuSiEx_union_ROC_BB,
                          MultiSuSiE_ROC_EU = MultiSuSiE_ROC_EU, MultiSuSiE_ROC_BB = MultiSuSiE_ROC_BB)
    method_names <- c("MESuSiE intersect WB","MESuSiE intersect BB",
                      "SuSiE WB","SuSiE BB",
                      "INFUSE intersect WB","INFUSE intersect BB",
                      "INFUSE union WB","INFUSE union BB",
                      "INFUSE-H WB","INFUSE-H BB",
                      "MESuSiE union WB","MESuSiE union BB",
                      "SuSiEx intersect WB","SuSiEx intersect BB",
                      "SuSiEx union WB","SuSiEx union BB",
                      "MultiSuSiE WB","MultiSuSiE BB")
    ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, length(set), h2, causal_num)


  }
}


#######################################
#
#  Loci performance summary
#
####################################### 
compare_table = compare_table[-1,]
summary_table <- compare_table %>%
  group_by(Category) %>%
  summarise(
    Count = n()  # Count of rows per categor
  )

summary_table
#######################################
#
#  Set Size and Power
#
####################################### 

all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
all_Set_data_dataframe<-all_Set_data_dataframe%>%mutate(Method = fct_relevel(Method, method_levels))
#levels(all_Set_data_dataframe$Method)<-c("MeSuSiE","MESuSiE_fun_","SuSiE")
cols <- c("Method","h2","causal_num")
set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
colnames(set_power_summary)[4]<-"Power_name"

set_power_summary<-set_power_summary%>%mutate(Method = fct_relevel(Method, method_levels))
#levels(set_power_summary$Method)<-c("MESuSiE intersect","SuSiE","MESuSiE_fun_")
set_power_summary
set_size_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars( Size),list(name = mean)))
colnames(set_size_summary)[4]<-" Size_name"
set_size_summary

#######################################
#
#   ROC data 
#
####################################### 
shared_all_ROC_data_dataframe<-ROC_list_to_data(shared_all_ROC_data_list)
either_all_ROC_data_dataframe<-ROC_list_to_data(either_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ROC_list_to_data(ancestry_all_ROC_data_list)
# Method names already carry the " WB"/" BB" suffix, so no recoding is needed.

########################################################################
#
#
#FDR Power based on threshold of 0.01,0.05,0.1,0.5
#
#
########################################################################

FDR_Power_shared<-FDR_Power(shared_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, method_levels))
FDR_Power_shared%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))


FDR_Power_either<-FDR_Power(either_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, method_levels))
FDR_Power_either%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

FDR_Power_ancestry<-FDR_Power(ancestry_all_ROC_data_dataframe)%>%
  mutate(Method = fct_relevel(Method, c(
    "SuSiE BB","SuSiE WB",
    "SuSiEx intersect BB","SuSiEx intersect WB",
    "SuSiEx union BB","SuSiEx union WB",
    "MultiSuSiE BB","MultiSuSiE WB",
    "MESuSiE intersect BB","MESuSiE intersect WB",
    "MESuSiE union BB","MESuSiE union WB",
    "INFUSE intersect BB","INFUSE intersect WB",
    "INFUSE union BB","INFUSE union WB",
    "INFUSE-H BB","INFUSE-H WB")))
FDR_Power_ancestry%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

####################################### 
#
#  Which signal targets exist in this scenario
#
#  Random base configurations can leave a particular causal type absent in a
#  finite group of loci. Empty target panels are skipped rather than displayed as
#  apparent method failure. Native AS-SVs themselves are never removed from the
#  EUR/AFR endpoints.
#######################################
has_targets <- function(roc_df, fdr_df) {
  if (is.null(roc_df) || !nrow(roc_df)) return(FALSE)
  if (all(is.na(roc_df$Power)) || all(!is.finite(roc_df$Power))) return(FALSE)
  # An identically-zero Power curve means the positive set was empty.
  if (all(roc_df$Power == 0, na.rm = TRUE)) return(FALSE)
  if (is.null(fdr_df) || !nrow(fdr_df)) return(FALSE)
  if (all(is.na(fdr_df$Power))) return(FALSE)
  TRUE
}
have_shared   <- has_targets(shared_all_ROC_data_dataframe, FDR_Power_shared)
have_ancestry <- has_targets(ancestry_all_ROC_data_dataframe, FDR_Power_ancestry)
cat("targets present -> shared:", have_shared, " ancestry:", have_ancestry, "\n")




#######################################
#
#   save the result
#
####################################### 

# Causal-only recovery, stratified by the exact causal origin. AS-CV and native
# AS-SV use the SAME truth-configuration score: Signal 1 -> WB/EUR, Signal 2 ->
# BB/AFR, Signal 3 -> shared. This table assesses recovery within each origin; the
# pooled ROC/calibration objects above assess classification against the full union.
origin_methods<-c("MESuSiE","SuSiE","SuSiEx_intersect","SuSiEx_union","MultiSuSiE",
                  "INFUSE","INFUSE_union","INFUSE_H","MESuSiE_union")
origin_labels<-c("MESuSiE intersect","SuSiE","SuSiEx intersect","SuSiEx union","MultiSuSiE",
                 "INFUSE intersect","INFUSE union","INFUSE-H","MESuSiE union")
origin_config_recovery<-rbindlist(lapply(seq_along(origin_methods),function(i){
  p<-origin_methods[i]
  score<-ifelse(data_all$Signal==1,data_all[[paste0(p,"_WB")]],
         ifelse(data_all$Signal==2,data_all[[paste0(p,"_BB")]],
         ifelse(data_all$Signal==3,data_all[[paste0(p,"_Shared")]],NA_real_)))
  data.table(causal_num=data_all$causal_num,h2=data_all$h2,locus=data_all$locus,
             RSID=data_all$RSID,Signal=data_all$Signal,
             causal_origin=data_all$causal_origin,Method=origin_labels[i],
             truth_config_PIP=score,
             score_semantics=if(grepl("^(SuSiEx|MultiSuSiE)",origin_labels[i]))
               "native overall PIP reused as truth-stratum recovery score" else
               if(origin_labels[i]=="SuSiE") "published ancestry-marginal heuristic" else
               "native SNP-level configuration PIP")
}))[Signal!=0]
origin_config_summary<-origin_config_recovery[,.(
  N=.N,mean_PIP=mean(truth_config_PIP),median_PIP=median(truth_config_PIP),
  prop_gt_05=mean(truth_config_PIP>0.5),prop_gt_09=mean(truth_config_PIP>0.9),
  prop_zero=mean(truth_config_PIP==0)),
  by=.(causal_num,h2,causal_origin,Method,score_semantics)]
fwrite(origin_config_summary,paste0(res_dir,"causal_origin_config_recovery.tsv"),sep="\t")

save(either_all_ROC_data_dataframe,shared_all_ROC_data_dataframe,
     ancestry_all_ROC_data_dataframe,FDR_Power_either,FDR_Power_shared,
     FDR_Power_ancestry,all_Set_data_dataframe,set_power_summary,time_all,data_all,
     origin_config_recovery,origin_config_summary,compare_table,summary_table,
     file=paste0(res_dir,res_out))

################################################################
#
#        Numbers to report
#
#
################################################################


##95% Credible set size/power/fdr
all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))
all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))


#Power with a FDR = 0.05 
FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))



FDR_Power_either%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))



#FDR and Power with a PIP threshold of 0.5	
shared_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
ancestry_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))




####################################### 
#
#  Layout adjustments for eight methods
#
#  The template's plotting functions and canvas sizes were tuned for five short
#  method names.  Eight longer names ("MESuSiE intersect", "MultiSuSiE", ...)
#  collide on the x axis and overflow a single-row legend, so the legend is
#  wrapped and the fonts/margins are relaxed.  Nothing in utility_assv.R changes.
#######################################
fix_legend <- function(p, nrow = 2) {
  p + guides(fill = guide_legend(nrow = nrow, byrow = TRUE),
             colour = guide_legend(nrow = nrow, byrow = TRUE),
             color = guide_legend(nrow = nrow, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6),
          legend.key.size = unit(3.2, "mm"),
          legend.box.margin = margin(t = -2, b = 0),
          legend.margin = margin(t = 0, b = 0))
}
safe_ggsave <- function(path, plot, ...) {
  if (is.null(plot)) { message("skipping figure (no data): ", basename(path)); return(invisible()) }
  ggsave(path, plot, ...)
}
fix_xaxis <- function(p) {
  p + scale_x_discrete(labels = availability_method_labels) +
    theme(axis.text.x = element_text(size = 5.6, angle = 30, hjust = 1))
}

#######################################
#
#  PIP calibration
#
####################################### 

####PIP Either

PIP_calibration_either<-create_obs_frq_byh2(data_all%>%filter(effect_visible)%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiEx_intersect_Either, SuSiEx_union_Either, MultiSuSiE_Either, INFUSE_Either, INFUSE_union_Either, INFUSE_H_Either, MESuSiE_union_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiEx_intersect_Either","SuSiEx_union_Either","MultiSuSiE_Either","INFUSE_Either","INFUSE_union_Either","INFUSE_H_Either","MESuSiE_union_Either"))
PIP_calibration_either<-PIP_calibration_either%>%mutate(Method = fct_recode(Method, "MESuSiE intersect" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiEx intersect" = "SuSiEx_intersect_Either", "SuSiEx union" = "SuSiEx_union_Either","MultiSuSiE"="MultiSuSiE_Either","INFUSE intersect" = "INFUSE_Either","INFUSE union" = "INFUSE_union_Either","INFUSE-H" = "INFUSE_H_Either","MESuSiE union" = "MESuSiE_union_Either"))%>%mutate(Method = fct_relevel(Method,method_levels))

levels(PIP_calibration_either$Method)<-gsub(" ","~",levels(PIP_calibration_either$Method))
levels(PIP_calibration_either$Method)<-sub("^MESuSiE~", 'INFUSE*"\u2212"~', levels(PIP_calibration_either$Method))
levels(PIP_calibration_either$Method)[levels(PIP_calibration_either$Method)=="INFUSE-H"]<-"INFUSE*\"-H\""
p_calibration_either_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either)


# plot_dir is already set at the top of this script for this scenario
if(make_figures) safe_ggsave(paste0(plot_dir,"pip_either_calibration.pdf"),p_calibration_either_byh2,height=220, width=400, units = "mm",dpi=500)


####PIP Shared	
# Calibrate shared-configuration PIPs against shared causal signals (Signal == 3).
PIP_calibration_Shared<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared, SuSiEx_intersect_Shared, SuSiEx_union_Shared, MultiSuSiE_Shared, INFUSE_Shared, INFUSE_union_Shared, INFUSE_H_Shared, MESuSiE_union_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","SuSiEx_intersect_Shared","SuSiEx_union_Shared","MultiSuSiE_Shared","INFUSE_Shared","INFUSE_union_Shared","INFUSE_H_Shared","MESuSiE_union_Shared"))
PIP_calibration_Shared<-PIP_calibration_Shared%>%mutate(Method = fct_recode(Method, "MESuSiE intersect" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiEx intersect" = "SuSiEx_intersect_Shared", "SuSiEx union" = "SuSiEx_union_Shared","MultiSuSiE"="MultiSuSiE_Shared","INFUSE intersect" = "INFUSE_Shared","INFUSE union" = "INFUSE_union_Shared","INFUSE-H" = "INFUSE_H_Shared","MESuSiE union" = "MESuSiE_union_Shared"))%>%mutate(Method = fct_relevel(Method,method_levels))

levels(PIP_calibration_Shared$Method)<-gsub(" ","~",levels(PIP_calibration_Shared$Method))
levels(PIP_calibration_Shared$Method)<-sub("^MESuSiE~", 'INFUSE*"\u2212"~', levels(PIP_calibration_Shared$Method))
levels(PIP_calibration_Shared$Method)[levels(PIP_calibration_Shared$Method)=="INFUSE-H"]<-"INFUSE*\"-H\""
if(have_shared){
  p_calibration_shared_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_Shared)
} else {
  p_calibration_shared_byh2<-NULL
  message("no identifiable shared targets; skipping the shared PIP-calibration panel")
}


###PIP ancestry
PIP_calibration_WB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_WB,SuSiE_WB,SuSiEx_intersect_WB,SuSiEx_union_WB,MultiSuSiE_WB,INFUSE_WB,INFUSE_union_WB,INFUSE_H_WB,MESuSiE_union_WB),c(1),c("MESuSiE_WB","SuSiE_WB","SuSiEx_intersect_WB","SuSiEx_union_WB","MultiSuSiE_WB","INFUSE_WB","INFUSE_union_WB","INFUSE_H_WB","MESuSiE_union_WB"))
PIP_calibration_BB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_BB,SuSiE_BB,SuSiEx_intersect_BB,SuSiEx_union_BB,MultiSuSiE_BB,INFUSE_BB,INFUSE_union_BB,INFUSE_H_BB,MESuSiE_union_BB),c(2),c("MESuSiE_BB","SuSiE_BB","SuSiEx_intersect_BB","SuSiEx_union_BB","MultiSuSiE_BB","INFUSE_BB","INFUSE_union_BB","INFUSE_H_BB","MESuSiE_union_BB"))

PIP_calibration_ancestry<-rbind(PIP_calibration_WB,PIP_calibration_BB)
anc_levels<-c("MESuSiE_WB","MESuSiE_BB","MESuSiE_union_WB","MESuSiE_union_BB","SuSiE_WB","SuSiE_BB","SuSiEx_intersect_WB","SuSiEx_intersect_BB","SuSiEx_union_WB","SuSiEx_union_BB","MultiSuSiE_WB","MultiSuSiE_BB","INFUSE_WB","INFUSE_BB","INFUSE_union_WB","INFUSE_union_BB","INFUSE_H_WB","INFUSE_H_BB")
anc_levels<-anc_levels[anc_levels%in%unique(as.character(PIP_calibration_ancestry$Method))]
PIP_calibration_ancestry<-PIP_calibration_ancestry%>%mutate(Method = factor(as.character(Method),levels = anc_levels))
# label_parsed needs each label to be a valid R expression, so build them with
# tildes and no spaces; the ancestry suffix is resolved before "_union".
anc_labels<-gsub("_WB$","~White~British",anc_levels)
anc_labels<-gsub("_BB$","~Black~British",anc_labels)
anc_labels<-gsub("_union","~union",anc_labels)
anc_labels<-gsub("^INFUSE_H~","INFUSE*\"-H\"~",anc_labels)
anc_labels<-sub("^MESuSiE~union~", 'INFUSE*"\u2212"~union~', anc_labels)
anc_labels<-sub("^MESuSiE~", 'INFUSE*"\u2212"~intersect~', anc_labels)
levels(PIP_calibration_ancestry$Method)<-anc_labels

if(have_ancestry){
  p_calibration_ancestry<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry)
} else {
  p_calibration_ancestry<-NULL
  message("no identifiable ancestry-specific targets; skipping the ancestry PIP-calibration panel")
}

################################################################
#
#        Plot
#
#
################################################################
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/common/plot_helpers.R"))
upper_limit<-round(all_Set_data_dataframe%>%filter(Method == "SuSiE",causal_num=="Num~Causal  == 5 ",h2=="~h^2 == 10^-4")%>%summarise(upper = quantile(Size,0.75))%>%pull(upper))
p_size_box<-Set_Size_fun(all_Set_data_dataframe%>%mutate(Size= log2(Size+1)),upper_limit =log2(upper_limit+1))+
  scale_fill_manual(values=method_colors,guide="none", labels = availability_method_labels)
p_size_box<-fix_xaxis(p_size_box+ ylab("log2(Set Size + 1)"))
p_size_box

p_power_bar<-fix_xaxis(Set_Power_fun(set_power_summary)+
                         scale_fill_manual(values=method_colors,guide="none", labels = availability_method_labels))
size_power<-p_size_box/p_power_bar+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 7,face="bold"))

if(make_figures) safe_ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=190, width=360, units = "mm",dpi=500)

##########################################################
#
#   Either Signal 
#  ROC | FDR Power 
#
##########################################################

###################
#
#ROC
#
###################
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,method_levels))
p_ROC_Either<-fix_legend(ROC_shared_fun(either_all_ROC_data_dataframe)+
                           scale_color_manual(values=method_colors, labels = availability_method_labels))
###################
#
#FDR&Power
#
###################

p_FDR_Power_Either<-fix_legend(FDR_Power_shared_fun(FDR_Power_either%>%filter(FDR!=0.01)%>%mutate(Method = fct_relevel(Method,method_levels)))+
                                 scale_fill_manual(values=method_colors, labels = availability_method_labels)+
                                 ylim(0,min(ceiling(max(FDR_Power_either$Power)*10)/10+0.1,1)))


ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
if(make_figures) safe_ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Either.pdf"),ROC_FDR_Power_Either,height=280, width=340, units = "mm",dpi=500)
##########################################################
#
# Shared Signal 
#  ROC | FDR Power 
#
##########################################################
###################
#
#ROC
# 
###################
if(have_shared){
  shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,method_levels))
  p_ROC_shared<-fix_legend(ROC_shared_fun(shared_all_ROC_data_dataframe)+
                             scale_color_manual(values=method_colors, labels = availability_method_labels))
  ###################
  #
  #FDR&Power
  #
  ###################
  p_FDR_Power_shared<-fix_legend(FDR_Power_shared_fun(FDR_Power_shared%>%filter(FDR!=0.01)%>%mutate(Method = fct_relevel(Method,method_levels)))+
                                   scale_fill_manual(values=method_colors, labels = availability_method_labels)+
                                   ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1)))

  ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
  if(make_figures) safe_ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=280, width=340, units = "mm",dpi=500)
} else {
  ROC_FDR_Power_shared <- NULL
  message("no identifiable shared targets in ", assv_scenario, "; skipping the Shared panel")
}


##########################################################
#
# Ancestry-specific Signal 
#  ROC | FDR Power | PIP calibration  
#
##########################################################
###################
#
#ROC
#
###################

if(have_ancestry){
  ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>% 
    mutate(Method = as.character(Method))
  # Method names can contain a space ("INFUSE union"), so split off the FINAL
  # token rather than the first: strsplit(" +") would give Ancestry = "union".
  anc_method <- ancestry_all_ROC_data_dataframe %>% pull(Method)
  ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%mutate(
    Method = sub(" +(WB|BB)$", "", anc_method),
    Ancestry = sub("^.* +", "", anc_method)
  )	%>%mutate(Method = fct_relevel(Method,method_levels),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

  p_ROC_ancestry<-fix_legend(ROC_ancestry_fun(ancestry_all_ROC_data_dataframe)+
                               scale_color_manual(values=method_colors, labels = availability_method_labels))

  ###################
  #
  #FDR&Power
  #
  ###################

  FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
  anc_method <- FDR_Power_ancestry %>% pull(Method)
  FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
    Method = sub(" +(WB|BB)$", "", anc_method),
    Ancestry = sub("^.* +", "", anc_method)
  )%>%mutate(Method = fct_relevel(Method,method_levels),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))


  p_FDR_Power_ancestry<-fix_legend(FDR_Power_ancestry_fun(FDR_Power_ancestry%>%filter(FDR!=0.01))+
                                     scale_fill_manual(values=method_colors, labels = availability_method_labels)+
                                     ylim(0, min(ceiling(max(FDR_Power_ancestry$Power)*10)/10+0.1,1)))
  ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(2, 1))
  if(make_figures) safe_ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=280, width=340, units = "mm",dpi=500)
} else {
  ROC_FDR_Power_ancestry <- NULL
  message("no identifiable ancestry-specific targets in ", assv_scenario, "; skipping the Ancestry panel")
}









##############################################################################################################################################################




##########################################################
#
#
#             Summarized result
#
#
##########################################################


p_size_box_raw = ggplot(data=all_Set_data_dataframe, aes(x = Method, y =log2(Size+1), fill=Method)) +
  geom_boxplot(aes(fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray") +
  scale_fill_manual(values=method_colors, guide = FALSE, labels = availability_method_labels) +
  ylab("log2(Set Size + 1)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)+ custom_theme()
p_size_box = fix_xaxis(p_size_box_raw)



p_power_bar = ggplot(data=set_power_summary%>%group_by(Method)%>%summarize(Power_name= mean(Power_name)),aes(x = Method, y =Power_name,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+
  scale_fill_manual(values=method_colors,guide=FALSE, labels = availability_method_labels)+
  ylab("Power")+xlab("")


p_power_bar = p_power_bar+geom_text(aes(x=Method,y=Power_name,label=round(Power_name,2)),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)
p_power_bar = p_power_bar +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)+theme_bw()+custom_theme()
# fix_xaxis MUST come after theme_bw()/custom_theme(): a complete theme like
# theme_bw() REPLACES the whole theme, so an axis.text.x setting made before it
# is silently discarded.  This is panel b of summarize.pdf -- the method names
# overlap unless the rotation is applied last.
p_power_bar = fix_xaxis(p_power_bar)



roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either_All, 
                      SuSiE_ROC_Either = SuSiE_ROC_Either_All, 

                      INFUSE_ROC_Either = INFUSE_ROC_Either_All,
                      INFUSE_union_ROC_Either = INFUSE_union_ROC_Either_All,
                      INFUSE_H_ROC_Either = INFUSE_H_ROC_Either_All,
                      MESuSiE_union_ROC_Either = MESuSiE_union_ROC_Either_All,
                      SuSiEx_intersect_ROC_Either = SuSiEx_intersect_ROC_Either_All,
                      SuSiEx_union_ROC_Either = SuSiEx_union_ROC_Either_All,
                      MultiSuSiE_ROC_Either = MultiSuSiE_ROC_Either_All)
# method_names must follow the ORDER of roc_data_list, not the display ordering
method_names <- c("MESuSiE intersect", "SuSiE", "INFUSE intersect","INFUSE union","INFUSE-H","MESuSiE union","SuSiEx intersect","SuSiEx union","MultiSuSiE")
either_all_ROC_data_dataframe_mean<- create_data_frame(roc_data_list, method_names, 600, h2, causal_num)

#	either_all_ROC_data_dataframe_mean<-either_all_ROC_data_dataframe%>%filter(h2=="~h^2 == 2%*%10^-4",causal_num == "Num~Causal  == 3 ")
p_roc_either_raw <- ggplot(either_all_ROC_data_dataframe_mean, aes(x = Power, y = 1 - FDR)) + 
  geom_line(aes(color = Method), size = 0.5) +
  scale_color_manual(values=method_colors, labels = availability_method_labels) + 
  xlab("Recall") + 
  ylab("Precision") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  geom_point(data = either_all_ROC_data_dataframe_mean %>% filter(Cutoff == 0.5), 
             aes(x = Power, y = 1 - FDR, group = Method, color = Method, shape = 'PIP > 0.5'), 
             size = 2.5, show.legend = c(color = FALSE)) +
  scale_shape_manual(name = ' ', values = c('PIP > 0.5' = 21)) + guides(color = guide_legend(order = 1), 
                                                                        fill = guide_legend(order = 1), 
                                                                        shape = guide_legend(order = 2))+custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
p_roc_either <- fix_legend(p_roc_either_raw)






FDR_Power_Combined<-rbind( FDR_Power_either%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Either"),
                           if(have_shared) FDR_Power_shared%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Shared") else NULL,
                           if(have_ancestry) FDR_Power_ancestry%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.05)%>%mutate(Cat = "Ancestry") else NULL)



FDR_Power_Combined$Method <- factor(
  FDR_Power_Combined$Method,
  levels = c(
    "SuSiE",
    "SuSiEx intersect",
    "SuSiEx union",
    "MultiSuSiE",
    "MESuSiE intersect",
    "INFUSE intersect",
    "INFUSE union",
    "INFUSE-H",
    "MESuSiE union"
  )
)
FDR_Power_Combined$Cat<-factor(	 FDR_Power_Combined$Cat,levels = c("Either","Shared","Ancestry"))
#FDR_Power_Combined$Method<-factor(FDR_Power_Combined$Method,levels = c("SuSiE","Paintor","Paintor_fun","MESuSiE intersect","MESuSiE_fun"))
p_FDR_Power_raw<-ggplot(  FDR_Power_Combined, aes(Cat, Power, fill = Method))+
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(x=Cat,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  scale_fill_manual(values=method_colors, labels = availability_method_labels)+
  ylab("Power")+xlab("FDR = 0.05")+ ylim(0,0.6)+
  theme_bw() + custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+theme(legend.position = 'bottom')



FDR_Power_Combined_1_percent<-rbind( FDR_Power_either%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                                       filter(FDR==0.1)%>%mutate(Cat = "Either"),
                                     if(have_shared) FDR_Power_shared%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                                       filter(FDR==0.1)%>%mutate(Cat = "Shared") else NULL,
                                     if(have_ancestry) FDR_Power_ancestry%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.1)%>%mutate(Cat = "Ancestry") else NULL)



FDR_Power_Combined_1_percent$Method <- factor(
  FDR_Power_Combined_1_percent$Method,
  levels = c(
    "SuSiE",
    "SuSiEx intersect",
    "SuSiEx union",
    "MultiSuSiE",
    "MESuSiE intersect",
    "INFUSE intersect",
    "INFUSE union",
    "INFUSE-H",
    "MESuSiE union"
  )
)
FDR_Power_Combined_1_percent$Cat<-factor(	 FDR_Power_Combined_1_percent$Cat,levels = c("Either","Shared","Ancestry"))
#FDR_Power_Combined$Method<-factor(FDR_Power_Combined$Method,levels = c("SuSiE","Paintor","Paintor_fun","MESuSiE intersect","MESuSiE_fun"))
p_FDR_Power_1_percent_raw<-ggplot(  FDR_Power_Combined_1_percent, aes(Cat, Power, fill = Method))+
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(x=Cat,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  scale_fill_manual(values=method_colors, labels = availability_method_labels)+
  ylab("Power")+xlab("FDR = 0.10")+ ylim(0,0.6)+
  theme_bw() + custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+theme(legend.position = 'bottom')
p_FDR_Power_1_percent<-fix_legend(p_FDR_Power_1_percent_raw)
p_FDR_Power<-fix_legend(p_FDR_Power_raw)


size_power<-p_size_box/p_power_bar
roc_fdr_power<-p_roc_either/p_FDR_Power+theme(legend.position = 'bottom')
p_total <- (size_power | roc_fdr_power) + 
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides = "collect") & 
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

p_total
if(make_figures) save(p_size_box,p_power_bar,p_roc_either,p_FDR_Power,p_FDR_Power_1_percent,file = paste0(plot_dir,"summarize_low_dim_assv.Rdata"))
if(make_figures) safe_ggsave(paste0(plot_dir,"summarize.pdf"),p_total,height=220, width=400, units = "mm",dpi=500)

p_total <- (size_power | roc_fdr_power)  +  
  plot_layout(guides = "collect") & 
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

p_total
if(make_figures) safe_ggsave(paste0(plot_dir,"summarize_notag.pdf"),p_total,height=220, width=400, units = "mm",dpi=500)
