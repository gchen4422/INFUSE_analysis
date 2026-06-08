library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)


###Further form Paintor credible set 

find_index_until_threshold <- function(vec, threshold) {
  # Sort the input vector in decreasing order and obtain the sorted indices
  sorted_indices <- order(vec, decreasing = TRUE)
  
  # Calculate the cumulative sum of the sorted vector
  cumsum_vec <- cumsum(vec[sorted_indices])
  
  # Find the index where the cumulative sum reaches or exceeds the threshold
  index <- which(cumsum_vec >= threshold)[1]
  
  # Return the original indices of the elements up to the found index
  sorted_indices[1:index]
}

Paintor_Partiton<-function(Paintor_res,COV_mat,threshold = 0.5){
  
  Paintor_subset<-Paintor_res
  COV_mat_subset<-COV_mat
  while(sum(Paintor_subset$Posterior_Prob)>0.95){
    
    local_index = which.max(Paintor_subset$Posterior_Prob)
    cor_SNP_index = which(COV_mat_subset[,local_index]>threshold)
    if(length(cor_SNP_index)!=0&sum(Paintor_subset$Posterior_Prob[cor_SNP_index])>0.95){
      
      local_subset = Paintor_subset[cor_SNP_index,]
      local_subset $Posterior_Prob = local_subset $Posterior_Prob/sum(local_subset $Posterior_Prob)
      row_index<-find_index_until_threshold(local_subset $Posterior_Prob,0.95)
      Paintor_res$Paintor_cs[Paintor_res$SNP%in%local_subset$SNP[row_index]]<-1
    }
    COV_mat_subset<-COV_mat_subset[-cor_SNP_index,-cor_SNP_index]
    Paintor_subset<-Paintor_subset[-cor_SNP_index,]
    
  }
  return(Paintor_res)
}

Paintor_credible_set<-function(Paintor_res,COV_mat){
  Paintor_res$Paintor_cs = 0 	
  if(sum(Paintor_res$Posterior_Prob)<1){
    Paintor_res$Paintor_cs = 0 	
  }else{
    
    Paintor_res<-Paintor_Partiton(Paintor_res,COV_mat,threshold = 0.5)
    
  }
  
  return(Paintor_res)
}


# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]



for (afr_file in afr_files) {
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  out_dir =paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"_ld_combined/summary_data/"),"paintor_all/")
  out_res_dir = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"_ld_combined/result/")
  
  for(i in 1:nrow(candidate_region)){
    EUR_COV<-as.matrix(fread(paste0(out_dir,"loci_",i,".LD1")))
    AFR_COV<-as.matrix(fread(paste0(out_dir,"loci_",i,".LD2")))
    max_COV<-pmax(abs(EUR_COV),abs(AFR_COV))
    # Paintor_res_sub_baseline<-fread(paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor"))
    # Paintor_res_sub_baseline_update<-Paintor_credible_set(Paintor_res_sub_baseline,max_COV )
    # fwrite(Paintor_res_sub_baseline_update,paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor.cs"))
     
    # Paintor_res_sub_20<-fread(paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor"))
    # Paintor_res_sub_20_update<-Paintor_credible_set(Paintor_res_sub_20,max_COV )
    # fwrite(Paintor_res_sub_20_update,paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor.cs"))
    
    Paintor_res_all<-fread(paste0(out_res_dir,"loci_",i,".all.fun.mcmc.paintor"))
    Paintor_res_all_update<-Paintor_credible_set(Paintor_res_all,max_COV )
    fwrite(Paintor_res_all_update,paste0(out_res_dir,"loci_",i,".all.fun.mcmc.paintor.cs"))
    
    Paintor_res<-fread(paste0(out_res_dir,"loci_",i,".mcmc.paintor"))
    Paintor_res_update<-Paintor_credible_set(Paintor_res,max_COV )
    fwrite(Paintor_res_update,paste0(out_res_dir,"loci_",i,".mcmc.paintor.cs"))
  }
  message(paste("Processed files for trait:", trait))
}









##################################################################################################
#
#		Step 2: Get the all result to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]



#annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/Summarized_Annotation_BP_updated_cadd.txt")
#table(annotation_file$brain_ind_eQTL)
all_trait_locus_summary<-c()
region_id = 0
for (afr_file in afr_files) {
  trait_name <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait_name, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  region_length_file<-candidate_region%>%summarise(min_count = min(width),max_count = max(width))
  region_length_file$Trait=trait_name
  cat(trait_name)
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"_ld_combined/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"_ld_combined/"),"summary_data/")

  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE<-which(paste0("loci_",all_gene_region_trait,"_res.RData")%in%list.files(path=result_dir))
  locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected))
  
  all_locus_summary<-c()
  for(region in locus_selected){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res.RData"))
    
    WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    if(trait_name == "LDL" |trait_name == "TC" ){
    WB = WB %>% select(-CHISQ)
    }

    
    #For PAINTOR result without annotation
    paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    

    #For PAINTOR result with all annotations
    paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    #For MESuSiE result without annotation
    locus_summary$MESuSiE_cs<-0
    if(length(unlist(MESuSiE_res$cs$cs))>0)locus_summary$MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either = MESuSiE_res$pip, MESuSiE_PIP_Shared = MESuSiE_res$pip_config[,3], MESuSiE_PIP_WB = MESuSiE_res$pip_config[,1],MESuSiE_PIP_BB = MESuSiE_res$pip_config[,2])
    

    #For MESuSiE result with all baseline annotation use glm
    locus_summary$MESuSiE_res_glmnet_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_glmnet_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_glmnet_annot_all_cs[unlist(MESuSiE_res_glmnet_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_glmnet = MESuSiE_res_glmnet_annot_all$pip, MESuSiE_PIP_Shared_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,2])
    
    #For MESuSiE result with all baseline annotation use mlk
    locus_summary$MESuSiE_res_mlk_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_mlk_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_mlk_annot_all_cs[unlist(MESuSiE_res_mlk_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_mlk = MESuSiE_res_mlk_annot_all$pip, MESuSiE_PIP_Shared_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,2])
    
    
    #For SuSiE result without annotation
    locus_summary$SuSiE_cs<-0
    if(length(unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs))))>0)locus_summary$SuSiE_cs[unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs)))]<-1
    locus_summary<-locus_summary%>%mutate(SuSiE_PIP = pmax(susie_EU$pip,susie_BB$pip), SuSiE_Shared = pmin(susie_EU$pip,susie_BB$pip),susie_EU = susie_EU$pip,susie_BB = susie_BB$pip)
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
    
  }
  ###Adding functional annotations to the data
  #annotation_file_order<-annotation_file[match(all_locus_summary$CHR_POS,annotation_file$SNP),]
  #all_locus_summary<-cbind(all_locus_summary,annotation_file_order)
  
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
  
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_sbp_pca_bic_modified_llk_tocheck_ld.RData")
save(res_all,file = res_name)







##################################For results without LD outlier correction###############################


for (afr_file in afr_files) {
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  out_dir =paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"_ld_nocontrol/summary_data/"),"paintor_all/")
  out_res_dir = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"_ld_nocontrol/result/")
  
  for(i in 1:nrow(candidate_region)){
    EUR_COV<-as.matrix(fread(paste0(out_dir,"loci_",i,".LD1")))
    AFR_COV<-as.matrix(fread(paste0(out_dir,"loci_",i,".LD2")))
    max_COV<-pmax(abs(EUR_COV),abs(AFR_COV))
    # Paintor_res_sub_baseline<-fread(paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor"))
    # Paintor_res_sub_baseline_update<-Paintor_credible_set(Paintor_res_sub_baseline,max_COV )
    # fwrite(Paintor_res_sub_baseline_update,paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor.cs"))
    
    # Paintor_res_sub_20<-fread(paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor"))
    # Paintor_res_sub_20_update<-Paintor_credible_set(Paintor_res_sub_20,max_COV )
    # fwrite(Paintor_res_sub_20_update,paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor.cs"))
    
    Paintor_res_all<-fread(paste0(out_res_dir,"loci_",i,".all.fun.mcmc.paintor"))
    Paintor_res_all_update<-Paintor_credible_set(Paintor_res_all,max_COV )
    fwrite(Paintor_res_all_update,paste0(out_res_dir,"loci_",i,".all.fun.mcmc.paintor.cs"))
    
    Paintor_res<-fread(paste0(out_res_dir,"loci_",i,".mcmc.paintor"))
    Paintor_res_update<-Paintor_credible_set(Paintor_res,max_COV )
    fwrite(Paintor_res_update,paste0(out_res_dir,"loci_",i,".mcmc.paintor.cs"))
  }
  message(paste("Processed files for trait:", trait))
}









##################################################################################################
#
#		Step 2: Get the all result to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]



#annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/Summarized_Annotation_BP_updated_cadd.txt")
#table(annotation_file$brain_ind_eQTL)
all_trait_locus_summary<-c()
region_id = 0
for (afr_file in afr_files) {
  trait_name <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait_name, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  region_length_file<-candidate_region%>%summarise(min_count = min(width),max_count = max(width))
  region_length_file$Trait=trait_name
  cat(trait_name)
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"_ld_nocontrol/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"_ld_nocontrol/"),"summary_data/")
  
  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE<-which(paste0("loci_",all_gene_region_trait,"_res.RData")%in%list.files(path=result_dir))
  locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected))
  
  all_locus_summary<-c()
  for(region in locus_selected){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res.RData"))
    
    WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    if(trait_name == "LDL" |trait_name == "TC" ){
      WB = WB %>% select(-CHISQ)
    }
    
    
    #For PAINTOR result without annotation
    paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    
    #For PAINTOR result with all annotations
    paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    #For MESuSiE result without annotation
    locus_summary$MESuSiE_cs<-0
    if(length(unlist(MESuSiE_res$cs$cs))>0)locus_summary$MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either = MESuSiE_res$pip, MESuSiE_PIP_Shared = MESuSiE_res$pip_config[,3], MESuSiE_PIP_WB = MESuSiE_res$pip_config[,1],MESuSiE_PIP_BB = MESuSiE_res$pip_config[,2])
    
    
    #For MESuSiE result with all baseline annotation use glm
    locus_summary$MESuSiE_res_glmnet_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_glmnet_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_glmnet_annot_all_cs[unlist(MESuSiE_res_glmnet_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_glmnet = MESuSiE_res_glmnet_annot_all$pip, MESuSiE_PIP_Shared_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_glmnet = MESuSiE_res_glmnet_annot_all$pip_config[,2])
    
    #For MESuSiE result with all baseline annotation use mlk
    locus_summary$MESuSiE_res_mlk_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_mlk_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_mlk_annot_all_cs[unlist(MESuSiE_res_mlk_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_mlk = MESuSiE_res_mlk_annot_all$pip, MESuSiE_PIP_Shared_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,2])
    
    
    #For SuSiE result without annotation
    locus_summary$SuSiE_cs<-0
    if(length(unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs))))>0)locus_summary$SuSiE_cs[unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs)))]<-1
    locus_summary<-locus_summary%>%mutate(SuSiE_PIP = pmax(susie_EU$pip,susie_BB$pip), SuSiE_Shared = pmin(susie_EU$pip,susie_BB$pip),susie_EU = susie_EU$pip,susie_BB = susie_BB$pip)
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
    
  }
  ###Adding functional annotations to the data
  #annotation_file_order<-annotation_file[match(all_locus_summary$CHR_POS,annotation_file$SNP),]
  #all_locus_summary<-cbind(all_locus_summary,annotation_file_order)
  
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
  
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_sbp_pca_bic_modified_llk_nocontrol_ld.RData")
save(res_all,file = res_name)




