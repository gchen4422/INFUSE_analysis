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
  out_dir =paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/summary_data/"),"paintor_all/")
  out_res_dir = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/result/")
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
#		Step 2: Get all BP result to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData")) ### use this as most recent BP results




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
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"summary_data/")
  
  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  #locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE_evo2<-which(paste0("loci_",all_gene_region_trait,"_res_evo2.RData")%in%list.files(path=result_dir))
  #locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected_MESuSiE_evo2))
  
  
  all_locus_summary<-c()
  for(region in locus_selected_MESuSiE_evo2){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res_evo2.RData"))
    
    locus_summary = res_all %>% filter(Trait == trait_name ) %>% filter(Region == region_id)
    evo2_delta_score= fread(paste0(sumstat_dir,"loci_",region,"_summ_scored.csv")) %>% select(SNP,evo2_delta_score)
    locus_summary <- locus_summary %>% left_join(evo2_delta_score, by = "SNP")   
    if(any(is.na(locus_summary$evo2_delta_score))) {
      warning("⚠️ Some SNPs in locus_summary did not find a matching Evo2 score!")
    }
    
    #WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    #BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    
    
    #For PAINTOR result without annotation
    #paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    #paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    #locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    
    #For PAINTOR result with all annotations
    #paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    #paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    

    # --- 1. Process MESuSiE: Transformed Evo2 ONLY (Model B) ---
    locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip,
      PIP_Shared_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip_config[,3],
      PIP_EUR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,1], # WB/EUR
      PIP_AFR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,2]  # BB/AFR
    )
    
    # --- 2. Process MESuSiE: Classic + Transformed Evo2 (Model C) ---
    locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Ensemble_infuse= MESuSiE_res_mlk_annot_all_log_evo2$pip,
      PIP_Shared_Ensemble_infuse = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,3],
      PIP_EUR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,1],
      PIP_AFR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,2]
    )
  
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
  }
  ###Adding functional annotations to the data
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all_new<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_bp_pca_bic_modified_llk_updated_evo2.RData")
save(res_all_new,file = res_name)



##################################################################################################
#
#		Step 2: Get all BMI result to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250318_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData")) ### use this as most recent BP results




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
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"summary_data/")
  
  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  #locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE_evo2<-which(paste0("loci_",all_gene_region_trait,"_res_evo2.RData")%in%list.files(path=result_dir))
  #locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected_MESuSiE_evo2))
  
  
  all_locus_summary<-c()
  for(region in locus_selected_MESuSiE_evo2){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res_evo2.RData"))
    
    locus_summary = res_all %>% filter(Trait == trait_name ) %>% filter(Region == region_id)
    evo2_delta_score= fread(paste0(sumstat_dir,"loci_",region,"_summ_scored.csv")) %>% select(SNP,evo2_delta_score)
    locus_summary <- locus_summary %>% left_join(evo2_delta_score, by = "SNP")   
    if(any(is.na(locus_summary$evo2_delta_score))) {
      warning("⚠️ Some SNPs in locus_summary did not find a matching Evo2 score!")
    }
    
    #WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    #BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    
    
    #For PAINTOR result without annotation
    #paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    #paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    #locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    
    #For PAINTOR result with all annotations
    #paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    #paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    
    # --- 1. Process MESuSiE: Transformed Evo2 ONLY (Model B) ---
    locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip,
      PIP_Shared_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip_config[,3],
      PIP_EUR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,1], # WB/EUR
      PIP_AFR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,2]  # BB/AFR
    )
    
    # --- 2. Process MESuSiE: Classic + Transformed Evo2 (Model C) ---
    locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Ensemble_infuse= MESuSiE_res_mlk_annot_all_log_evo2$pip,
      PIP_Shared_Ensemble_infuse = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,3],
      PIP_EUR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,1],
      PIP_AFR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,2]
    )
    
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
  }
  ###Adding functional annotations to the data
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all_new<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_bmi_pca_bic_modified_llk_evo2.RData")
save(res_all_new,file = res_name)


##################################################################################################
#
#		Step 2: Get all BMI result with thread set to 1 to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250318_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData")) ### use this as most recent BP results




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
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"summary_data/")
  
  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  #locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE_evo2<-which(paste0("loci_",all_gene_region_trait,"_res_evo2_return_test.RData")%in%list.files(path=result_dir))
  #locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected_MESuSiE_evo2))
  
  
  all_locus_summary<-c()
  for(region in locus_selected_MESuSiE_evo2){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res_evo2_return_test.RData"))
    
    locus_summary = res_all %>% filter(Trait == trait_name ) %>% filter(Region == region_id)
    evo2_delta_score= fread(paste0(sumstat_dir,"loci_",region,"_summ_scored.csv")) %>% select(SNP,evo2_delta_score)
    locus_summary <- locus_summary %>% left_join(evo2_delta_score, by = "SNP")   
    if(any(is.na(locus_summary$evo2_delta_score))) {
      warning("⚠️ Some SNPs in locus_summary did not find a matching Evo2 score!")
    }
    
    #WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    #BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    
    
    #For PAINTOR result without annotation
    #paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    #paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    #locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    
    #For PAINTOR result with all annotations
    #paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    #paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    
    # --- 1. Process MESuSiE: Transformed Evo2 ONLY (Model B) ---
    locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip,
      PIP_Shared_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip_config[,3],
      PIP_EUR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,1], # WB/EUR
      PIP_AFR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,2]  # BB/AFR
    )
    
    # --- 2. Process MESuSiE: Classic + Transformed Evo2 (Model C) ---
    locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Ensemble_infuse= MESuSiE_res_mlk_annot_all_log_evo2$pip,
      PIP_Shared_Ensemble_infuse = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,3],
      PIP_EUR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,1],
      PIP_AFR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,2]
    )
    
    
    #For MESuSiE result without annotation
    locus_summary$MESuSiE_cs<-0
    if(length(unlist(MESuSiE_res$cs$cs))>0)locus_summary$MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either = MESuSiE_res$pip, MESuSiE_PIP_Shared = MESuSiE_res$pip_config[,3], MESuSiE_PIP_WB = MESuSiE_res$pip_config[,1],MESuSiE_PIP_BB = MESuSiE_res$pip_config[,2])
    
    
    #For MESuSiE result with all baseline annotation use mlk
    locus_summary$MESuSiE_res_mlk_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_mlk_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_mlk_annot_all_cs[unlist(MESuSiE_res_mlk_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_mlk = MESuSiE_res_mlk_annot_all$pip, MESuSiE_PIP_Shared_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,2])
    
    
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
  }
  ###Adding functional annotations to the data
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all_new<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_bmi_pca_bic_modified_llk_evo2_rerun_test.RData")
save(res_all_new,file = res_name)



##################################################################################################



res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"



load(paste0(res_dir,"res_bmi_pca_bic_modified_llk_evo2.RData"))   ### use this as most recent BMI results
res_all = res_all_new
rm(res_all_new)
#res_all = res_all %>% filter(Trait == "SBP")


load(paste0(res_dir,"res_bmi_pca_bic_modified_llk_evo2_rerun_test.RData"))   ### use this as most recent BMI results

library(dplyr)
library(tibble)

cols_to_compare <- c(
  # Evo2 only
  "MESuSiE_res_mlk_annot_log_evo2_cs",
  "PIP_Either_Evo2_Only_infuse",
  "PIP_Shared_Evo2_Only_infuse",
  "PIP_EUR_Evo2_Only_infuse",
  "PIP_AFR_Evo2_Only_infuse",
  
  # Classic + Evo2
  "MESuSiE_res_mlk_annot_all_log_evo2_cs",
  "PIP_Either_Ensemble_infuse",
  "PIP_Shared_Ensemble_infuse",
  "PIP_EUR_Ensemble_infuse",
  "PIP_AFR_Ensemble_infuse",
  
  # MESuSiE baseline
  "MESuSiE_cs",
  "MESuSiE_PIP_Either",
  "MESuSiE_PIP_Shared",
  "MESuSiE_PIP_WB",
  "MESuSiE_PIP_BB",
  
  # MESuSiE + all classic annotations (MLK)
  "MESuSiE_res_mlk_annot_all_cs",
  "MESuSiE_PIP_Either_all_mlk",
  "MESuSiE_PIP_Shared_all_mlk",
  "MESuSiE_PIP_WB_all_mlk",
  "MESuSiE_PIP_BB_all_mlk"
)

cols_to_compare <- intersect(cols_to_compare, intersect(names(res_all), names(res_all_new)))

compare_one_col <- function(new_vec, old_vec, col_name, tol = 1e-8) {
  stopifnot(length(new_vec) == length(old_vec))
  
  diff_vec <- new_vec - old_vec
  abs_diff <- abs(diff_vec)
  valid <- !is.na(abs_diff)
  
  # avoid NA correlation when one vector is constant
  corr_val <- suppressWarnings(
    if (sd(new_vec, na.rm = TRUE) == 0 || sd(old_vec, na.rm = TRUE) == 0) {
      NA_real_
    } else {
      cor(new_vec, old_vec, use = "complete.obs")
    }
  )
  
  tibble(
    column = col_name,
    n = length(new_vec),
    n_na_new = sum(is.na(new_vec)),
    n_na_old = sum(is.na(old_vec)),
    n_both_non_na = sum(!is.na(new_vec) & !is.na(old_vec)),
    identical_exact = identical(new_vec, old_vec),
    identical_within_tol = all(abs_diff[valid] < tol),
    n_changed = sum(abs_diff[valid] >= tol),
    prop_changed = mean(abs_diff[valid] >= tol),
    mean_diff = mean(diff_vec, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    median_abs_diff = median(abs_diff, na.rm = TRUE),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    correlation = corr_val
  )
}

diff_summary <- bind_rows(lapply(cols_to_compare, function(cc) {
  compare_one_col(res_all_new[[cc]], res_all[[cc]], cc, tol = 1e-8)
}))

print(diff_summary)




col_name <- "MESuSiE_PIP_Either_all_mlk"

pip_diff_tbl <- tibble(
  row_id = res_all$SNP,
  CHR_POS = res_all$CHR_POS,
  old_pip = res_all[[col_name]],
  new_pip = res_all_new[[col_name]],
  diff = res_all_new[[col_name]] - res_all[[col_name]],
  abs_diff = abs(diff),
  old_cs = res_all$MESuSiE_res_mlk_annot_all_cs,
  new_cs = res_all_new$MESuSiE_res_mlk_annot_all_cs
)

# top changed SNPs
top_snp_changes <- pip_diff_tbl %>%
  arrange(desc(abs_diff))

print(top_snp_changes, n = 50)




##################################################################################################
#
#		Step 2: Get all BP result with thread set to 1 to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData")) ### use this as most recent BP results




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
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait_name,"/"),"summary_data/")
  
  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  #locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  locus_selected_MESuSiE_evo2<-which(paste0("loci_",all_gene_region_trait,"_res_evo2_return_test.RData")%in%list.files(path=result_dir))
  #locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected_MESuSiE_evo2))
  
  
  all_locus_summary<-c()
  for(region in locus_selected_MESuSiE_evo2){
    region_id = region_id + 1 
    load(paste0(result_dir,"loci_",region,"_res_evo2_return_test.RData"))
    
    locus_summary = res_all %>% filter(Trait == trait_name ) %>% filter(Region == region_id)
    evo2_delta_score= fread(paste0(sumstat_dir,"loci_",region,"_summ_scored.csv")) %>% select(SNP,evo2_delta_score)
    locus_summary <- locus_summary %>% left_join(evo2_delta_score, by = "SNP")   
    if(any(is.na(locus_summary$evo2_delta_score))) {
      warning("⚠️ Some SNPs in locus_summary did not find a matching Evo2 score!")
    }
    
    #WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    #BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    
    
    
    #For PAINTOR result without annotation
    #paintor_res<-fread(paste0(result_dir,"loci_",region,".mcmc.paintor.cs"),header=T)
    #paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-INFO)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    #locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    
    #For PAINTOR result with all annotations
    #paintor_res_all<-fread(paste0(result_dir,"loci_",region,".all.fun.mcmc.paintor.cs"),header=T)
    #paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    
    # --- 1. Process MESuSiE: Transformed Evo2 ONLY (Model B) ---
    locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip,
      PIP_Shared_Evo2_Only_infuse = MESuSiE_res_mlk_annot_log_evo2$pip_config[,3],
      PIP_EUR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,1], # WB/EUR
      PIP_AFR_Evo2_Only_infuse    = MESuSiE_res_mlk_annot_log_evo2$pip_config[,2]  # BB/AFR
    )
    
    # --- 2. Process MESuSiE: Classic + Transformed Evo2 (Model C) ---
    locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs <- 0
    if(length(unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)) > 0) {
      locus_summary$MESuSiE_res_mlk_annot_all_log_evo2_cs[unlist(MESuSiE_res_mlk_annot_all_log_evo2$cs$cs)] <- 1
    }
    
    locus_summary <- locus_summary %>% mutate(
      PIP_Either_Ensemble_infuse= MESuSiE_res_mlk_annot_all_log_evo2$pip,
      PIP_Shared_Ensemble_infuse = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,3],
      PIP_EUR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,1],
      PIP_AFR_Ensemble_infuse    = MESuSiE_res_mlk_annot_all_log_evo2$pip_config[,2]
    )
    
    #For MESuSiE result without annotation
    locus_summary$MESuSiE_cs<-0
    if(length(unlist(MESuSiE_res$cs$cs))>0)locus_summary$MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either = MESuSiE_res$pip, MESuSiE_PIP_Shared = MESuSiE_res$pip_config[,3], MESuSiE_PIP_WB = MESuSiE_res$pip_config[,1],MESuSiE_PIP_BB = MESuSiE_res$pip_config[,2])
    
    
    #For MESuSiE result with all baseline annotation use mlk
    locus_summary$MESuSiE_res_mlk_annot_all_cs<-0
    if(length(unlist(MESuSiE_res_mlk_annot_all$cs$cs))>0)locus_summary$MESuSiE_res_mlk_annot_all_cs[unlist(MESuSiE_res_mlk_annot_all$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_mlk = MESuSiE_res_mlk_annot_all$pip, MESuSiE_PIP_Shared_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,3], MESuSiE_PIP_WB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,1],MESuSiE_PIP_BB_all_mlk = MESuSiE_res_mlk_annot_all$pip_config[,2])
    
    
    
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)
  }
  ###Adding functional annotations to the data
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
res_all_new<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_bp_pca_bic_modified_llk_updated_evo2_rerun_test.RData")
save(res_all_new,file = res_name)




##################################################################################################



res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"



load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated_evo2.RData"))   ### use this as most recent BP results
res_all = res_all_new
rm(res_all_new)
#res_all = res_all %>% filter(Trait == "SBP")


load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated_evo2_rerun_test.RData"))   ### use this as most recent BP results

library(dplyr)
library(tibble)

cols_to_compare <- c(
  # Evo2 only
  "MESuSiE_res_mlk_annot_log_evo2_cs",
  "PIP_Either_Evo2_Only_infuse",
  "PIP_Shared_Evo2_Only_infuse",
  "PIP_EUR_Evo2_Only_infuse",
  "PIP_AFR_Evo2_Only_infuse",
  
  # Classic + Evo2
  "MESuSiE_res_mlk_annot_all_log_evo2_cs",
  "PIP_Either_Ensemble_infuse",
  "PIP_Shared_Ensemble_infuse",
  "PIP_EUR_Ensemble_infuse",
  "PIP_AFR_Ensemble_infuse",
  
  # MESuSiE baseline
  "MESuSiE_cs",
  "MESuSiE_PIP_Either",
  "MESuSiE_PIP_Shared",
  "MESuSiE_PIP_WB",
  "MESuSiE_PIP_BB",
  
  # MESuSiE + all classic annotations (MLK)
  "MESuSiE_res_mlk_annot_all_cs",
  "MESuSiE_PIP_Either_all_mlk",
  "MESuSiE_PIP_Shared_all_mlk",
  "MESuSiE_PIP_WB_all_mlk",
  "MESuSiE_PIP_BB_all_mlk"
)

cols_to_compare <- intersect(cols_to_compare, intersect(names(res_all), names(res_all_new)))

compare_one_col <- function(new_vec, old_vec, col_name, tol = 1e-8) {
  stopifnot(length(new_vec) == length(old_vec))
  
  diff_vec <- new_vec - old_vec
  abs_diff <- abs(diff_vec)
  valid <- !is.na(abs_diff)
  
  # avoid NA correlation when one vector is constant
  corr_val <- suppressWarnings(
    if (sd(new_vec, na.rm = TRUE) == 0 || sd(old_vec, na.rm = TRUE) == 0) {
      NA_real_
    } else {
      cor(new_vec, old_vec, use = "complete.obs")
    }
  )
  
  tibble(
    column = col_name,
    n = length(new_vec),
    n_na_new = sum(is.na(new_vec)),
    n_na_old = sum(is.na(old_vec)),
    n_both_non_na = sum(!is.na(new_vec) & !is.na(old_vec)),
    identical_exact = identical(new_vec, old_vec),
    identical_within_tol = all(abs_diff[valid] < tol),
    n_changed = sum(abs_diff[valid] >= tol),
    prop_changed = mean(abs_diff[valid] >= tol),
    mean_diff = mean(diff_vec, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    median_abs_diff = median(abs_diff, na.rm = TRUE),
    max_abs_diff = max(abs_diff, na.rm = TRUE),
    correlation = corr_val
  )
}

diff_summary <- bind_rows(lapply(cols_to_compare, function(cc) {
  compare_one_col(res_all_new[[cc]], res_all[[cc]], cc, tol = 1e-10)
}))

print(diff_summary)




col_name <- "MESuSiE_PIP_Either_all_mlk"

pip_diff_tbl <- tibble(
  row_id = res_all$SNP,
  CHR_POS = res_all$CHR_POS,
  old_pip = res_all[[col_name]],
  new_pip = res_all_new[[col_name]],
  diff = res_all_new[[col_name]] - res_all[[col_name]],
  abs_diff = abs(diff),
  old_cs = res_all$MESuSiE_res_mlk_annot_all_cs,
  new_cs = res_all_new$MESuSiE_res_mlk_annot_all_cs
)

# top changed SNPs
top_snp_changes <- pip_diff_tbl %>%
  arrange(desc(abs_diff))

print(top_snp_changes, n = 50)


