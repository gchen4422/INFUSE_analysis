library(data.table)
library(dplyr)
library(rtracklayer)
library(gwascat)
library(IRanges)
library(GenomicRanges)
library(GenomicScores)
library(tidyr)
wrk_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes"
system(paste0("mkdir -p ",wrk_dir))


# Loop over num_causal, LD_BLOCK, and h2_num




annot_name_updated = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common")
                       
annot_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")
                  
annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/Annotation/Summarized_Annotation.txt")

for(num_causal in 1:3){
    for(h2_num in 1:2){
      
      summary_stat_list <- list()
      
      for(LD_BLOCK in 1:100){
      
      
      
      cat("Running simulation for num_causal =", num_causal,
          "LD_BLOCK =", LD_BLOCK,
          "h2_num =", h2_num, "\n")
      
      # Set working directories based on current num_causal
      wrk_dir   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_", num_causal, "/")
      system(paste0("mkdir -p ", wrk_dir))
      data_dir<-paste0(wrk_dir,"summary_data/w2_shared_75/")
      system(paste0("mkdir -p ", data_dir))
      out_dir<- paste0(wrk_dir, "out/")
      system(paste0("mkdir -p ", out_dir))
      
      # Set genotype and LD directories
      geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
      LD_dir_eur   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_", num_causal, "/summary_data/")
      
      geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
      LD_dir_afr   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_", num_causal, "/summary_data/")
      
      # Set annotation directory
      annot_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")
      
      ##############################################
      #    Format MESuSiE Input
      ##############################################
      
      # Read the summary statistic file for this iteration
      zfile <- read.table(paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num), header = TRUE)

      annotation_file_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_",LD_BLOCK,".annot")
      annot_file<-fread(annotation_file_name) %>% select(c(1:3,all_of(annot_name_updated))) %>% mutate(RSID = paste0(CHR,":",BP)) %>% select(-c(1:3))
      
      summary_stat_combined <- full_join(
        zfile, 
        annot_file, 
        by = c("RSID")      )  %>% mutate(Region = LD_BLOCK, h2 = h2_num)
      
      summary_stat_list[[LD_BLOCK]] <- summary_stat_combined
      
      
      # Read LD matrices
      #EU_cov <- as.matrix(fread(paste0(LD_dir_afr, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD1")))
      #colnames(EU_cov) <- zfile$RSID
      #BB_cov <- as.matrix(fread(paste0(LD_dir_afr, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD2")))
      #colnames(BB_cov) <- zfile$RSID
      #R_mat_list <- list("EUR" = EU_cov, "AFR" = BB_cov)
      } 
      summary_stat_all <- bind_rows(summary_stat_list)
      fwrite(summary_stat_all,paste0(out_dir,"CAUSAL_",num_causal,"_h2_",h2_num,"_REF_ALL"),col.names = T,row.names = F,quote=F)
    } 
} 

rm(list = ls())









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
      local_subset$Posterior_Prob = local_subset$Posterior_Prob/sum(local_subset $Posterior_Prob)
      row_index<-find_index_until_threshold(local_subset $Posterior_Prob,0.95)
      Paintor_res$Paintor_cs[Paintor_res$RSID%in%local_subset$RSID[row_index]]<-1
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



for(num_causal in 1:3){
  for(h2_num in 1:2){
      
    wrk_dir   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_", num_causal, "/")
    system(paste0("mkdir -p ", wrk_dir))
    data_dir  <- paste0(wrk_dir, "summary_data/w2_startover/")
    system(paste0("mkdir -p ", data_dir))
    result_dir<- paste0(wrk_dir, "result/w2_modified_llk/")
    system(paste0("mkdir -p ", result_dir))
    out_dir<- paste0(wrk_dir, "out/")
    system(paste0("mkdir -p ", out_dir))
    
  for(LD_BLOCK in 1:100){
    EUR_COV<-as.matrix(fread(paste0(data_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".LD1")))
    AFR_COV<-as.matrix(fread(paste0(data_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".LD2")))
    max_COV<-pmax(abs(EUR_COV),abs(AFR_COV))
    # Paintor_res_sub_baseline<-fread(paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor"))
    # Paintor_res_sub_baseline_update<-Paintor_credible_set(Paintor_res_sub_baseline,max_COV )
    # fwrite(Paintor_res_sub_baseline_update,paste0(out_res_dir,"loci_",i,".sub_baseline.fun.mcmc.paintor.cs"))
    
    # Paintor_res_sub_20<-fread(paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor"))
    # Paintor_res_sub_20_update<-Paintor_credible_set(Paintor_res_sub_20,max_COV )
    # fwrite(Paintor_res_sub_20_update,paste0(out_res_dir,"loci_",i,".sub_20.fun.mcmc.paintor.cs"))
    
    Paintor_res_all<-fread(paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".fun.mcmc.paintor"))
    Paintor_res_all_update<-Paintor_credible_set(Paintor_res_all,max_COV )
    fwrite(Paintor_res_all_update,paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".fun.mcmc.paintor.cs"))
    
    Paintor_res<-fread(paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".mcmc.paintor"))
    Paintor_res_update<-Paintor_credible_set(Paintor_res,max_COV )
    fwrite(Paintor_res_update,paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".mcmc.paintor.cs"))
    }
  }
}








##################################################################################################
#
#		Step 2: Get the all result to be analyzed
#
##################################################################################################
annotation_file<-fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_3/out/CAUSAL_3_h2_1_REF_ALL") %>% select(c(3,9:15))

for(num_causal in 1:3){
  for(h2_num in 1:2){
    
  all_trait_locus_summary<-c()
  region_id = 0
  wrk_dir   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_", num_causal, "/")
  system(paste0("mkdir -p ", wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/w2_startover/")
  system(paste0("mkdir -p ", data_dir))
  result_dir<-paste0(wrk_dir,"result/w2_modified_llk/")
  system(paste0("mkdir -p ", result_dir))
  out_dir<- paste0(wrk_dir, "out/sparse_pca_modified_llk_updated/")
  system(paste0("mkdir -p ", out_dir))
  
  
  
  ###Check all results are available
  #locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor.cs")%in%list.files(path=result_dir))
  #locus_selected_MESuSiE<-which(paste0("loci_",all_gene_region_trait,"_res.RData")%in%list.files(path=result_dir))
  #locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  #cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected))
  
  all_locus_summary<-c()
  for(LD_BLOCK in 1:100){
    
    region_id = region_id + 1 
    load(paste0(result_dir,"Res_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData"))
    
    locus_summary <- read.table(paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num), header = TRUE)
    

    # if(any(WB$SNP!=paintor_res$SNP)==TRUE){
    #   cat(trait_name,",",region)
    # }
    
    #For PAINTOR result without annotation
    paintor_res<-fread(paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".mcmc.paintor.cs"),header=T)
    paintor_res<-paintor_res%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    #locus_summary<-WB%>%dplyr::rename(N_1 = N,Beta_WB = BETA,zscore_WB = Z, MAF_WB = MAF)%>%select(-PVAL,-SE,-MISSINGRATE)%>%mutate(Region = region_id,MAF_BB = BB$MAF, N_2=BB$N,zscore_BB = BB$Z,Beta_BB = BB$BETA, PHENONAME = trait_name) %>% dplyr::rename(Trait = PHENONAME)
    locus_summary<-locus_summary%>%mutate(Paintor_cs = paintor_res$Paintor_cs, Paintor_PIP = paintor_res$Paintor_PIP)
    
    #For PAINTOR result with all annotations
    paintor_res_all<-fread(paste0(result_dir,paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),".fun.mcmc.paintor.cs"),header=T)
    paintor_res_all<-paintor_res_all%>%dplyr::rename(Paintor_PIP =Posterior_Prob )
    locus_summary<-locus_summary%>%mutate(Paintor_cs_all = paintor_res_all$Paintor_cs, Paintor_PIP_all = paintor_res_all$Paintor_PIP)
    
    #For MESuSiE result without annotation
    locus_summary$MESuSiE_cs<-0
    if(length(unlist(MESuSiE_res$cs$cs))>0)locus_summary$MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either = MESuSiE_res$pip, MESuSiE_PIP_Shared = MESuSiE_res$pip_config[,3], MESuSiE_PIP_WB = MESuSiE_res$pip_config[,1],MESuSiE_PIP_BB = MESuSiE_res$pip_config[,2])
    
    #For MESuSiE result with all baseline annotation use glm
    locus_summary$MESuSiE_res_glmnet_annot_subset_cs<-0
    if(length(unlist(MESuSiE_res_glmnet_annot_subset$cs$cs))>0)locus_summary$MESuSiE_res_glmnet_annot_subset_cs[unlist(MESuSiE_res_glmnet_annot_subset$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_glmnet = MESuSiE_res_glmnet_annot_subset$pip, MESuSiE_PIP_Shared_all_glmnet = MESuSiE_res_glmnet_annot_subset$pip_config[,3], MESuSiE_PIP_WB_all_glmnet = MESuSiE_res_glmnet_annot_subset$pip_config[,1],MESuSiE_PIP_BB_all_glmnet = MESuSiE_res_glmnet_annot_subset$pip_config[,2])
    
    #For MESuSiE result with all baseline annotation use mlk
    locus_summary$MESuSiE_res_mlk_annot_subset_cs<-0
    if(length(unlist(MESuSiE_res_mlk_annot_subset$cs$cs))>0)locus_summary$MESuSiE_res_mlk_annot_subset_cs[unlist(MESuSiE_res_mlk_annot_subset$cs$cs)]<-1
    locus_summary<-locus_summary%>%mutate(MESuSiE_PIP_Either_all_mlk = MESuSiE_res_mlk_annot_subset$pip, MESuSiE_PIP_Shared_all_mlk = MESuSiE_res_mlk_annot_subset$pip_config[,3], MESuSiE_PIP_WB_all_mlk = MESuSiE_res_mlk_annot_subset$pip_config[,1],MESuSiE_PIP_BB_all_mlk = MESuSiE_res_mlk_annot_subset$pip_config[,2])
    
    
    #For SuSiE result without annotation
    locus_summary$SuSiE_cs<-0
    if(length(unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs))))>0)locus_summary$SuSiE_cs[unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs)))]<-1
    locus_summary<-locus_summary%>%mutate(SuSiE_PIP = pmax(susie_EU$pip,susie_BB$pip), SuSiE_Shared = pmin(susie_EU$pip,susie_BB$pip),susie_EU = susie_EU$pip,susie_BB = susie_BB$pip)
    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(LD_BLOCK)
    
  }
  ###Adding functional annotations to the data
  annotation_file_order<-annotation_file[match(all_locus_summary$RSID,annotation_file$RSID),]
  all_locus_summary<-cbind(all_locus_summary,annotation_file_order)
  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)
  
  all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
  res_all<-all_trait_locus_summary
  
  res_name<-paste0(out_dir,"res_simulation_h2_",h2_num,"_pca_bic.RData")
  save(res_all,file = res_name)
  
  }
}



  
