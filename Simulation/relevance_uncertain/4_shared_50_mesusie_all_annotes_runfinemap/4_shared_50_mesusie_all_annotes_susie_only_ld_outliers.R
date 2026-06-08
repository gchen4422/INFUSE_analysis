library(Rcpp)
library(dplyr)
library(data.table)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/modified_kriging_rss.R")


args <- as.numeric(commandArgs(TRUE))

LD_BLOCK <- args[1]
h2_num <- args[2]
num_causal <- args[3]
flip_prop <- 0.01
rho_num <- args[4]

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
data_dir_flipped = paste0(data_dir,"flipped_",flip_prop,"_rho_",rho_num,"/")
result_dir<-paste0(wrk_dir,"result/w2_modified_llk_ld_outliers_rho_",rho_num,"/")
system(paste0("mkdir -p ",result_dir))


geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
LD_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")

geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
LD_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")


#baselineLF UKB2.2 annoations
#annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/annot/")  # we take the full set candidate of annotations

annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")  # we take the full set candidate of annotations
#Deepsea annoations
#annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/deepsea/annotation_combined")
##############################################
#
#		Format MESuSiE Input
#
##############################################
zfile<-read.table(paste0(
  data_dir_flipped,
  "CAUSAL_", num_causal,
  "_LOCI_",  LD_BLOCK,
  "_h2_",    h2_num,
  "_flipped_",flip_prop,
  "_rho_",   rho_num,
  ".txt"
),header=T)

summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )
summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 

EU_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID
R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)


joint_df <- kriging_rss_joint(
  summary_stat_1$Z, summary_stat_2$Z,
  EU_cov, BB_cov,
  median(summary_stat_1$N), median(summary_stat_2$N)
)


idx_flag_joint = ifelse(joint_df$plot$plot_env$idx_flag_joint > nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint - nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint)
joint_flagged <- sort(unique(c(
  joint_df$plot$plot_env$idx_flag1,
  joint_df$plot$plot_env$idx_flag2,
  idx_flag_joint)))


keep_idx <- seq_along(summary_stat_1$Z)


if (length(joint_flagged) > 0) {
  # drop the flagged indices
  keep_idx <- keep_idx[-joint_flagged]
  
  # subset your summary‐stat tables and LD matrices
  summary_stat_1_sub <- summary_stat_1[keep_idx, ]
  summary_stat_2_sub <- summary_stat_2[keep_idx, ]
  EU_cov_sub        <- EU_cov[keep_idx, keep_idx]
  BB_cov_sub        <- BB_cov[keep_idx, keep_idx]
  
  # build your “_sub” lists
  summary_stat_sd_list_sub <- list(
    EUR = summary_stat_1_sub,
    AFR = summary_stat_2_sub
  )
  R_mat_list_sub <- list(
    EUR = EU_cov_sub,
    AFR = BB_cov_sub
  )
} else {
  # nothing flagged → just reuse the originals
  summary_stat_sd_list_sub <- summary_stat_sd_list
  R_mat_list_sub         <- R_mat_list
}


##############################################
#
#		Run MESuSiE
#
##############################################

library(MESuSiE)
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))
start.time<-Sys.time()
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10)
MESuSiE.time<-Sys.time()-start.time
MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")


# 
# #annot_file_subset<-annot_file_subset%>%select(c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me1_Trynka_common","Coding_UCSC_common"))
# MESuSiE_inf_res<-MESuSiE_inf(R_mat_list,summary_stat_sd_list,estimate_residual_variance =T,estimate_tau =T,tau = c(0,0),max_iter=10,annot = as.matrix(annot_file_subset))


library(susieR)
start.time<-Sys.time()
susie_EU<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F)
susie_BB<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F)
susie.time<-Sys.time()-start.time
susie.time<-as.numeric(susie.time,units = "mins")

MESuSiE_name<-paste0(result_dir,"Res_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".RData")


library(Rcpp)
library(nloptr)
library(glmnet)
sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6.cpp")
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6.R")
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_utility.R")
#sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_working.cpp")
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_working.R")
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_utility_working.R")
#sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_working_bfgs.cpp")
#sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_susie.cpp")
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_susie.R")
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_pca.R")
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_pca_bic.R")

#load("/scratch/negishi/chen4422/hapnest/Simulation_update/deepsea/annotation_combined/selected_10_annotes.Rdata")
#load("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/select_annot.Rdata")

#annot_name_updated = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common","Coding_UCSC_common","Coding_UCSC.flanking.500_common","CTCF_Hoffman_common","Enhancer_Andersson_common","DHS_peaks_Trynka_common")
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))
annot_file_subset<-fread(paste0(annot_dir,"/Region_",LD_BLOCK,".annot")) #%>% select(all_of(annot_name_updated))
annot_file_subset_flipped = annot_file_subset[keep_idx,]
#annot_file_subset$Coding = 1
#colnames(annot_file) = gsub("[^a-zA-Z0-9_]", ".", colnames(annot_file))
#annot_file_subset = annot_file%>% select(all_of(candidate))

# MESuSiE_res_glmnet_annot_subset<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                               ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                               annot = as.matrix(annot_file_subset), 
#                                               annot_method = "glmnet")

# MESuSiE_res_glmnet_annot_subset_optim<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                                     ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                                     annot = as.matrix(annot_file_subset), 
#                                                     annot_method = "glmnet",est_annot_prior = "optim")


MESuSiE_res_mlk <- meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
                                           ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
                                           annot = as.matrix(annot_file_subset), 
                                           annot_method = "mlk")


MESuSiE_res_mlk_annot_subset<-meSuSie_core(R_mat_list_sub,summary_stat_sd_list_sub,L=10,residual_variance=NULL,prior_weights=NULL,
                                           ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
                                           annot = as.matrix(annot_file_subset_flipped), 
                                           annot_method = "mlk")

# MESuSiE_res_mlk_annot_subset_optim<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                                  ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                                  annot = as.matrix(annot_file_subset), 
#                                                  annot_method = "mlk",est_annot_prior = "optim")



save(susie_EU,susie_BB,MESuSiE_res,MESuSiE_res_mlk,
     MESuSiE_res_mlk_annot_subset,
     file = MESuSiE_name )




