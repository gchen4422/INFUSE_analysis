library(Rcpp)
library(dplyr)
library(data.table)
#library(RhpcBLASctl)
#blas_set_num_threads(1)
#omp_set_num_threads(1)


#cat("BLAS threads:", blas_get_num_procs(), "\n")
#cat("OMP threads:", omp_get_max_threads(), "\n")

args <- commandArgs(TRUE)
trait <- "DBP"         # trait as a string
LD_BLOCK <- as.numeric(args[1])  # LD_BLOCK as a numeric value


candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
#system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
result_dir<-paste0(wrk_dir,"result/")
out_dir<-paste0(wrk_dir,"out/")
#system(paste0("mkdir -p ",data_dir))
#system(paste0("mkdir -p ",result_dir))
#system(paste0("mkdir -p ",out_dir))
#for(i in 1:nrow(candidate_region)){

##############################################
#
#		Format MESuSiE Input
#
##############################################
#zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
summary_stat_eur = fread(paste0(data_dir,"eur_loci_",LD_BLOCK,"_summ"))
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "BETA")] <- "Beta"
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "SE")] <- "Se"
summary_stat_afr = fread(paste0(data_dir,"afr_loci_",LD_BLOCK,"_summ"))
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "BETA")] <- "Beta"
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "SE")] <- "Se"
summary_stat_sd_list = list("EUR" = summary_stat_eur,"AFR"=summary_stat_afr ) 

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
#colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
#colnames(BB_cov)<-zfile$RSID
R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)
##############################################
#
#		Run MESuSiE
#
##############################################

library(MESuSiE)
meSuSie_core_ori = meSuSie_core
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))
start.time<-Sys.time()
MESuSiE_res<-meSuSie_core_ori(R_mat_list,summary_stat_sd_list,L=10)
MESuSiE.time<-Sys.time()-start.time
MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")




#library(susieR)
#start.time<-Sys.time()
#susie_EU<-susie_rss(summary_stat_eur$Z,EU_cov,check_prior=F)
#susie_BB<-susie_rss(summary_stat_afr$Z,BB_cov,check_prior=F)
#susie.time<-Sys.time()-start.time
#susie.time<-as.numeric(susie.time,units = "mins")


library(Rcpp)
library(nloptr)
library(glmnet)
sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6.cpp")
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_pca_bic_updated_evo2.R")
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_utility.R")
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))



### Use all existing annotations

annot_file_all<-fread(paste0(data_dir,"loci_",LD_BLOCK,".annot")) %>% select(-c(1:4)) %>% select(-GTEx_eQTL_MaxCPP_common)
evo2_delta_score= fread(paste0(data_dir,"loci_",LD_BLOCK,"_summ_scored.csv")) %>% select(c(1:3),evo2_delta_score)

# Check if they all have the exact same number of rows
cat("--- Row Counts ---\n")
row_counts <- c(
  EUR   = nrow(summary_stat_eur),
  AFR   = nrow(summary_stat_afr),
  Evo2  = nrow(evo2_delta_score),
  Annot = nrow(annot_file_all)
)
print(row_counts)

# Check if the SNP order matches perfectly (using EUR as the baseline)
Evo2_matches_EUR <- identical(summary_stat_eur$SNP, evo2_delta_score$SNP)

if (Evo2_matches_EUR) {
  
  evo2_vector <- abs(evo2_delta_score$evo2_delta_score)
  log_evo2 <- log10(evo2_vector + 1e-9)
  transformed_evo2 <- (log_evo2 - min(log_evo2)) / (max(log_evo2) - min(log_evo2))
  
  annot_file_all_evo2 <- cbind(annot_file_all, abs_evo2 = evo2_vector)
  annot_file_all_log_evo2 <- cbind(annot_file_all, transformed_evo2 = transformed_evo2)
  
  
  message("✅ Success: 'annot_file_all_evo' and 'evo2_vector' created using absolute Evo2 scores.")
} else {
  stop("🚨 ERROR: SNP order mismatch detected between Evo2 and EUR summary stats. Pipeline halted to prevent data corruption.")
}

# ---------------------------------------------------------
# MODEL A: Classic Annotations Only
# ---------------------------------------------------------
MESuSiE_res_mlk_annot_all <- infuse_core(R_mat_list, summary_stat_sd_list, L=10, residual_variance=NULL, prior_weights=NULL,
                                         ancestry_weight=ancestry_weight, optim_method ="optim", estimate_residual_variance=F, max_iter=100,
                                         annot = as.matrix(annot_file_all), 
                                         annot_method = "mlk")

# ---------------------------------------------------------
# MODEL B: Transformed Evo2 Only (Zero-Shot)
# ---------------------------------------------------------
MESuSiE_res_mlk_annot_log_evo2 <- infuse_core(R_mat_list, summary_stat_sd_list, L=10, residual_variance=NULL, prior_weights=NULL,
                                              ancestry_weight=ancestry_weight, optim_method ="optim", estimate_residual_variance=F, max_iter=100,
                                              annot = as.matrix(data.frame(Transformed_Evo2 = transformed_evo2)), 
                                              annot_method = "mlk")

# ---------------------------------------------------------
# MODEL C: Classic Annotations + Transformed Evo2 (Ensemble)
# ---------------------------------------------------------
MESuSiE_res_mlk_annot_all_log_evo2 <- infuse_core(R_mat_list, summary_stat_sd_list, L=10, residual_variance=NULL, prior_weights=NULL,
                                                  ancestry_weight=ancestry_weight, optim_method ="optim", estimate_residual_variance=F, max_iter=100,
                                                  annot = as.matrix(annot_file_all_log_evo2), 
                                                  annot_method = "mlk")

# ---------------------------------------------------------
# SAVE THE CORRECT RESULTS
# ---------------------------------------------------------
MESuSiE_name <- paste0(result_dir, "loci_", LD_BLOCK, "_res_evo2_return_test.RData")

save(MESuSiE_res,MESuSiE_res_mlk_annot_all,MESuSiE_res_mlk_annot_log_evo2, 
     MESuSiE_res_mlk_annot_all_log_evo2,
     file = MESuSiE_name)




#start.time<-Sys.time()
#MESuSiE_res_mlk_annot_all<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                        ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                        annot = as.matrix(annot_file_all), 
#                                        annot_method = "mlk")
#MESuSiE_res_mlk_annot_all.time<-Sys.time()-start.time
#MESuSiE_res_mlk_annot_all.time<-as.numeric(MESuSiE_res_mlk_annot_all.time,units = "mins")


#MESuSiE_res_mlk_annot_evo2<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                             ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                             annot = as.matrix(evo2_vector), 
#                                             annot_method = "mlk")

#MESuSiE_res_mlk_annot_log_evo2<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                         ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                         annot = as.matrix(transformed_evo2), 
#                                         annot_method = "mlk")



#MESuSiE_res_mlk_annot_all<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                        ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                        annot = as.matrix(annot_file_all), 
#                                        annot_method = "mlk")

#MESuSiE_res_mlk_annot_all_evo2<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                                 ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                                 annot = as.matrix(annot_file_all_evo2), 
#                                                 annot_method = "mlk")

#MESuSiE_res_mlk_annot_all_log_evo2<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                             ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                             annot = as.matrix(annot_file_all_log_evo2), 
#                                             annot_method = "mlk")


#MESuSiE_res_glmnet_annot_all<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
#                                           ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
#                                           annot = as.matrix(annot_file_all), 
#                                           annot_method = "glmnet")
#MESuSiE_res_glmnet_annot_all.time<-Sys.time()-start.time
#MESuSiE_res_glmnet_annot_all.time<-as.numeric(MESuSiE_res_glmnet_annot_all.time,units = "mins")




# MESuSiE_res_mlk_annot$pip[zfile$Signal!=0]
# MESuSiE_res_mlk$pip[zfile$Signal!=0]
# MESuSiE_res_mlk_annot$pip[zfile$Signal!=0]

##############################################
#
#		Run Paintor
#
##############################################


#For baseline all annotations
# system(paste0("mkdir -p ",paste0(data_dir,"paintor_all/")))
# paintor_dir<-paste0(data_dir,"paintor_all/")
# annotation_paintor_name<-paste0(paintor_dir,"loci_",LD_BLOCK,".annotations")
# write.table(annot_file_all,annotation_paintor_name,col.names = T,row.names = F,quote=F, sep=" ")
# 
# input_file_name<-paste0(paintor_dir,"loci_",LD_BLOCK,".input")
# write.table(paste0("loci_",LD_BLOCK),input_file_name,col.names = F,row.names = F,quote=F)
# #paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
# sumstats = cbind(summary_stat_eur[,c(1,8)],summary_stat_afr[,c(8)])
# colnames(sumstats) = c("SNP","zscore_1","zscore_2")
# fwrite(sumstats,paste0(paintor_dir,"loci_",LD_BLOCK),col.names = T,row.names = F,quote=F,sep=" ")
# 
# 
# # Concatenate annotations into a single string
# annotations = colnames(annot_file_all)
# annotations_part <- paste0("-annotations ", paste(annotations, collapse = ","))
# 
# 
# ld_EU_paintor_name<-paste0(paintor_dir,"loci_",LD_BLOCK,".LD1")
# ld_BB_paintor_name<-paste0(paintor_dir,"loci_",LD_BLOCK,".LD2")
# fwrite(EU_cov,ld_EU_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
# fwrite(BB_cov,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
# 
# #paintor_suffix<-paste0(result_dir,"loci_",LD_BLOCK,"_res")
# paintor_suffix<-paste0("loci_",LD_BLOCK)
# input_dir<-paintor_dir
# 
# setwd(input_dir)
# 
# # Full command with annotations part separated
# start.time<-Sys.time()
# system(paste0(
#   "~/gwas_software/PAINTOR_V3.0/PAINTOR -input -input ", input_file_name,
#   " -Zhead zscore_1,zscore_2 -LDname LD1,LD2 -in ", input_dir,
#   " -out ", result_dir, " -mcmc ", annotations_part,
#   " -Gname ", paste0(paintor_suffix, "_all.enrich"),
#   " -Lname ", paste0(paintor_suffix, "_all_bayes_factor"),
#   " -RESname ", paste0("all.fun.mcmc.paintor")
# ))
# Paintor_all.time<-Sys.time()-start.time
# Paintor_all.time<-as.numeric(Paintor_all.time,units = "mins")
# 
# start.time<-Sys.time()
# system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))
# Paintor.time<-Sys.time()-start.time
# Paintor.time<-as.numeric(Paintor.time,units = "mins")

#Remove LD files
#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))

##############################################
#
#		Save time used
#
##############################################

# time_out_name<-paste0(result_dir,"time_loci_",LD_BLOCK,".txt")
# time_vec<-c(MESuSiE.time,susie.time,
#                 MESuSiE_res_glmnet_annot_all.time, MESuSiE_res_mlk_annot_all.time,
#                 Paintor_all.time, Paintor.time)
# write.table(time_vec,time_out_name,col.names = F,row.names = F,quote=F)

