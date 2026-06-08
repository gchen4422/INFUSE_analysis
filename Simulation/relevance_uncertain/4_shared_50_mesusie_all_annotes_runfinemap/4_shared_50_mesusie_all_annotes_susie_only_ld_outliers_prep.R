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

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_",num_causal,"/")
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

#annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")  # we take the full set candidate of annotations
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


zfile_sub = zfile[-joint_flagged,]


z_file_out <- paste0(
  data_dir_flipped,
  "CAUSAL_", num_causal,
  "_LOCI_",  LD_BLOCK,
  "_h2_",    h2_num,
  "_flipped_",flip_prop,
  "_rho_",   rho_num,
  "_extracted_mesusie.txt"
)

write.table(zfile_sub,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  



