library(Rcpp)
library(dplyr)
args <- as.numeric(commandArgs(TRUE))

LD_BLOCK <- args[1]
h2_num <- args[2]
num_causal <- args[3]
flip_prop <- 0.01
rho_num <- args[4]



library(data.table)

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


annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")

#load("select_annot_20.Rdata")
##############################################
#
#		Format  Input
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
#annot_name_updated = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common","Coding_UCSC_common","Coding_UCSC.flanking.500_common","CTCF_Hoffman_common","Enhancer_Andersson_common","DHS_peaks_Trynka_common")
annot_file_subset<-fread(paste0(annot_dir,"Region_",LD_BLOCK,".annot"))#%>% select(all_of(annot_name_updated))
##############################################
#
#		Run Paintor
#
##############################################
annotation_paintor_name<-paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".txt.annotations")
write.table(annot_file_subset,annotation_paintor_name,col.names = T,row.names = F,quote=F, sep=" ")

input_file_name<-paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".input")
write.table(paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".txt"),input_file_name,col.names = F,row.names = F,quote=F)


paintor_dir<-data_dir_flipped
paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num)


ld_EU_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".txt.LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num,".txt.LD2")
fwrite(EU_cov,ld_EU_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
fwrite(BB_cov,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")

paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_flipped_",flip_prop, "_rho_",rho_num)
input_dir<-paintor_dir

setwd(input_dir)





annotations = colnames(annot_file_subset)
# Concatenate annotations into a single string
annotations_part <- paste0("-annotations ", paste(annotations, collapse = ","))

# Full command with annotations part separated
 system(paste0(
   "~/gwas_software/PAINTOR_V3.0/PAINTOR -input -input ", input_file_name,
   " -Zhead zscore_1,zscore_2 -LDname LD1,LD2 -in ", input_dir,
   " -out ", result_dir, " -mcmc ", annotations_part,
   " -Gname ", paste0(paintor_suffix, ".enrich"),
   " -Lname ", paste0(paintor_suffix, "_bayes_factor"),
   " -RESname ", paste0("fun.mcmc.paintor")
 ))



system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))

#Remove LD files
#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))
