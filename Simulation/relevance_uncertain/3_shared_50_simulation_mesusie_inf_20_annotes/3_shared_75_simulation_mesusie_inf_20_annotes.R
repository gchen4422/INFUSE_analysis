# Loading necessary libraries
library(data.table)
library(snpStats)
library(mvtnorm)
library(dplyr)

# Setting constants and initializing parameters
#args <- as.numeric(commandArgs(TRUE))
#LD_BLOCK <- args[1]
num_EU <- 300000
num_BB <- 300000

geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
LD_dir_eur <- geno_dir_eur

geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
LD_dir_afr <- geno_dir_afr

annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/annot/")


# Function to process genotype data
process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}



for (LD_BLOCK in 1:100){
  
  
  
  # Seed initialization
  set.seed(1128)
  seed_LOCIrator <- sample.int(1000, 100)[LD_BLOCK]
  set.seed(seed_LOCIrator)
  
  
  
  
  # Reading and processing EU data
  outfile_EU <- paste0(geno_dir_eur, "/loci_", LD_BLOCK)
  EU_plink <- read.plink(paste0(outfile_EU, ".bed"))
  #order_CHR_POS <- EU_plink$map %>% mutate(CHR_POS = paste0(chromosome, "_", position)) %>% pull(CHR_POS)
  EU_plink_geno <- as(EU_plink$genotypes, "numeric")
  EU_plink_geno <- process_geno(EU_plink_geno)
  
  ld_EU_name <- paste0(LD_dir_eur, "/loci_", LD_BLOCK, ".ld")
  EU_cov <- as.matrix(fread(ld_EU_name))
  
  # Reading and processing BB data
  outfile_BB <- paste0(geno_dir_afr, "/loci_", LD_BLOCK)
  BB_plink <- read.plink(paste0(outfile_BB, ".bed"))
  BB_plink$map$CHR_POS = paste0(BB_plink$map$chromosome,":",BB_plink$map$position)
  BB_plink_geno <- as(BB_plink$genotypes, "numeric")
  colnames(BB_plink_geno) = BB_plink$map$CHR_POS
  BB_plink_geno <- process_geno(BB_plink_geno)
  
  ld_BB_name <- paste0(LD_dir_afr, "/loci_", LD_BLOCK, ".ld")
  BB_cov<-as.matrix(fread(ld_BB_name))
  colnames(BB_cov) = BB_plink$map$CHR_POS
  
  # Check SNP position mismatch 
  if(any(colnames(EU_cov)!=EU_plink$map$snp.name)|any(colnames(BB_cov)!=BB_plink$map$snp.name)){
    cat("SNP mismatch of geno and LD reference")
  }
  # Checking allele match
  if (any(EU_plink$map$allele.1 != BB_plink$map$allele.1)) {
    cat("Allele mismatch of the two ancestries")
  }
  
  
  
  # Add functional annotation information to simulation
  load("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/select_annot.Rdata")
  annotation_file_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_",LD_BLOCK,".annot")
  annot_file<-fread(annotation_file_name) 
  # annot_file<-annot_file%>%filter(SNP%in%EU_plink$map$snp.name)%>%arrange(match(SNP,EU_plink$map$snp.name))
  # annot_file_subset<-annot_file%>%select(c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common","Conserved_Vertebrate_phastCons46way.flanking.500_common","Conserved_Primate_phastCons46way_common","MAF_Adj_Predicted_Allele_Age_common","Repressed_Hoffman_common",
  #                                          "MAF_Adj_LLD_AFR_common","DGF_ENCODE_common","DGF_ENCODE.flanking.500_common","H3K4me3_Trynka.flanking.500_common","WeakEnhancer_Hoffman.flanking.500_common","MAF_Adj_ASMC_common","Nucleotide_Diversity_10kb_common","GTEx_eQTL_MaxCPP_common", "BLUEPRINT_DNA_methylation_MaxCPP_common","H3K9ac_Trynka_common","CpG_Content_50kb_common"))
  
  annot_file_subset<-annot_file%>%select(all_of(annot_name_updated))
  enrichment_beta<-c(rep(2,5),rep(0,15))
  
  #annot_file_subset<-annot_file%>%select(c("Promoter_UCSC_lowfreq","synonymous_common","H3K9ac_Trynka_common","SuperEnhancer_Hnisz_lowfreq","UTR_3_UCSC_lowfreq"))
  #enrichment_beta<-rep(2,5)
  
  true_prior_prob<-exp(as.matrix(as.matrix(annot_file_subset)%*%enrichment_beta))/sum(exp(as.matrix(annot_file_subset)%*%enrichment_beta))
  true_prior_prob = drop(true_prior_prob)
  save(true_prior_prob, file = paste0("true_prior/weight_",LD_BLOCK,".Rdata"))
  
  write.table(annot_file_subset,paste0(annot_dir,"Region_",LD_BLOCK,".annot"),quote=F,row.names=F,col.names=T,sep=" ")
  
  # Non-zero for infinitetimal effect 
  h2_small_EU = 0/nrow(EU_cov)
  h2_small_BB = 0/nrow(BB_cov)
  
  
  for(num_causal in 1:3){
    
    
    wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_",num_causal,"/")
    system(paste0("mkdir -p ",wrk_dir))
    data_dir<-paste0(wrk_dir,"summary_data/w2_shared_75/")
    result_dir<-paste0(wrk_dir,"result/")
    out_dir<-paste0(wrk_dir,"out/")
    system(paste0("mkdir -p ",data_dir))
    system(paste0("mkdir -p ",result_dir))
    system(paste0("mkdir -p ",out_dir))
    
    ##sample number of ancestry-specific causal SNP, and then sample from SNP list
    ##Better way can be sample causal SNP list first, and sample 0,1 to define as EU/BB specific
    #num_shared_causal_SNP<-c(1,3,5)[num_causal]
    #num_unique_causal_SNP<-c(1,3,5)[num_causal]
    num_causal_SNP <- c(2,3,5)[num_causal]
    causal_SNP_config <- sample(c(1,2,3),num_causal_SNP,replace = TRUE,prob = c(0.125,0.125,0.75))
    causal_SNP_list<-BB_plink$map$snp.name[sample(seq(1,nrow(BB_cov)),num_causal_SNP,replace=F,true_prior_prob)]
    
    print(LD_BLOCK)
    for(h2_num in 1:2){
      
      h2= c(1e-4,2e-4)[h2_num] 
      signal<-beta_EUR_all<-beta_AFR_all<-rep(0,nrow(EU_cov))  
      
      for(causal_index in 1:num_causal_SNP){
        
        if(causal_SNP_config[causal_index]==1){
          
          EU_plink_causal<-matrix(EU_plink_geno[,colnames(EU_plink_geno)%in%causal_SNP_list[causal_index]],ncol=1)
          beta_EUR_causal<-rnorm(1,0,1)
          beta_EUR_causal_scale<-beta_EUR_causal*sqrt(h2/drop(var(EU_plink_causal*beta_EUR_causal)))
          beta_EUR_all[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-beta_EUR_causal_scale
          
          signal[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-1
        }else if(causal_SNP_config[causal_index]==2){
          
          BB_plink_causal<-matrix(BB_plink_geno[,colnames(BB_plink_geno)%in%causal_SNP_list[causal_index]],ncol=1)
          beta_AFR_causal<-rnorm(1,0,1)
          beta_AFR_causal_scale<-beta_AFR_causal*sqrt(h2/drop(var(BB_plink_causal*beta_AFR_causal)))
          beta_AFR_all[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-beta_AFR_causal_scale
          
          signal[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-2
          
        }else if(causal_SNP_config[causal_index]==3){
          
          EU_plink_causal<-matrix(EU_plink_geno[,colnames(EU_plink_geno)%in%causal_SNP_list[causal_index]],ncol=1)
          BB_plink_causal<-matrix(BB_plink_geno[,colnames(BB_plink_geno)%in%causal_SNP_list[causal_index]],ncol=1)
          beta_shared<-rmvnorm(1,mean=rep(0,2),sigma=matrix(c(1,0.8,0.8,1),ncol=2,nrow=2))
          beta_EUR_causal_scale<-beta_shared[1]*sqrt(h2/drop(var(EU_plink_causal*beta_shared[1])))
          beta_AFR_causal_scale<-beta_shared[2]*sqrt(h2/drop(var(BB_plink_causal*beta_shared[2])))
          beta_EUR_all[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-beta_EUR_causal_scale
          beta_AFR_all[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-beta_AFR_causal_scale
          
          signal[which(EU_plink$map$snp.name%in%causal_SNP_list[causal_index])]<-3
        }
      }
      beta_EUR_small<-rnorm(nrow(EU_cov),0,h2_small_EU)
      beta_EUR_marginal<-as.matrix(EU_cov)%*%(beta_EUR_all+beta_EUR_small)
      
      beta_AFR_small<-rnorm(nrow(BB_cov),0,h2_small_BB)
      beta_AFR_marginal<-as.matrix(BB_cov)%*%(beta_AFR_all + beta_AFR_small)
      
      h2_causal_EUR<-drop(var(EU_plink_geno%*%beta_EUR_all))
      h2_causal_AFR<-drop(var(BB_plink_geno%*%beta_AFR_all))
      y_null_EU<-rnorm(nrow(EU_plink_geno),0,sqrt(1-h2_causal_EUR-h2_small_EU*nrow(EU_cov)))
      y_null_EU<-y_null_EU-mean(y_null_EU)
      y_null_BB<-rnorm(nrow(BB_plink_geno),0,sqrt(1-h2_causal_AFR-h2_small_BB*nrow(EU_cov)))
      y_null_BB<-y_null_BB-mean(y_null_BB)
      err_beta_EU<-t(EU_plink_geno)%*%y_null_EU/nrow(EU_plink_geno)
      err_beta_BB<-t(BB_plink_geno)%*%y_null_BB/nrow(BB_plink_geno)
      
      err_beta_EU_scale<-sqrt(nrow(EU_plink_geno)/num_EU)*err_beta_EU
      err_beta_BB_scale<-sqrt(nrow(BB_plink_geno)/num_BB)*err_beta_BB
      
      z_EU<-(beta_EUR_marginal+err_beta_EU_scale)*sqrt(num_EU)
      z_BB<-(beta_AFR_marginal+err_beta_BB_scale)*sqrt(num_BB)
      
      z_file<-data.frame("CHR" = as.numeric(EU_plink$map$chr[1]),"POS" =EU_plink$map$position ,"RSID" =EU_plink$map$snp.name ,"zscore_1"=z_EU,"zscore_2"=z_BB,"Signal"=signal,"N_1" = num_EU,"N_2" = num_BB)
      
      z_file_out<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
      write.table(z_file,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
      
      #ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
      #ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
      #fwrite(EU_cov,ld_EU_paintor_name,sep=" ",col.names=F)
      #fwrite(BB_cov,ld_BB_paintor_name,sep=" ",col.names=F)
    }
  }
  
}
