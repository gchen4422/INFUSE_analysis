library(data.table)
library(MESuSiE)
library(susieR)
library(snpStats)
library(dplyr)

# args <- as.numeric(commandArgs(TRUE))
# i <- args[1]
# h2_num <- args[2]
# num_causal <- args[3]


preprocess_data <- function(data_matrix) {
  # Replace missing values with column means and scale columns
  data_matrix <- apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x)
    return(x)
  })
  
  return(data_matrix)
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
  # eur_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/processed/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  # afr_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/processed/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  # chr = 0
  
  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/risk_loci_ld_eur_",trait,"/")
  LD_dir_eur = geno_dir_eur
  
  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/risk_loci_ld_afr_",trait,"/")
  LD_dir_afr = geno_dir_afr
  
  annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/risk_loci_annot_",trait,"/")
  
  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))
  
  
  #for(i in 1:20){
  for(i in 1:nrow(candidate_region)){
  #for(i in 1:5){
    
    snplist_ori = fread(paste0(geno_dir_afr,"loci_",i),header = F)%>% pull(V1)
    snp_to_exclude_path <- paste0(geno_dir_afr, "loci_", i, "_snp_to_exclude")
    
    # Check if the file is empty
    if (file.info(snp_to_exclude_path)$size > 0) {
      # Read the file if it's not empty
      snp_to_exclude <- fread(snp_to_exclude_path, header = FALSE) %>% pull(V1)
      snplist_filtered <- setdiff(snplist_ori, snp_to_exclude)
    } else {
      # Handle the case where the file is empty
      snplist_filtered  = snplist_ori
    }
    
    
    ##Check for LD mismatch, remove mismatched SNPs
    idx_to_remove = c()
    while(!(identical(idx_to_remove, integer(0)))){
      
      
      
      sumstats_afr = afr_sumstats_sig %>% filter(SNP %in% snplist_filtered)
      sumstats_eur = eur_sumstats_sig %>% filter(SNP %in% snplist_filtered)

      
      eur_plink <- read.plink(paste0(geno_dir_eur,"loci_",i,".bed"))
      eur_plink_geno <- as(eur_plink$genotypes, "numeric")
      unique_value_cols_eur <- apply(eur_plink_geno, 2, function(col) length(unique(col)) == 1)
      unique_value_col_indices_eur <- which(unique_value_cols_eur)
      eur_mat = as.matrix(eur_plink_geno)
      eur_mat <- eur_mat[, colnames(eur_mat) %in% snplist_filtered]
      eur_mat = preprocess_data(eur_mat)
      cov_matrix_eur <- cov2cor(crossprod(eur_mat))
      EU_cov <- (cov_matrix_eur + t(cov_matrix_eur)) / 2

      
      afr_plink <- read.plink(paste0(geno_dir_afr,"loci_",i,".bed"))
      afr_plink_geno <- as(afr_plink$genotypes, "numeric")
      unique_value_cols_afr <- apply(afr_plink_geno, 2, function(col) length(unique(col)) == 1)
      unique_value_col_indices_afr <- which(unique_value_cols_afr)
      afr_mat = as.matrix(afr_plink_geno)
      afr_mat <- afr_mat[, colnames(afr_mat) %in% snplist_filtered]
      afr_mat = preprocess_data(afr_mat)
      cov_matrix_afr <- cov2cor(crossprod(afr_mat))
      BB_cov <- (cov_matrix_afr + t(cov_matrix_afr)) / 2
      
      
      EU_diagnostic <- kriging_rss(sumstats_eur$Z, EU_cov, n = median(sumstats_eur$N))
      idx_to_remove_eur = EU_diagnostic$plot$plot_env$idx
      BB_diagnostic <- kriging_rss(sumstats_afr$Z, BB_cov, n = median(sumstats_eur$N))
      idx_to_remove_afr = BB_diagnostic$plot$plot_env$idx
      idx_to_remove = sort(union(idx_to_remove_eur,idx_to_remove_afr))
      
      if (!(identical(idx_to_remove, integer(0)))){
        snplist_filtered = snplist_filtered[-idx_to_remove]
      }
      
      gc()
      
    }
    
    

    fwrite(sumstats_afr,file = paste0(data_dir,"afr_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    fwrite(sumstats_eur,file = paste0(data_dir,"eur_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    
    fwrite(BB_cov, paste0(data_dir,"afr_loci_",i,".ld"),sep =" ")
    fwrite(EU_cov, paste0(data_dir,"eur_loci_",i,".ld"),sep =" ")
    
    annot_list = fread(paste0(annot_dir,"loci_",i,".annot"))
    annot_filtered = annot_list %>% filter(SNP %in% snplist_filtered)
    write.table(annot_filtered,file = paste0(data_dir,"loci_",i,".annot"), quote = FALSE, row.names=FALSE, col.names = T)
  }
  
  message(paste("Processed files for trait:", trait))
  
}