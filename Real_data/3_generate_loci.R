library(data.table)
library(dplyr)
library(snpStats)
library(MESuSiE)



preprocess_data <- function(data_matrix) {
  # Replace missing values with column means and scale columns
  data_matrix <- apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x)
    return(x)
  })
  
  return(data_matrix)
}

# Function to compute covariance matrix and write it to a file
compute_cov <- function(data_matrix, file_name) {
  # Compute the covariance matrix
  cov_matrix <- cov2cor(crossprod(data_matrix))
  
  # Ensure the covariance matrix is symmetric
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  
  # Write the covariance matrix to a file
  fwrite(cov_matrix, file_name,sep =" ")
}


# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

#afr_files = afr_files[-2] # ANX does not have risk loci
#eur_files = eur_files[-2] # ANX does not have risk loci

#afr_files = afr_files[-c(1:6)] 
#eur_files = eur_files[-c(1:6)] 

for (afr_file in afr_files[c(1,4)]) {

trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))
colnames(candidate_region)[1] = "CHR"
colnames(candidate_region)[2] = "MinBP"
colnames(candidate_region)[3] = "MaxBP"
eur_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
afr_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
chr = 0

  for(i in 1:nrow(candidate_region)){
    
    start = candidate_region[i,]$MinBP
    end = candidate_region[i,]$MaxBP
    loci_first = eur_ref[(eur_ref$V4 >= start) & (eur_ref$V4 <= end)]
    loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
    loci_snp = loci$V2
    snplist = paste("loci",i, sep = "_")
    dir.create(paste0("risk_loci_ld_eur_", trait))
    write.table(loci_snp,file = paste0(paste0("risk_loci_ld_eur_", trait,"/"),snplist), quote = FALSE, row.names=FALSE, col.names = F)
    ldname = paste0(paste0("risk_loci_ld_eur_", trait,"/loci_"), i)
    command = paste("~/gwas_software/plink2.0/plink2 --bfile ",paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")," --extract",paste0(paste0("risk_loci_ld_eur_", trait,"/"),snplist)," --allow-no-sex --make-bed --out",ldname)
    t1 <- system(command, intern = TRUE)
    
    eur_plink <- read.plink(paste0(ldname,".bed"))
    eur_plink_geno <- as(eur_plink$genotypes, "numeric")
    unique_value_cols_eur <- apply(eur_plink_geno, 2, function(col) length(unique(col)) == 1)
    unique_value_col_indices_eur <- which(unique_value_cols_eur)
    eur_mat = as.matrix(eur_plink_geno)
    eur_mat = preprocess_data(eur_mat)
    compute_cov(eur_mat, paste0(ldname,".ld"))

    loci_first = afr_ref[(afr_ref$V4 >= start) & (afr_ref$V4 <= end)]
    loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
    loci_snp = loci$V2
    dir.create(paste0("risk_loci_ld_afr_", trait))
    write.table(loci_snp,file = paste0(paste0("risk_loci_ld_afr_", trait,"/"),snplist), quote = FALSE, row.names=FALSE, col.names = F)
    ldname = paste0(paste0("risk_loci_ld_afr_", trait,"/loci_"), i)
    command = paste("~/gwas_software/plink2.0/plink2 --bfile ",paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")," --extract",paste0(paste0("risk_loci_ld_afr_", trait,"/"),snplist)," --allow-no-sex --make-bed --out",ldname)
    t1 <- system(command, intern = TRUE)
    
    afr_plink <- read.plink(paste0(ldname,".bed"))
    afr_plink_geno <- as(afr_plink$genotypes, "numeric")
    unique_value_cols_afr <- apply(afr_plink_geno, 2, function(col) length(unique(col)) == 1)
    unique_value_col_indices_afr <- which(unique_value_cols_afr)
    afr_mat = as.matrix(afr_plink_geno)
    afr_mat = preprocess_data(afr_mat)
    compute_cov(afr_mat, paste0(ldname,".ld"))
    

    snp_to_exclude = c(names(unique_value_col_indices_eur),names(unique_value_col_indices_afr))
    write.table(snp_to_exclude,file = paste0(paste0("risk_loci_ld_afr_", trait,"/"),snplist,"_snp_to_exclude"), quote = FALSE, row.names=FALSE, col.names = F)
    
    
    if(chr != loci$V1[1]){
      annotation_file_name<-paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/baselineLF2.2.UKB.",loci$V1[1],".annot.gz")
      annotation_file<-fread(annotation_file_name)
    }
    
    chr = loci$V1[1]
    #loci_first = annotation_file[(annotation_file$BP >= start) & (annotation_file$BP <= end)]
    #loci = loci_first[loci_first$CHR == candidate_region[i,]$CHR,]
    snplist_extract = fread(paste0(paste0("risk_loci_ld_eur_", trait,"/"),snplist),header = F)
    loci_sub = annotation_file %>% filter(SNP %in% snplist_extract$V1) 
    #Remove potential duplication, i.e. snp in same location but with different rsid
    loci_sub = loci_sub[!duplicated(loci_sub$BP), ]
    dir.create(paste0("risk_loci_annot_", trait))
    write.table(loci_sub,file = paste0(paste0("risk_loci_annot_", trait,"/"),snplist,".annot"), quote = FALSE, row.names=FALSE, col.names = T)
    
  }

message(paste("Processed files for trait:", trait))

}

