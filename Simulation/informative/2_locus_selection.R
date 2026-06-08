library(data.table)
library(dplyr)
library(snpStats)
library(MESuSiE)


# Read locus info and filter by region
selected.100.region<-read.table("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/simulation_locus_info.txt",header=T)
#candidate_region <- filter(selected.100.region, sim_region == region)
candidate_region = selected.100.region

eur_ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/eursel50k_nofilter_nodup.bim")
afr_ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr.bim")

#Find snps present in both EUR and AFR ancestry
afr_ref = afr_ref %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5) 
eur_ref = eur_ref %>% mutate(CHR_POS = V2,REF = V6,ALT = V5)




#Find multi-allelic snps in African population
multi_allelic = afr_ref$CHR_POS[which(duplicated(afr_ref$CHR_POS) == T)]
eur_biallelic = eur_ref %>% filter(!(V2 %in% multi_allelic))
afr_biallelic = afr_ref %>% filter(!(CHR_POS %in% multi_allelic))

common_SNP_bim = find_common_snps(eur_biallelic, afr_biallelic)

eur_biallelic <- eur_biallelic %>%
  filter(CHR_POS %in% common_SNP_bim)

afr_biallelic <- afr_biallelic %>%
  filter(CHR_POS %in% common_SNP_bim)


eur_ref = eur_biallelic
afr_ref = afr_biallelic



#identify snps to be flipped
flip_indices <- which(eur_ref$V5 == afr_ref$V6 & eur_ref$V6 == afr_ref$V5)

afr_ref$flip = 0
afr_ref$flip[flip_indices] = 1


#write.table(cbind(afr_ref$V2,eur_ref$REF), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/afr_ref_snp_update.txt", quote = FALSE, row.names=FALSE, col.names = F)
#The below equivalent to flip allele in summary statistics when running fine-mapping
#command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/afr_ref_snp_update.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changed_map_afr")
#t1 <- system(command, intern = TRUE)




# Function to preprocess data by replacing missing values with column means 
# and scaling columns to have zero mean and unit variance
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



chr = 0
for(i in 1:nrow(candidate_region)){
  
  start = candidate_region[i,]$MinBP
  end = candidate_region[i,]$MaxBP  
  loci_first = eur_ref[(eur_ref$V4 >= start) & (eur_ref$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  snplist = paste("loci",i, sep = "_")
  write.table(loci_snp,file = paste0("risk_loci_ld_eur/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_eur/loci", i,sep = "_")
  command = paste("~/gwas_software/plink2/plink --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/eursel50k_nofilter_nodup --extract",paste0("risk_loci_ld_eur/",snplist)," --allow-no-sex --make-bed --out",ldname)
  t1 <- system(command, intern = TRUE)
  
  #Check LD and summary stats allele consistency
  info = fread(paste0("risk_loci_ld_eur/loci_",i,".bim"))
  eur_plink <- read.plink(paste0(ldname,".bed"))
  eur_plink_geno <- as(eur_plink$genotypes, "numeric")
  eur_mat = as.matrix(eur_plink_geno)
  eur_mat = preprocess_data(eur_mat)
  
  compute_cov(eur_mat, paste0(ldname,".ld"))
  #save(eur_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  
  loci_first = afr_ref[(afr_ref$V4 >= start) & (afr_ref$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  write.table(loci_snp,file = paste0("risk_loci_ld_afr/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_afr/loci", i,sep = "_")
  command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changed_map_afr --extract",paste0("risk_loci_ld_afr/",snplist)," --allow-no-sex --make-bed --out",ldname)
  #command = paste("~/gwas_software/plink2/plink --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changed_map_afr --extract",paste0("risk_loci_ld_afr/",snplist)," --allow-no-sex --make-bed --out",ldname)
  
  t1 <- system(command, intern = TRUE)
  
  #Check LD and summary stats allele consistency
  info = fread(paste0("risk_loci_ld_afr/loci_",i,".bim"))
  afr_plink <- read.plink(paste0(ldname,".bed"))
  afr_plink_geno <- as(afr_plink$genotypes, "numeric")
  afr_mat = as.matrix(afr_plink_geno)
  afr_mat = preprocess_data(afr_mat)
  
  compute_cov(afr_mat, paste0(ldname,".ld"))
  #save(afr_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  

  if(chr != loci$V1[1]){
  annotation_file_name<-paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/baselineLF2.2.UKB.",loci$V1[1],".annot.gz")
  annotation_file<-fread(annotation_file_name)
  }
  chr = loci$V1[1]
  loci_first = annotation_file[(annotation_file$BP >= start) & (annotation_file$BP <= end)]
  loci = loci_first[loci_first$CHR == candidate_region[i,]$CHR,] %>% mutate(CHR_POS = paste0(CHR,":",BP))
  snplist_extract = fread(paste0("risk_loci_ld_eur/loci_",i),header = F)
  loci_sub = loci %>% filter(CHR_POS %in% snplist_extract$V1) %>% select(-CHR_POS)
  #Remove potential duplication, i.e. snp in same location but with different rsid
  loci_sub = loci_sub[!duplicated(loci_sub$BP), ]
  write.table(loci_sub,file = paste0("risk_loci_annot/",snplist,".annot"), quote = FALSE, row.names=FALSE, col.names = T)
  
}


