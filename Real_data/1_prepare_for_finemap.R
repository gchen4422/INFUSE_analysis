library(data.table)
library(tidyverse)
source("utility.R")

ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_afr_maf_unrelated_biallelic_nodup_extracted.bim")
colnames(ref) = c("CHR","SNP","CM","POS","ALT","REF")
ref <- GWAS_QC(ref, 0.01)  # Do maf qc, remove strand Ambiguous Variants, Remove Multi-allelic Variants, Remove Indel and MHC region SNPs


test = fread("../eur_w_ld_chr/w_hm3.snplist")
ref_hm3 = ref %>% filter(SNP %in% test$SNP) %>% select(-3)

fwrite(ref_hm3,file = "../eur_w_ld_chr/hm3_1kg38_snplist",
        quote = FALSE, row.names = FALSE, col.names = TRUE,sep = " ")


# fwrite(
#   ref[,c(2)], 
#   file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_snplist", 
#   quote = FALSE, row.names = FALSE, col.names = FALSE
# )


annot = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/annot_combine/first_three_columns.annot")

# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

# Ensure corresponding AFR/EUR files are matched
for (afr_file in afr_files[2]) {
  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  #trait <- sub("*._afr_.*", "", afr_file)
  trait <- sub("_.*", "", afr_file)
  trait = sub("^\\.\\./", "", trait)
  # Find the matching EUR file
  eur_file <- eur_files[grep(trait, eur_files)]
  
  if (length(eur_file) == 0) {
    message(paste("No EUR file found for AFR file:", afr_file))
    next
  }
  
  # Read AFR summary stats
  afrdata <- fread(afr_file)
  colnames(afrdata)[which(colnames(afrdata) == "BP")] <- "POS"
  colnames(afrdata)[which(colnames(afrdata) == "A1")] <- "ALT"  # This is the effective allele from hail table allele2, see https://saigegit.github.io/SAIGE-doc/docs/single_step2.html
  colnames(afrdata)[which(colnames(afrdata) == "A2")] <- "REF"  # Reference allele
  colnames(afrdata)[which(colnames(afrdata) == "FRQ")] <- "MAF" # With respect to A1
  #afrdata$PVAL = 10^(-1*afrdata$PVAL) #For regenie output give p as -log10p
  
  # Perform QC and filtering
  afrdata_qced <- GWAS_QC(afrdata, 0.01)  # Do maf qc, remove strand Ambiguous Variants, Remove Multi-allelic Variants, Remove Indel and MHC region SNPs
  afrdata_qced_tidy <- afrdata_qced %>% 
    filter(CHR %in% 1:22) %>% 
    #filter(MISSINGRATE <= 0.01) %>% 
    filter(CHR_POS %in% ref$CHR_POS) %>% arrange(CHR,POS) %>% select(-CHISQ)
  
  # Read EUR summary stats
  eurdata <- fread(eur_file)
  colnames(eurdata)[which(colnames(eurdata) == "BP")] <- "POS"
  colnames(eurdata)[which(colnames(eurdata) == "A1")] <- "ALT"
  colnames(eurdata)[which(colnames(eurdata) == "A2")] <- "REF"
  colnames(eurdata)[which(colnames(eurdata) == "FRQ")] <- "MAF" 
  # lambda <- median(eurdata$CHISQ) / qchisq(0.5, df = 1)                #Perform genomic control for eur HDL and TG summary stats
  # eurdata$CHISQ_gc <- eurdata$CHISQ / lambda
  # eurdata$PVAL <- pchisq(eurdata$CHISQ_gc, df = 1, lower.tail = FALSE)                               
  #eurdata$PVAL = 10^(-1*eurdata$PVAL) #For regenie output give p as -log10p
  
  
  
  # Perform QC and filtering
  eurdata_qced <- GWAS_QC(eurdata, 0.01) 
  eurdata_qced_tidy <- eurdata_qced %>% 
    filter(CHR %in% 1:22) %>% 
    #filter(MISSINGRATE <= 0.01) %>% 
    filter(CHR_POS %in% ref$CHR_POS)  %>% arrange(CHR,POS) %>% select(-CHISQ)
  
  # Find common SNPs
  common_snp <- find_common_snps(eurdata_qced_tidy, afrdata_qced_tidy)
  
  # Filter for common SNPs
  eurdata_qced_tidy_tofinemap <- eurdata_qced_tidy %>% filter(CHR_POS %in% common_snp)
  eurdata_qced_tidy_tofinemap$SNP = ref$SNP[match(eurdata_qced_tidy_tofinemap$CHR_POS,ref$CHR_POS)]
  afrdata_qced_tidy_tofinemap <- afrdata_qced_tidy %>% filter(CHR_POS %in% common_snp)
  
  # Perform allele flipping if necessary
  afrdata_qced_tidy_tofinemap <- allele_flip(eurdata_qced_tidy_tofinemap, afrdata_qced_tidy_tofinemap)
  afrdata_qced_tidy_tofinemap$SNP = ref$SNP[match(afrdata_qced_tidy_tofinemap$CHR_POS,ref$CHR_POS)]
  
  
  eurdata_qced_tidy_tofinemap = eurdata_qced_tidy_tofinemap %>% filter(SNP %in% annot$SNP)
  afrdata_qced_tidy_tofinemap = afrdata_qced_tidy_tofinemap %>% filter(SNP %in% annot$SNP)
  
  
  common_snp_ref <- find_common_snps(eurdata_qced_tidy_tofinemap, ref)
  
  eurdata_qced_tidy_tofinemap = eurdata_qced_tidy_tofinemap %>% filter(CHR_POS %in% common_snp_ref)
  afrdata_qced_tidy_tofinemap = afrdata_qced_tidy_tofinemap %>% filter(CHR_POS %in% common_snp_ref)
  
  # Write results
  fwrite(
    eurdata_qced_tidy_tofinemap, 
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"), 
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  
  fwrite(
    afrdata_qced_tidy_tofinemap, 
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"), 
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  
  fwrite(
    eurdata_qced_tidy_tofinemap[,c(1,6)], 
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_GRCh38_combined_format_noMHC_to_finemap_ref_flip"), 
    quote = FALSE, row.names = FALSE, col.names = FALSE,sep = " "
  )
  
  fwrite(
    eurdata_qced_tidy_tofinemap[,c(1)], 
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_GRCh38_combined_format_noMHC_to_finemap_snplist"), 
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  
  message(paste("Processed files for trait:", trait))
}



for (afr_file in afr_files) {
  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  # For EUR
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_eur_maf_unrelated_biallelic_nodup_extracted"
  snplist = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_GRCh38_combined_format_noMHC_to_finemap_snplist")
  ref_flip = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_GRCh38_combined_format_noMHC_to_finemap_ref_flip")
  outfile = paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")
  command = paste0("~/gwas_software/plink2.0/plink2 --bfile ",bfile," --extract ",snplist ," --ref-allele ",ref_flip," 2 1 --make-bed --out ",outfile)
  t1 <- system(command, intern = TRUE)
  
  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  # For AFR
  bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_afr_maf_unrelated_biallelic_nodup_extracted"
  outfile = paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")
  command = paste0("~/gwas_software/plink2.0/plink2 --bfile ",bfile," --extract ",snplist ," --ref-allele ",ref_flip," 2 1 --make-bed --out ",outfile)
  t1 <- system(command, intern = TRUE)
  message(paste("Processed files for trait:", trait))
}

