library(data.table)
library(tidyverse)
library(GenomicRanges)



# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]

for (afr_file in afr_files) {
  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  # For EUR
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv")) %>% filter(PVAL < 5e-8)
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv")) %>% filter(PVAL < 5e-8)
  
  significant_snps = rbind(afr_sumstats_sig[,1:3],eur_sumstats_sig[,1:3])
  significant_snps = significant_snps %>% arrange(CHR,POS)
  
  significant_snps <- significant_snps %>%
    mutate(
      start = POS - 500000,
      end = POS + 500000
    )
  
  # Convert significant SNPs to a GRanges object
  granges_snps <- GRanges(
    seqnames = as.factor(significant_snps$CHR),
    ranges = IRanges(start = significant_snps$start, end = significant_snps$end),
    SNP = significant_snps$SNP
  )
  
  # Merge overlapping regions
  merged_ranges <-  GenomicRanges::reduce(granges_snps)
  
  # Convert the merged genomic ranges back to a data frame
  merged_ranges_df <- as.data.frame(merged_ranges)
  
  fwrite(
    merged_ranges_df, 
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/", trait, "_risk_loci"), 
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  message(paste("Processed files for trait:", trait))
}
