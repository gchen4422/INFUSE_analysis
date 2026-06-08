library(data.table)
library(dplyr)
library(rtracklayer)
library(gwascat)
library(IRanges)
library(GenomicRanges)
library(GenomicScores)
library(tidyr)
wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/"
system(paste0("mkdir -p ",wrk_dir))


# List all files
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250417_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]



for (afr_file in afr_files){
  
  
  
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))
  data_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/"),"summary_data/")
  out_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/"),"out/")
  
  summary_stat_list_eur = list()
  summary_stat_list_afr = list()
  summary_stat_list <- list()
  
  for (i in 1:nrow(candidate_region)) {
    
    # Read the European and African summary statistics
    summary_stat_eur <- fread(paste0(data_dir, "eur_loci_", i, "_summ"))
    summary_stat_eur$PHENONAME = trait
    summary_stat_afr <- fread(paste0(data_dir, "afr_loci_", i, "_summ"))
    summary_stat_afr$PHENONAME = trait
    
    # Combine European and African data
    summary_stat_combined <- full_join(
      summary_stat_eur, 
      summary_stat_afr, 
      by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME"), 
      suffix = c("_eur", "_afr")
    ) %>%
      mutate(Region = i) %>%
      select(c(
        "SNP", "CHR", "POS", "CHR_POS", "PHENONAME", "Region", 
        "ALT", "REF", "MAF_eur", "MAF_afr", "BETA_eur", "BETA_afr", 
        "Z_eur", "Z_afr", "N_eur", "N_afr"
      )) %>% mutate(PHENONAME = trait)
    
    summary_stat_eur_save <- summary_stat_eur %>% mutate(Region = i) %>%  mutate(PHENONAME = trait)
    summary_stat_afr_save <- summary_stat_afr %>% mutate(Region = i) %>%  mutate(PHENONAME = trait)
    
    # Append the combined data frame to the list
    summary_stat_list_eur[[i]] <- summary_stat_eur_save
    summary_stat_list_afr[[i]] <- summary_stat_afr_save
    summary_stat_list[[i]] <- summary_stat_combined
  }
  
  summary_stat_all_eur <- bind_rows(summary_stat_list_eur)
  summary_stat_all_afr <- bind_rows(summary_stat_list_afr)
  summary_stat_all <- bind_rows(summary_stat_list)
  
  fwrite(summary_stat_all_eur,paste0(out_dir,trait,"_REF_EUR"),col.names = T,row.names = F,quote=F)
  fwrite(summary_stat_all_afr,paste0(out_dir,trait,"_REF_AFR"),col.names = T,row.names = F,quote=F)
  fwrite(summary_stat_all,paste0(out_dir,trait,"_REF_ALL"),col.names = T,row.names = F,quote=F)
}



summary_stat_list_all_loci <- list()

for (afr_file in afr_files){
  

  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap/", trait, "_risk_loci"))
  out_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/"),"out/")
  
  summary_stat_trait = fread(paste0(out_dir,trait,"_REF_EUR")) %>%mutate(ID = CHR_POS)
  summary_stat_list_all_loci[[afr_file]] <- summary_stat_trait

  

}

summary_stat_all_combined <- bind_rows(summary_stat_list_all_loci)
fwrite(summary_stat_all_combined,paste0(wrk_dir,"EU_all_combined_bp_updated"),col.names = T,row.names = F,quote=F)
#fwrite(summary_stat_all_combined,paste0(wrk_dir,"EU_all_combined_bmi"),col.names = T,row.names = F,quote=F)





###############For annotation#####################
all_dat<-fread(paste0(wrk_dir,"EU_all_combined_bp_updated"))
#all_dat<-fread(paste0(wrk_dir,"EU_all_combined_bp"))
all_dat<-all_dat%>%distinct(CHR_POS,.keep_all = TRUE)%>%mutate(ID = CHR_POS)

vcf_file<-all_dat%>%select(CHR,POS,ID,REF,ALT,CHR_POS)
annotation_out_name<-paste0(wrk_dir,"EU_AA_bp_updated.vcf")
#annotation_out_name<-paste0(wrk_dir,"EU_AA_bp.vcf")
write.table(vcf_file,annotation_out_name,col.names=T,row.names=F,quote=F)
system(paste0("sed -i '1s/^/#/' ",annotation_out_name))
system(paste0("sed -i '1i###fileformat=VCFv4.0 \' ",annotation_out_name))
most_severe_out<-paste0(wrk_dir,"Most_Severe")



loftee_out<-paste0(wrk_dir,"Loftee")
loftee_out_HC<-paste0(wrk_dir,"Loftee_HC")
loftee_out_LC<-paste0(wrk_dir,"Loftee_LC")
hc_cmd<-paste0("cat ",loftee_out,"|grep LoF=HC>",loftee_out_HC)
system(hc_cmd)
lc_cmd<-paste0("cat ",loftee_out,"|grep LoF=LC>",loftee_out_LC)
system(lc_cmd)
hc_annotation<-read.table(loftee_out_HC,header=F)
lc_annotation<-read.table(loftee_out_LC,header=F)
hc_snp<-unique(hc_annotation$V1)
lc_snp<-unique(lc_annotation$V1)
most_severe_annotation<-fread(most_severe_out,skip=42)
colnames(most_severe_annotation)[1]<-"variant"

missense_snp<-c(hc_snp,lc_snp,most_severe_annotation%>%filter(Consequence=="missense_variant")%>%pull(variant))
synonymous_snp<-most_severe_annotation%>%filter(!(variant%in%missense_snp),Consequence=="synonymous_variant")%>%pull(variant)
utr_5_snp<-most_severe_annotation%>%filter(!(variant%in%c(missense_snp,synonymous_snp)),Consequence=="5_prime_UTR_variant")%>%pull(variant)
utr_3_snp<-most_severe_annotation%>%filter(!(variant%in%c(missense_snp,synonymous_snp,utr_5_snp)),Consequence=="3_prime_UTR_variant")%>%pull(variant)



Promotor_CRE_all<-c()
for(chr in 1:22){
  Promotor_CRE_all<-rbind(Promotor_CRE_all,fread(paste0(wrk_dir,"LDSC_base_annotation/baselineLD.",chr,".annot.gz"))) #This annotation is on build 38
}
Promotor_CRE_all<-Promotor_CRE_all%>%mutate(CHR_POS=paste0(CHR,"_",BP))%>%filter(CHR_POS%in%most_severe_annotation$variant)
promotor_snp<-Promotor_CRE_all%>%filter(!(CHR_POS%in%c(missense_snp,synonymous_snp,utr_5_snp,utr_3_snp)),Promoter_UCSC.extend.500==1)%>%pull(CHR_POS)
cre_snp<-Promotor_CRE_all%>%filter(!(CHR_POS%in%c(missense_snp,synonymous_snp,utr_5_snp,utr_3_snp,promotor_snp)),H3K27ac_Hnisz.extend.500==1,DHS_Trynka.extend.500==1)%>%pull(CHR_POS)




annotation_data<-data.frame("SNP"=most_severe_annotation$variant)
annotation_data<-annotation_data%>%mutate(missense = ifelse(SNP%in%missense_snp,1,0),synonymous = ifelse(SNP%in%synonymous_snp,1,0),utr_5 = ifelse(SNP%in%utr_5_snp,1,0),utr_3 = ifelse(SNP%in%utr_3_snp,1,0),promotor = ifelse(SNP%in%promotor_snp,1,0),CRE = ifelse(SNP%in%cre_snp,1,0))

apply(annotation_data[,2:ncol(annotation_data)],2,table)

#### eQTL processing
liver_ind_gtex_eqtl<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/Liver.v8.independent_eqtls.txt.gz")

blood_ind_gtex_eqtl<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/Whole_Blood.v8.independent_eqtls.txt.gz")

brain_ind_gtex_eqtl<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/Brain_Frontal_Cortex_BA9.v8.independent_eqtls.txt.gz")



format_eqtl<-function(dat,CHR_POS){
  variant_gh38<-do.call(rbind,lapply(dat$variant_id,function(x)unlist(strsplit(x,"_"))))
  df <- data.frame(cbind(variant_gh38[,1], variant_gh38[,2], variant_gh38[,2]))
  colnames(df) <- c('chr', 'start', 'end')
  df$chr<-gsub("chr","",df$chr)
  df$chr_pos<-paste0(df$chr,"_",df$start)
  gwas_eqtl_intersect<-intersect(df$chr_pos,CHR_POS)
  return(gwas_eqtl_intersect)
}

liver_ind_eQTL<-format_eqtl(liver_ind_gtex_eqtl,annotation_data$SNP)
blood_ind_eQTL<-format_eqtl(blood_ind_gtex_eqtl,annotation_data$SNP)
brain_ind_eQTL<-format_eqtl(brain_ind_gtex_eqtl,annotation_data$SNP)


annotation_data<-annotation_data%>%mutate(liver_ind_eQTL = ifelse(SNP%in%liver_ind_eQTL,1,0), blood_ind_eQTL = ifelse(SNP%in%blood_ind_eQTL,1,0), brain_ind_eQTL = ifelse(SNP%in%brain_ind_eQTL,1,0))
#annotation_data<-annotation_data%>%mutate(pituitary_ind_eQTL = ifelse(SNP%in%pituitary_ind_eQTL,1,0), muscle_skeletal_ind_eQTL = ifelse(SNP%in%muscle_skeletal_ind_eQTL,1,0), cells_cultured_fibroblasts_ind_eQTL = ifelse(SNP%in%cells_cultured_fibroblasts_ind_eQTL,1,0))

annotation_data<-annotation_data%>%separate(SNP, into = c("CHR", "POS"), sep = "_", convert = TRUE,remove = FALSE)

#####Phylop and Phast scores
phylop <- getGScores("phyloP100way.UCSC.hg38")
Grange_input<-GRanges(seqnames=paste0("chr",annotation_data$CHR), IRanges(start=annotation_data$POS, width=rep(1,nrow(annotation_data))))
phylop_summary<-gscores(phylop, Grange_input,scores.only=T)

annotation_data$phylop<-phylop_summary$default
annotation_data$phylop[is.na(annotation_data$phylop)]<-"NA"

annotation_data$GERP4<-Promotor_CRE_all$GERP.RSsup4[match(annotation_data$SNP,Promotor_CRE_all$CHR_POS)]
annotation_data$GERP<-Promotor_CRE_all$GERP.NS[match(annotation_data$SNP,Promotor_CRE_all$CHR_POS)]
write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_bp_updated.txt"),col.names=T,row.names=F,quote=F,sep=" ")













wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/"

format_eqtl<-function(dat,CHR_POS){
  variant_gh38<-do.call(rbind,lapply(dat$variant_id,function(x)unlist(strsplit(x,"_"))))
  df <- data.frame(cbind(variant_gh38[,1], variant_gh38[,2], variant_gh38[,2]))
  colnames(df) <- c('chr', 'start', 'end')
  df$chr<-gsub("chr","",df$chr)
  df$chr_pos<-paste0(df$chr,"_",df$start)
  gwas_eqtl_intersect<-intersect(df$chr_pos,CHR_POS)
  return(gwas_eqtl_intersect)
}

annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_BP_updated.txt"))

############################

# For brain eqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  "Brain_Amygdala",
  "Brain_Cortex",
  "Brain_Putamen_basal_ganglia",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Caudate_basal_ganglia",
  "Brain_Hippocampus",
  "Brain_Substantia_nigra",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Hypothalamus",
  "Brain_Cerebellum",
  "Brain_Nucleus_accumbens_basal_ganglia"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))
  
  # Read the eQTL file
  eqtl_data <- fread(file_path)
  
  # Format the eQTL data using your SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_brain_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))

annotation_data<-annotation_data%>%mutate(brain_ind_eQTL = ifelse(SNP%in%all_brain_eqtl_snps,1,0))
table(annotation_data$brain_ind_eQTL)


#annotation_data = annotation_data %>% select(-Human_Promoter_Villar_ExAC.flanking.500_common,-brain_eQTL,-Human_Promoter_Villar_ExAC_common)
############################

        # For heart eqtl

############################

# List of brain tissue names (without file extension)
tissues <- c(
  "Heart_Atrial_Appendage",
  "Heart_Left_Ventricle"
  )

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))
  
  # Read the eQTL file
  eqtl_data <- fread(file_path)
  
  # Format the eQTL data using your SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_heart_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(heart_ind_eQTL = ifelse(SNP%in%all_heart_eqtl_snps,1,0))
table(annotation_data$heart_ind_eQTL)
which(is.na(annotation_data))


############################

# For artery eqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  "Artery_Aorta",
  "Artery_Coronary"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))
  
  # Read the eQTL file
  eqtl_data <- fread(file_path)
  
  # Format the eQTL data using your SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_artery_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(artery_ind_eQTL = ifelse(SNP%in%all_artery_eqtl_snps,1,0))
table(annotation_data$artery_ind_eQTL)


############################

# For brain sqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  "Brain_Amygdala",
  "Brain_Cortex",
  "Brain_Putamen_basal_ganglia",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Caudate_basal_ganglia",
  "Brain_Hippocampus",
  "Brain_Substantia_nigra",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Hypothalamus",
  "Brain_Cerebellum",
  "Brain_Nucleus_accumbens_basal_ganglia"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/GTEx_Analysis_v8_sQTL_independent/"

# Create an empty list to store formatted sqtl results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_sqtls.txt.gz"))
  
  # Read the sqtl file
  sqtl_data <- fread(file_path)
  
  # Format the sqtl data using your SNP annotation data
  formatted_sqtls[[tissue]] <- format_eqtl(sqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_brain_sqtl_snps <- unique(unlist(lapply(formatted_sqtls, function(x) x)))

annotation_data<-annotation_data%>%mutate(brain_ind_sqtl = ifelse(SNP%in%all_brain_sqtl_snps,1,0))
table(annotation_data$brain_ind_sqtl)



############################

# For heart sqtl

############################

# List of brain tissue names (without file extension)
tissues <- c(
  #"Heart_Atrial_Appendage",
  "Heart_Left_Ventricle"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/GTEx_Analysis_v8_sQTL_independent/"

# Create an empty list to store formatted sqtl results for each tissue
formatted_sqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_sqtls.txt.gz"))
  
  # Read the sqtl file
  sqtl_data <- fread(file_path)
  
  # Format the sqtl data using your SNP annotation data
  formatted_sqtls[[tissue]] <- format_eqtl(sqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_heart_sqtl_snps <- unique(unlist(lapply(formatted_sqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(heart_ind_sqtl = ifelse(SNP%in%all_heart_sqtl_snps,1,0))
table(annotation_data$heart_ind_sqtl)



############################

# For artery sqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  #"Artery_Aorta",
  "Artery_Coronary"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/GTEx_Analysis_v8_sQTL_independent/"

# Create an empty list to store formatted sqtl results for each tissue
formatted_sqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_sqtls.txt.gz"))
  
  # Read the sqtl file
  sqtl_data <- fread(file_path)
  
  # Format the sqtl data using your SNP annotation data
  formatted_sqtls[[tissue]] <- format_eqtl(sqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_artery_sqtl_snps <- unique(unlist(lapply(formatted_sqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(artery_ind_sqtl = ifelse(SNP%in%all_artery_sqtl_snps,1,0))
table(annotation_data$artery_ind_sqtl)



############################

# For liver sqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  
      "Liver"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/GTEx_Analysis_v8_sQTL_independent/"

# Create an empty list to store formatted sqtl results for each tissue
formatted_sqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_sqtls.txt.gz"))
  
  # Read the sqtl file
  sqtl_data <- fread(file_path)
  
  # Format the sqtl data using your SNP annotation data
  formatted_sqtls[[tissue]] <- format_eqtl(sqtl_data, annotation_data$SNP)
  
  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}


all_liver_sqtl_snps <- unique(unlist(lapply(formatted_sqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(liver_ind_sqtl = ifelse(SNP%in%all_liver_sqtl_snps,1,0))
table(annotation_data$liver_ind_sqtl)




##Update for rsid

ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_afr_maf_unrelated_biallelic_nodup_extracted.bim")
colnames(ref) = c("CHR","SNP","CM","POS","ALT","REF")
ref$CHR_POS = paste0(ref$CHR,"_",ref$POS)
annotation_data$RSID = ref$SNP[match(annotation_data$SNP,ref$CHR_POS)]


#write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_bmi.txt"),col.names=T,row.names=F,quote=F,sep=" ")

write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_bp_updated"),col.names=T,row.names=F,quote=F,sep=" ")


snplist_annot = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/annot_combine/combined_baselineLF2.2.UKB.annot")

wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/"
annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_bp_updated.txt"))

baseline_annot_bp = snplist_annot %>% filter(SNP %in% annotation_data$RSID) %>% arrange(CHR,BP) %>% distinct(SNP,.keep_all = TRUE)

baseline_annot_bp <- baseline_annot_bp[match(annotation_data$RSID, baseline_annot_bp$SNP), ]

all(baseline_annot_bp$SNP == annotation_data$RSID)



write.table(baseline_annot_bp,paste0(wrk_dir,"/Baseline_Annotation_BP_updated.txt"),col.names=T,row.names=F,quote=F,sep=" ")

#write.table(baseline_annot_bp,paste0(wrk_dir,"/Baseline_Annotation_BMI.txt"),col.names=T,row.names=F,quote=F,sep=" ")

#write.table(baseline_annot_bp,paste0(wrk_dir,"/Baseline_Annotation_Lipids.txt"),col.names=T,row.names=F,quote=F,sep=" ")


baseline_annot_bp = fread(paste0(wrk_dir,"/Baseline_Annotation_BP_updated.txt")) %>% select(-c(1:4))  %>% select(-contains("lowfreq"))
#all(baseline_annot_bp$SNP == annotation_data$RSID)
annotation_data <- annotation_data %>% mutate(across(c(phylop, GERP4, GERP), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
annotation_data_combined = cbind(annotation_data,baseline_annot_bp) %>% select(-RSID)
write.table(annotation_data_combined,paste0(wrk_dir,"/Summarized_Annotation_BP_updated.txt"),col.names=T,row.names=F,quote=F,sep=" ")



baseline_annot_lipids = fread(paste0(wrk_dir,"/Baseline_Annotation_Lipids.txt")) %>% select(-c(1:4))  %>% select(-contains("lowfreq"))
annotation_data_combined = cbind(annotation_data,baseline_annot_lipids) %>% select(-RSID)
write.table(annotation_data_combined,paste0(wrk_dir,"/Summarized_Annotation_lipids.txt"),col.names=T,row.names=F,quote=F,sep=" ")





wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/"
annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation.txt"))
baseline_annot_bp = fread(paste0(wrk_dir,"/Baseline_Annotation_BP.txt")) %>% select(-c(1:4))  %>% select(-contains("lowfreq"))
annotation_data_combined = cbind(annotation_data,baseline_annot_bp) %>% select(-RSID)
write.table(annotation_data_combined,paste0(wrk_dir,"/Summarized_Annotation.txt"),col.names=T,row.names=F,quote=F,sep=" ")






####Add CADD scores####

wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/"
annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_BP_updated_cadd.txt"))

#files <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/cadd/EU_AA_bmi.part", 1:10, "_cadd.tsv.gz")
files <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/cadd/EU_AA_bp_updated.part", 1:10, "_cadd.tsv.gz")

df_list <- lapply(files, function(f) {
  fread(f, header = TRUE, sep = "\t")
})

combined_df <- do.call(rbind, df_list)
colnames(combined_df)[1] = "CHR"
combined_df  = combined_df %>% mutate(SNP = paste0(CHR,"_",Pos)) %>% arrange(CHR,Pos)

annotation_data = annotation_data  %>% arrange(CHR,POS)

all(annotation_data$SNP == combined_df$SNP)

annotation_data$CADD = combined_df$PHRED
fwrite(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_BP_updated_cadd.txt"),col.names=T,row.names=F,quote=F,sep=" ")



annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/Summarized_Annotation_BP_updated.txt")

annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_BP_updated_cadd.txt"))
