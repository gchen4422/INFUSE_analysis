library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)
library(ggbreak)
library(DescTools)
library(coin)
library(susieR)
library(ggrepel)
library(stringr)

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"


#load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
#res_bp <- res_all
#res_bp <- res_bp %>% select(-all_of(colnames(res_bp)[grepl("eQTL", colnames(res_bp)) & !grepl("GTEx_eQTL_MaxCPP_common", colnames(res_bp))]))

# Load BMI data with specific variable name  
#load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData"))
#res_bmi <- res_all
#res_bmi <- res_bmi %>% select(-all_of(colnames(res_bmi)[grepl("eQTL", colnames(res_bmi)) & !grepl("GTEx_eQTL_MaxCPP_common", colnames(res_bmi))]))


#res_all <- rbind(res_bp, res_bmi)
#save(res_all, file = paste0(res_dir, "res_all_pca_bic_modified_llk.RData"))


load(paste0(res_dir,"res_all_pca_bic_modified_llk.RData"))
res_all$Region = ifelse(res_all$Trait == "BMI", res_all$Region+ 84, res_all$Region)
res_all %>% dplyr::group_by(Trait) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
############################################################################################################
#
#   Stats to report
#
#############################################################################################################


#res_all = res_all %>% filter(Trait == "DBP")


# Count unique regions
unique_regions <- res_all %>% summarise(n_unique_regions = n_distinct(Region))
print(unique_regions)

snp_summary <- res_all %>%group_by(Region) %>% summarise(NSNP = n()) %>%summarise(summary(NSNP))
print(snp_summary)




# Calculate position length by region and summarize
pos_length_summary <- res_all %>%group_by(Region) %>%summarise(max_pos = max(as.numeric(POS)),min_pos = min(as.numeric(POS)),pos_length = as.numeric(max_pos) - as.numeric(min_pos)) %>%summarise(summary(pos_length / (1024 * 1024)))
print(pos_length_summary)



# Count unique regions by trait
unique_regions_by_trait <- res_all %>%group_by(Trait) %>%summarise(n_unique_regions = n_distinct(Region))
print(unique_regions_by_trait)




# 95% credible set size 
region_cs<-res_all%>%group_by(Region) %>% summarise(across(c("MESuSiE_cs", "SuSiE_cs","Paintor_cs","Paintor_cs_all","MESuSiE_res_mlk_annot_all_cs","MESuSiE_res_glmnet_annot_all_cs"), ~ sum(.x, na.rm = TRUE)))
region_cs%>%summarise(sum(MESuSiE_cs!=0),sum(SuSiE_cs!=0),sum(Paintor_cs!=0),sum(Paintor_cs_all!=0),sum(MESuSiE_res_mlk_annot_all_cs !=0),sum(MESuSiE_res_glmnet_annot_all_cs !=0))
res_all %>%
  group_by(Region) %>%
  summarise(
    across(c("MESuSiE_cs", "SuSiE_cs", "Paintor_cs","Paintor_cs_all","MESuSiE_res_mlk_annot_all_cs", "MESuSiE_res_glmnet_annot_all_cs"), 
           ~ sum(.x, na.rm = TRUE), .names = "sum_{col}")) %>%
  summarise(
    MESuSiE_cs_median = median(sum_MESuSiE_cs[sum_MESuSiE_cs != 0], na.rm = TRUE),
    SuSiE_cs_median = median(sum_SuSiE_cs[sum_SuSiE_cs != 0], na.rm = TRUE),
    Paintor_cs_median = median(sum_Paintor_cs[sum_Paintor_cs != 0], na.rm = TRUE),
    Paintor_cs_all_median = median(sum_Paintor_cs_all[sum_Paintor_cs_all != 0], na.rm = TRUE),
    MESuSiE_res_mlk_annot_all_cs_median = median(sum_MESuSiE_res_mlk_annot_all_cs[sum_MESuSiE_res_mlk_annot_all_cs != 0], na.rm = TRUE),
    MESuSiE_res_glmnet_annot_all_cs_median = median(sum_MESuSiE_res_glmnet_annot_all_cs[sum_MESuSiE_res_glmnet_annot_all_cs != 0], na.rm = TRUE)
    #MESuSiE_res_mlk_annot_20_baseline_cs_median = median(sum_MESuSiE_res_mlk_annot_20_baseline_cs[sum_MESuSiE_res_mlk_annot_20_baseline_cs != 0], na.rm = TRUE),
    #MESuSiE_res_glmnet_annot_20_baseline_cs_median = median(sum_MESuSiE_res_glmnet_annot_20_baseline_cs[sum_MESuSiE_res_glmnet_annot_20_baseline_cs != 0], na.rm = TRUE)
  )
#MESuSiE identified set overlaps with SuSiE, Paintor credible set
res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_cs==1))
res_all%>%summarise(overlap_prop = sum(MESuSiE_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_cs==1))
#MESuSiE with annot mlk identified set overlaps with SuSiE, Paintor credible set
res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_res_mlk_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_mlk_annot_all_cs==1))
res_all%>%summarise(overlap_prop = sum(MESuSiE_res_mlk_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_mlk_annot_all_cs==1))
#MESuSiE with annot glm identified set overlaps with SuSiE, Paintor credible set
res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_res_glmnet_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_glmnet_annot_all_cs==1))
res_all%>%summarise(overlap_prop = sum(MESuSiE_res_glmnet_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_glmnet_annot_all_cs==1))

# Median |Z-scores| of SNPs within credible set to quantile the strength of detected SNPs
MESuSiE_median_z<-res_all%>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
#MESuSiE_mlk_10_median_z<-res_all%>%filter(MESuSiE_res_mlk_annot_10_baseline_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
#MESuSiE_glm_10_median_z<-res_all%>%filter(MESuSiE_res_glmnet_annot_10_baseline_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
MESuSiE_mlk_median_z<-res_all%>%filter(MESuSiE_res_mlk_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
MESuSiE_glm_median_z<-res_all%>%filter(MESuSiE_res_glmnet_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
SuSiE_median_z<-res_all%>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
Paintor_median_z<-res_all%>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
Paintor_all_median_z<-res_all%>%filter(Paintor_cs_all==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
#Paintor_20_median_z<-res_all%>%filter(Paintor_cs_sub_20==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)


cat(round(c(MESuSiE_median_z,MESuSiE_mlk_median_z,MESuSiE_glm_median_z,SuSiE_median_z,Paintor_median_z,Paintor_all_median_z),2))
# eQTL enrichment by traitf
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
res_all<-res_all%>%mutate(utr_comb_flanking = ifelse((UTR_3_UCSC.flanking.500_common+UTR_5_UCSC.flanking.500_common)>0,1,0))

#ann_col_name<-c("missense", "synonymous", "utr_comb", "promotor", "CRE","brain_ind_eQTL","non_synonymous_common","Conserved_LindbladToh_common","GTEx_eQTL_MaxCPP_common","H3K4me3_Trynka_common","utr_comb_flanking","BLUEPRINT_H3K27acQTL_MaxCPP_common","BLUEPRINT_H3K4me1QTL_MaxCPP_common","BLUEPRINT_DNA_methylation_MaxCPP_common")
#ann_col_name<-c( "non_synonymous_common","synonymous", "promotor","utr_comb","Conserved_LindbladToh_common","H3K4me3_Trynka_common","GTEx_eQTL_MaxCPP_common","brain_ind_eQTL")
ann_col_name<-c( "non_synonymous_common","synonymous", "promotor","utr_comb","Conserved_LindbladToh_common","H3K4me3_Trynka_common","GTEx_eQTL_MaxCPP_common")

# Functions for calculating fold enrichment
calc_fold_enrichment <- function(df, cs_col, ann_col_name) {
  df %>%
    group_by(Region) %>%
    filter(sum(!!sym(cs_col)) != 0) %>%
    group_by(Trait, !!sym(cs_col)) %>%
    summarise(across(ann_col_name, ~ sum(.x, na.rm = TRUE) / n())) %>%
    group_by(Trait) %>%
    summarise(across(ann_col_name, ~ .x[!!sym(cs_col) == 1] / .x[!!sym(cs_col) == 0]))
}

MESuSiE_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_cs", ann_col_name)
#MESuSiE_mlk_10_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_mlk_annot_10_baseline_cs", ann_col_name)
#MESuSiE_glm_10_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_glmnet_annot_10_baseline_cs", ann_col_name)
MESuSiE_mlk_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_mlk_annot_all_cs", ann_col_name)
MESuSiE_glm_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_glmnet_annot_all_cs", ann_col_name)
SuSiE_PIP_ann <- calc_fold_enrichment(res_all, "SuSiE_cs", ann_col_name)
Paintor_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs", ann_col_name)
Paintor_all_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs_all", ann_col_name)
#Paintor_20_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs_sub_20", ann_col_name)



# Combine results
Trait_CS_enrichment <- bind_rows(
  MESuSiE_PIP_ann %>% mutate(Method = "MESuSiE"),
  #MESuSiE_mlk_10_PIP_ann %>% mutate(Method = "MESuSiE_mlk_10"),
  #MESuSiE_glm_10_PIP_ann %>% mutate(Method = "MESuSiE_glm_10"),
  MESuSiE_mlk_PIP_ann %>% mutate(Method = "MESuSiE_mlk"),
  MESuSiE_glm_PIP_ann %>% mutate(Method = "MESuSiE_glm"),
  SuSiE_PIP_ann %>% mutate(Method = "SuSiE"),
  Paintor_PIP_ann %>% mutate(Method = "Paintor"),
  #Paintor_10_PIP_ann %>% mutate(Method = "Paintor_annot"),
  Paintor_all_PIP_ann %>% mutate(Method = "Paintor_annot_all"),
  
  
) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))%>%
  #dplyr::select(Trait,brain_ind_eQTL ,Method )%>%dplyr::rename(eQTL = brain_ind_eQTL)
  dplyr::select(Trait,Method)

# Pivot to long format
Trait_CS_enrichment_long <- Trait_CS_enrichment %>%
  pivot_longer(cols = -c(Method, Trait), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))

Trait_CS_enrichment_long%>%group_by(Method)%>%summarise(min(Prop),max(Prop),mean(Prop))


# Enrichment of 95% credible set without by trait
calc_fold_enrichment_marginal<-function(df, cs_col, ann_col_name) {
  df %>%group_by(Region) %>%
    filter(sum(!!sym(cs_col)) != 0) %>%
    group_by( !!sym(cs_col)) %>%
    summarise(across(ann_col_name, ~ sum(.x, na.rm = TRUE) / n())) %>%
    summarise(across(ann_col_name, ~ .x[!!sym(cs_col) == 1] / .x[!!sym(cs_col) == 0]))
}
MESuSiE_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_cs", ann_col_name)
#MESuSiE_mlk_10_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_res_mlk_annot_10_baseline_cs", ann_col_name)
#MESuSiE_glm_10_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_res_glmnet_annot_10_baseline_cs", ann_col_name)
MESuSiE_mlk_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_res_mlk_annot_all_cs", ann_col_name)
MESuSiE_glm_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_res_glmnet_annot_all_cs", ann_col_name)
SuSiE_PIP_ann <- calc_fold_enrichment_marginal(res_all, "SuSiE_cs", ann_col_name)
Paintor_PIP_ann <- calc_fold_enrichment_marginal(res_all, "Paintor_cs", ann_col_name)
#Paintor_10_PIP_ann <- calc_fold_enrichment_marginal(res_all, "Paintor_cs_sub_baseline", ann_col_name)
Paintor_all_PIP_ann <- calc_fold_enrichment_marginal(res_all, "Paintor_cs_all", ann_col_name)

CS_enrichment <- bind_rows(
  MESuSiE_PIP_ann %>% mutate(Method = "MESuSiE"),
  #MESuSiE_mlk_10_PIP_ann %>% mutate(Method = "MESuSiE_mlk_10"),
  #MESuSiE_glm_10_PIP_ann %>% mutate(Method = "MESuSiE_glm_10"),
  MESuSiE_mlk_PIP_ann %>% mutate(Method = "MESuSiE_mlk"),
  MESuSiE_glm_PIP_ann %>% mutate(Method = "MESuSiE_glm"),
  SuSiE_PIP_ann %>% mutate(Method = "SuSiE"),
  Paintor_PIP_ann %>% mutate(Method = "Paintor"),
  Paintor_all_PIP_ann %>% mutate(Method = "Paintor_annot_all"),
  #Paintor_10_PIP_ann %>% mutate(Method = "Paintor_annot")
) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))%>% 
  dplyr::rename(Missense = missense ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,eQTL = brain_ind_eQTL)


# Pivot to long format
CS_enrichment_long <- CS_enrichment %>%
  pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all"))) %>%
  mutate(Cat = factor(Cat, levels = c("Missense", "Synonymous", "UTR", "Promotor", "CRE","eQTL")))%>%
  mutate(Prop = round(Prop, 2))


# CS_enrichment_long%>%filter(Cat == "Missense")
# CS_enrichment_long%>%filter(Cat == "Synonymous")
# CS_enrichment_long%>%filter(Cat == "UTR")
# CS_enrichment_long%>%filter(Cat == "Promotor")
# CS_enrichment_long%>%filter(Cat == "CRE")
# CS_enrichment_long%>%filter(Cat == "eQTL")


MESuSiE_mlk_highest_cat<-CS_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_mlk")%>%pull(Cat)
CS_enrichment_long%>%filter(Method == "MESuSiE_mlk",Cat %in% MESuSiE_mlk_highest_cat)



MESuSiE_glm_highest_cat<-CS_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_glm")%>%pull(Cat)
CS_enrichment_long%>%filter(Method == "MESuSiE_glm",Cat %in% MESuSiE_glm_highest_cat)

# Enrichment of top 100 signal
res_all<-res_all%>%mutate(Paintor_Signal = ifelse(Paintor_PIP>0.5,1,0))
#res_all<-res_all%>%mutate(Paintor_10_Signal = ifelse(Paintor_PIP_sub_baseline>0.5,1,0))
res_all<-res_all%>%mutate(Paintor_all_Signal = ifelse(Paintor_PIP_all>0.5,1,0))
res_all<-res_all%>%mutate(SuSiE_Signal = case_when(
  susie_EU>0.5&susie_BB >0.5~1,
  susie_EU>0.5&susie_BB<0.5~2,
  susie_EU<0.5&susie_BB>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(MESuSiE_Signal = case_when(
  MESuSiE_PIP_Shared>0.5~1,
  MESuSiE_PIP_WB>0.5~2,
  MESuSiE_PIP_BB>0.5~3,
  .default =0))
# res_all<-res_all%>%mutate(MESuSiE_MLK_10_Signal = case_when(
#   MESuSiE_PIP_Shared_10_sub_mlk>0.5~1,
#   MESuSiE_PIP_WB_10_sub_mlk>0.5~2,
#   MESuSiE_PIP_BB_10_sub_mlk>0.5~3,
#   .default =0))
# res_all<-res_all%>%mutate(MESuSiE_GLM_10_Signal = case_when(
#   MESuSiE_PIP_Shared_10_sub_glmnet>0.5~1,
#   MESuSiE_PIP_WB_10_sub_glmnet>0.5~2,
#   MESuSiE_PIP_BB_10_sub_glmnet>0.5~3,
#   .default =0))
res_all<-res_all%>%mutate(MESuSiE_MLK_Signal = case_when(
  MESuSiE_PIP_Shared_all_mlk>0.5~1,
  MESuSiE_PIP_WB_all_mlk>0.5~2,
  MESuSiE_PIP_BB_all_mlk>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(MESuSiE_GLM_Signal = case_when(
  MESuSiE_PIP_Shared_all_glmnet>0.5~1,
  MESuSiE_PIP_WB_all_glmnet>0.5~2,
  MESuSiE_PIP_BB_all_glmnet>0.5~3,
  .default =0))

top_N_signal = res_all %>%filter(MESuSiE_Signal != 0, MESuSiE_PIP_Either > 0.5) %>%nrow()
top_N_signal = res_all %>%filter(MESuSiE_PIP_Either_all_mlk > 0.5) %>%nrow()

bg_an<-res_all%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/(n()-top_N_signal)))
MESuSiE_Signal_ann<-res_all%>%filter(MESuSiE_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#MESuSiE_MLK_10_Signal_ann<-res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_10_sub_mlk))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_10_sub_mlk)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#MESuSiE_GLM_10_Signal_ann<-res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_10_sub_glmnet))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_10_sub_glmnet)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
MESuSiE_MLK_Signal_ann<-res_all%>%filter(MESuSiE_MLK_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_all_mlk))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
MESuSiE_GLM_Signal_ann<-res_all%>%filter(MESuSiE_GLM_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_all_glmnet))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_glmnet)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
SuSiE_Signal_ann<-res_all%>%filter(SuSiE_Signal!=0)%>% arrange(desc(SuSiE_PIP))%>%top_n(n = top_N_signal, wt = SuSiE_PIP)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
Paintor_Signal_ann<-res_all%>%filter(Paintor_Signal!=0) %>% arrange(desc(Paintor_PIP))%>%top_n(n = top_N_signal, wt = Paintor_PIP)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#Paintor_Signal_10_ann<-res_all%>%filter(Paintor_10_Signal!=0) %>% arrange(desc(Paintor_PIP_sub_baseline))%>%top_n(n = top_N_signal, wt = Paintor_PIP_sub_baseline)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
Paintor_Signal_all_ann<-res_all%>%filter(Paintor_all_Signal!=0) %>% arrange(desc(Paintor_PIP_all))%>%top_n(n = top_N_signal, wt = Paintor_PIP_all)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an

Signal_enrichment <- bind_rows(
  MESuSiE_Signal_ann %>% mutate(Method = "MESuSiE"),
  #MESuSiE_MLK_10_Signal_ann %>% mutate(Method = "MESuSiE_mlk_10"),
  #MESuSiE_GLM_10_Signal_ann %>% mutate(Method = "MESuSiE_glm_10"),
  MESuSiE_MLK_Signal_ann %>% mutate(Method = "MESuSiE_mlk"),
  MESuSiE_GLM_Signal_ann %>% mutate(Method = "MESuSiE_glm"),
  SuSiE_Signal_ann %>% mutate(Method = "SuSiE"),
  Paintor_Signal_ann %>% mutate(Method = "Paintor"),
  #Paintor_Signal_10_ann %>% mutate(Method = "Paintor_annot")
  Paintor_Signal_all_ann %>% mutate(Method = "Paintor_annot_all"),
) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))%>% 
  dplyr::rename(Missense = missense ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,eQTL = brain_ind_eQTL)
  
# Pivot to long format
Signal_enrichment_long <- Signal_enrichment %>%
  pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE", "Paintor","Paintor_annot_all"))) %>%
  mutate(Cat = factor(Cat, levels = c("Missense", "Synonymous", "UTR", "Promotor", "CRE","eQTL")))%>%
  mutate(Prop = round(Prop, 2))


MESuSiE_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE")%>%pull(Cat)
Signal_enrichment_long%>%filter(Method == "MESuSiE",Cat %in% MESuSiE_highest_cat)


MESuSiE_mlk_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_mlk")%>%pull(Cat)
Signal_enrichment_long%>%filter(Method == "MESuSiE_mlk",Cat %in% MESuSiE_mlk_highest_cat)



MESuSiE_glm_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_glm")%>%pull(Cat)
Signal_enrichment_long%>%filter(Method == "MESuSiE_glm",Cat %in% MESuSiE_glm_highest_cat)


# Number of Signal with PIP > 0.5

##Number of signals
res_all%>%group_by(MESuSiE_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(MESuSiE_MLK_10_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(MESuSiE_GLM_10_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(MESuSiE_MLK_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(MESuSiE_GLM_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(SuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Paintor_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(Paintor_10_Signal)%>%summarise(N_signl = n())


res_all%>%group_by(Trait,MESuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,MESuSiE_MLK_10_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,MESuSiE_GLM_10_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,MESuSiE_MLK_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,MESuSiE_GLM_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,SuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,Paintor_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(Trait,Paintor_10_Signal)%>%summarise(N_signl = n())

# Proportion of shared signal
res_all%>%filter(MESuSiE_Signal!=0)%>%group_by(MESuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
#res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>%group_by(MESuSiE_MLK_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
#res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>%group_by(MESuSiE_GLM_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_Signal!=0)%>%group_by(MESuSiE_MLK_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_Signal!=0)%>%group_by(MESuSiE_GLM_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(SuSiE_Signal!=0)%>%group_by(SuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))

res_all%>%filter(MESuSiE_Signal!=0)%>%group_by(Trait,MESuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>%group_by(Trait,MESuSiE_MLK_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>%group_by(Trait,MESuSiE_GLM_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_Signal!=0)%>%group_by(Trait,MESuSiE_MLK_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_Signal!=0)%>%group_by(Trait,MESuSiE_GLM_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(SuSiE_Signal!=0)%>%group_by(Trait,SuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))

###Proportion of shared signal lowering threshold for SuSiE
res_all%>%group_by(Trait)%>%summarise(SuSiE_Shared_n = sum(SuSiE_Shared>0.1),SuSiE_WB_n = sum(susie_EU>0.1&susie_BB<0.1),SuSiE_BB_n = sum(susie_EU<0.1&susie_BB>0.1))%>%summarise(SuSiE_Shared_n/(SuSiE_Shared_n+SuSiE_WB_n+SuSiE_BB_n))
res_all%>%summarise(SuSiE_Shared_n = sum(SuSiE_Shared>0.1),SuSiE_WB_n = sum(susie_EU>0.1&susie_BB<0.1),SuSiE_BB_n = sum(susie_EU<0.1&susie_BB>0.1))%>%summarise(SuSiE_Shared_n/(SuSiE_Shared_n+SuSiE_WB_n+SuSiE_BB_n))




################Signal number by trait, not running for MESuSiE + Annotation#######################
Signal_number<- res_all%>%group_by(Trait)%>%
  summarise(Paintor_Either_n = sum(Paintor_Signal!=0),
            SuSiE_Shared_n = sum(SuSiE_Signal==1),
            SuSiE_WB_n = sum(SuSiE_Signal==2),
            SuSiE_BB_n = sum(SuSiE_Signal==3),
            MESuSiE_Shared_n = sum(MESuSiE_Signal==1),
            MESuSiE_WB_n = sum(MESuSiE_Signal==2),
            MESuSiE_BB_n = sum(MESuSiE_Signal==3))
Signal_number<-Signal_number%>%
  pivot_longer(cols = -c(Trait), names_to = "Cat", values_to = "Num")%>%
  separate(Cat, into = c("Method", "Signal"), sep = "_", extra = "merge") %>%
  mutate(
    Method = case_when(
      str_detect(Method, "MESuSiE") ~ "MESuSiE",
      str_detect(Method, "SuSiE") ~ "SuSiE",
      str_detect(Method, "Paintor") ~ "Paintor",
      TRUE ~ Method
    ),
    Signal = case_when(
      str_detect(Signal, "BB_n") ~ "AFR",
      str_detect(Signal, "WB_n") ~ "EUR",
      str_detect(Signal, "Shared_n") ~ "Shared",
      str_detect(Signal, "Either") ~ "Either",
      TRUE ~ Signal
    )
  )
data.frame(Signal_number)

# Correlation of Shared and ancestry specific signal
# Ancestry-specific signal MESuSiE unique/SuSiE unique/Both
cor_ancestry_MESuSiE <- res_all %>%
  filter(MESuSiE_Signal%in%c(2,3)) %>% filter(!(SuSiE_Signal%in%c(2,3)))%>%
  mutate(Method = "MESuSiE") %>%
  dplyr::select(Beta_WB, Beta_BB, Trait, Method)
cor_ancestry_SuSiE <- res_all %>%
  filter(SuSiE_Signal%in%c(2,3)) %>%filter(!(MESuSiE_Signal%in%c(2,3)))%>%
  mutate(Method = "SuSiE") %>%
  dplyr::select(Beta_WB, Beta_BB, Trait, Method)
cor_ancestry_both <- res_all %>%filter(MESuSiE_Signal%in%c(2,3),SuSiE_Signal%in%c(2,3))%>%
  mutate(Method = "Both") %>%
  dplyr::select(Beta_WB, Beta_BB, Trait, Method)
cor_ancestry <- rbind(cor_ancestry_MESuSiE, cor_ancestry_SuSiE,cor_ancestry_both)
cor_ancestry<-cor_ancestry%>%mutate(Method = factor(Method,levels = c("MESuSiE","SuSiE","Both")))	
###Ancestry specific signal marginal correlations
cor_ancestry %>%
  group_by(Method, Trait) %>%
  summarise(cor(Beta_WB, Beta_BB))
res_all %>%filter(MESuSiE_Signal%in%c(2,3))%>%
  summarise(cor(Beta_WB, Beta_BB))
cor_ancestry  %>%
  group_by(Method) %>%
  summarise(cor(Beta_WB, Beta_BB))
###Shared signal marginal correlations
res_all %>%filter(MESuSiE_Signal==1) %>%group_by(Trait) %>%summarise(cor(Beta_WB, Beta_BB))
res_all %>%filter(MESuSiE_Signal==1) %>%summarise(cor(Beta_WB, Beta_BB))
#	Test Correlation Difference ancestry vs shared
library(boot)
function_cor <- function(data, i){
  d2 <- data[i,] 
  return(cor(d2$Beta_WB, d2$Beta_BB))
}
set.seed(1128)
ancestry_signal_data<-res_all%>%filter(MESuSiE_Signal%in%c(2,3))
bootstrap_correlation_ancestry <- boot(ancestry_signal_data,function_cor,R=1000)
shared_signal_data<-res_all%>%filter(MESuSiE_Signal==1)
bootstrap_correlation_shared <- boot(shared_signal_data,function_cor,R=1000)
res_all%>%filter(MESuSiE_Signal==1)%>%group_by(Trait)%>%summarise(margin_cor = cor(Beta_WB,Beta_BB))
res_all%>%filter(MESuSiE_Signal==1)%>%summarise(margin_cor = cor(Beta_WB,Beta_BB))
cor_dif<-mean(bootstrap_correlation_shared$t-bootstrap_correlation_ancestry$t)
cor_dif_se<-sqrt(var(bootstrap_correlation_shared$t-bootstrap_correlation_ancestry$t))
cor_dif_P<-2*pnorm(-abs(cor_dif/cor_dif_se))

# MAF Diff
maf_shared<-data.frame(res_all)%>%filter(MESuSiE_Signal==1)%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
maf_ancestry<-data.frame(res_all)%>%filter(MESuSiE_Signal%in%c(2,3))%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
wilcox.test(maf_ancestry, maf_shared, alternative = "two.sided")	

#alternative hypothesis: true location shift is not equal to 0
median(maf_shared) 
median(maf_ancestry) 



maf_shared<-data.frame(res_all)%>%filter(MESuSiE_MLK_Signal == 1)%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
maf_ancestry<-data.frame(res_all)%>%filter(MESuSiE_MLK_Signal %in%c(2,3))%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
wilcox.test(maf_ancestry, maf_shared, alternative = "two.sided")	

median(maf_shared)
median(maf_ancestry)  





#####################################################################
#
#     Plot directory and theme
#
####################################################################
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/"
system(paste0("mkdir -p ",plot_dir))
custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
    strip.background = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7, face="bold"),
    plot.title = element_text(size=7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}
################################################
#
#		Set Size/Z-score/eQTL 
#
#
###############################################
################################################
#
#		Set SiZe Part
#
###############################################		
###Median set size by Trait
all_sets_info<-data.frame(res_all%>%group_by(Trait,Region) %>% summarise(across(c("MESuSiE_cs","MESuSiE_res_mlk_annot_all_cs","MESuSiE_res_glmnet_annot_all_cs","SuSiE_cs","Paintor_cs","Paintor_cs_all"), ~ sum(.x, na.rm = TRUE))))%>%filter(MESuSiE_cs!=0, , MESuSiE_res_mlk_annot_all_cs!=0,MESuSiE_res_glmnet_annot_all_cs != 0,SuSiE_cs!=0,Paintor_cs!=0,Paintor_cs_all !=0) #Median Set Size across all locus
all_sets_info_long<-all_sets_info%>%pivot_longer(!(Trait|Region), names_to = "Method", values_to = "Count")
all_sets_info_long$Method<-factor(all_sets_info_long$Method,levels=c("MESuSiE_cs","MESuSiE_res_mlk_annot_all_cs","MESuSiE_res_glmnet_annot_all_cs","SuSiE_cs","Paintor_cs","Paintor_cs_all"))
levels(all_sets_info_long$Method)<-c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot")

# all_sets_info_long_save <- all_sets_info_long %>%
#   mutate(causal = num_causal,
#          h2 = h2_num)


p_set = ggplot(data =all_sets_info_long,aes(x = Trait, y=Count,fill=Method))+geom_boxplot(aes(x = Trait,fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray")+scale_fill_manual(values=c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot"="#d17f45"),guide=FALSE)
p_set =p_set + theme_bw() + xlab("") +ylab("Set Size")+coord_cartesian(ylim=c(0,300))
p_set= p_set+custom_theme()

################################################
#
#		Z-score Part
#
###############################################		
MESuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))
MESuSiE_mlk_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_res_mlk_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
MESuSiE_glm_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_res_glmnet_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
SuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(SuSiE_cs==1)%>%summarise(zmax =median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
Paintor_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
#Paintor_10_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs_sub_baseline==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
Paintor_all_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs_all==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB ))))%>%pull(zmax)
set_size_z_info<-data.frame(cbind(MESuSiE_cs_Z,MESuSiE_mlk_cs_Z,MESuSiE_glm_cs_Z,SuSiE_cs_Z,Paintor_cs_Z,Paintor_all_cs_Z))
colnames(set_size_z_info)<-c("Trait",c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot"))
set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(Trait), names_to = "Method", values_to = "Z")%>%mutate(Method = factor(Method, levels=c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot")))

p_z = ggplot(data = set_size_z_info_long,aes(x = Trait, y=Z,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#023e8a", "MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot"="#d17f45"))
p_z = p_z + geom_text(label = round(set_size_z_info_long$Z,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)
p_z = p_z + theme_bw() + xlab("") +ylab("Median |Z|")+ ylim(0,max(round(set_size_z_info_long$Z,2)+1))
p_z = p_z +custom_theme()
################################################
#
#		eQTL enrichment 
#
#
###############################################	
p_eQTL <- ggplot(Trait_CS_enrichment_long, aes(x = Trait, y = Prop, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +scale_fill_manual(values = c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(,label = round(Trait_CS_enrichment_long$Prop,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  xlab("") + ylab("eQTL Fold Enrichment") + ylim(0,max(round(Trait_CS_enrichment_long$Prop))+1)+
  theme_bw() + custom_theme()

p_out<-p_set/p_z/p_eQTL+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1.5,1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 6,face="bold"))
p_out
ggsave(paste0(plot_dir,"set_size_eqtl_bp.pdf"),p_out,width=140,height =150,dpi=500,units='mm')
################################################################################
#
#
#     Functional Annotation enrichment for 95% credible set SNPS
#
#
################################################################################
p_signal <- ggplot(data = CS_enrichment_long,aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c("MESuSiE"="#023e8a", "MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 5/14*5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("") + ylab("Fold Enrichment Credible Set") +ylim(0,round(max(CS_enrichment_long$Prop))+1) +
  theme_bw() + custom_theme()
p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 6,face="bold"))
ggsave(paste0(plot_dir,"Enrichment_Comparison_Set_CS.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')

################################################################################
#
#
#     Functional Annotation enrichment for top 100 PIP SNPs
#
#
################################################################################	
p_signal <- ggplot(data = Signal_enrichment_long,aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 5/14*5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("") + ylab("Fold Enrichment Top Signal") +ylim(0,round(max(Signal_enrichment_long$Prop))+1) +
  theme_bw() + custom_theme()
p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 7,face="bold"))
ggsave(paste0(plot_dir,"Enrichment_Comparison_Set_Top.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')


############################################################################
#
#
#              Proportion of Signal Plot
#
#
############################################################################

Signal_number<- res_all%>%group_by(Trait)%>%
  summarise(Paintor_Either_n = sum(Paintor_Signal!=0),
            #Paintor10_Either_n = sum(Paintor_10_Signal!=0),
            SuSiE_Shared_n = sum(SuSiE_Signal==1),
            SuSiE_WB_n = sum(SuSiE_Signal==2),
            SuSiE_BB_n = sum(SuSiE_Signal==3),
            MESuSiE_Shared_n = sum(MESuSiE_Signal==1),
            MESuSiE_WB_n = sum(MESuSiE_Signal==2),
            MESuSiE_BB_n = sum(MESuSiE_Signal==3),
            #MESuSiEmlk10_Shared_n = sum(MESuSiE_MLK_10_Signal==1),
            #MESuSiEmlk10_WB_n = sum(MESuSiE_MLK_10_Signal==2),
            #MESuSiEmlk10_BB_n = sum(MESuSiE_MLK_10_Signal==3),
            #MESuSiEglm10_Shared_n = sum(MESuSiE_GLM_10_Signal==1),
            #MESuSiEglm10_WB_n = sum(MESuSiE_GLM_10_Signal==2),
            #MESuSiEglm10_BB_n = sum(MESuSiE_GLM_10_Signal==3),
            MESuSiEmlk_Shared_n = sum(MESuSiE_MLK_Signal==1),
            MESuSiEmlk_WB_n = sum(MESuSiE_MLK_Signal==2),
            MESuSiEmlk_BB_n = sum(MESuSiE_MLK_Signal==3),
            MESuSiEglm_Shared_n = sum(MESuSiE_GLM_Signal==1),
            MESuSiEglm_WB_n = sum(MESuSiE_GLM_Signal==2),
            MESuSiEglm_BB_n = sum(MESuSiE_GLM_Signal==3))

Signal_number<-Signal_number%>%
  pivot_longer(cols = -c(Trait), names_to = "Cat", values_to = "Num")%>%
  separate(Cat, into = c("Method", "Signal"), sep = "_", extra = "merge") %>%
  mutate(
    Method = case_when(
      #str_detect(Method, "MESuSiEmlk10") ~ "MESuSiE_mlk_10",
      #str_detect(Method, "MESuSiEglm10") ~ "MESuSiE_glm_10",
      str_detect(Method, "MESuSiEmlk") ~ "MESuSiE_mlk",
      str_detect(Method, "MESuSiEglm") ~ "MESuSiE_glm",
      str_detect(Method, "MESuSiE") ~ "MESuSiE",
      str_detect(Method, "SuSiE") ~ "SuSiE",
      #str_detect(Method, "Paintor10") ~ "Paintor_annot",
      str_detect(Method, "Paintor") ~ "Paintor",
      TRUE ~ Method
    ),
    Signal = case_when(
      str_detect(Signal, "BB_n") ~ "AFR",
      str_detect(Signal, "WB_n") ~ "EUR",
      str_detect(Signal, "Shared_n") ~ "Shared",
      str_detect(Signal, "Either") ~ "Either",
      TRUE ~ Signal
    )
  )
Signal_number<-Signal_number%>%group_by(Method,Trait)%>%mutate(prop = Num/sum(Num)*100,ypos = cumsum(prop)- 0.5*prop)%>%mutate(label = paste0(Signal," ",Num))
Signal_number<-Signal_number%>%mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE", "Paintor")),Signal = factor(Signal, levels = c("EUR","AFR","Shared","Either")))

signal_num_plot<-ggplot(Signal_number, aes(x="", y=prop, fill=Signal)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = label),  size=7/14*5,color="white") +
  scale_fill_manual(values=c("#ABDB9F","#F2C1B6","#6162B0","gray"))+
  facet_grid(vars(Trait),vars(Method),labeller=label_parsed)+theme(
    strip.text.x = element_text(size = 7,face="bold"),
    strip.text.y = element_text(size = 7,face="bold"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())
ggsave(paste0(plot_dir,"signal_prop_pie_bmi.pdf"),signal_num_plot,dpi=500,height=130,width  =250,unit='mm')

###########################################################################################
#
#
#         Correlation of ancestry-specific signal by MESuSiE and SuSiE
#
#
#########################################################################################

betaplot<- ggscatter(cor_ancestry, x = "Beta_WB", y = "Beta_BB",
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "red", fill = "lightgray",size = 0.5), # Customize reg. line
                     conf.int = FALSE, # Add confidence interval
                     size = 0.5
)+ stat_cor(method = "pearson",aes(label = ..r.label..),size = 5/14*5)+xlab("UKBB")+ylab("GLGC")+
  theme_bw() + custom_theme()+
  facet_wrap(vars(Method),ncol=3)
ggsave(paste0(plot_dir,"beta_cor.pdf"),betaplot,dpi=500,height=60,width=90,units="mm")


#################################################################
#
#
#               MAF and Conservative Score Part
#
#
##################################################################
################################################
#
#		Phylop/GERP Score Part
#
#
################################################		

res_all<- res_all %>%
  mutate(
    MESuSiE_Conservative = case_when(
      MESuSiE_PIP_WB>0.5|MESuSiE_PIP_BB>0.5 ~ "Ancestry-specific",
      MESuSiE_PIP_Shared>0.5 ~ "Shared",
      TRUE ~ "None"
    ),
    MESuSiE_mlk_Conservative = case_when(
      MESuSiE_PIP_WB_all_mlk>0.5|MESuSiE_PIP_BB_all_mlk>0.5 ~ "Ancestry-specific",
      MESuSiE_PIP_Shared_all_mlk>0.5 ~ "Shared",
      TRUE ~ "None"
    ),
    MESuSiE_glm_Conservative = case_when(
      MESuSiE_PIP_WB_all_glmnet>0.5|MESuSiE_PIP_BB_all_glmnet>0.5 ~ "Ancestry-specific",
      MESuSiE_PIP_Shared_all_glmnet>0.5 ~ "Shared",
      TRUE ~ "None"
    ),
    SuSiE_Conservative = case_when(
      (susie_EU>0.5&susie_BB<0.5)|(susie_EU>0.5&susie_BB<0.5) ~ "Ancestry-specific",
      SuSiE_Shared>0.5 ~ "Shared",
      TRUE ~ "None"
    )
  )%>%
  mutate(
    
    MESuSiE_Conservative = factor(MESuSiE_Conservative, levels = c("None", "Ancestry-specific", "Shared")),
    MESuSiE_mlk_Conservative = factor(MESuSiE_mlk_Conservative, levels = c("None", "Ancestry-specific", "Shared")),
    MESuSiE_glm_Conservative = factor(MESuSiE_glm_Conservative, levels = c("None", "Ancestry-specific", "Shared")),
    SuSiE_Conservative = factor(SuSiE_Conservative, levels = c("None", "Ancestry-specific", "Shared"))
  )



res_all%>%group_by(MESuSiE_mlk_Conservative)%>%summarise(mean(abs(phylop),na.rm=T)) 
#The absolute values of phyloP scores represent –log p-values under a null hypothesis of neutral evolution.
res_all%>%group_by(MESuSiE_mlk_Conservative)%>%summarise(mean(abs(GERP),na.rm=T)) 
res_all%>%group_by(MESuSiE_mlk_Conservative)%>%summarise(mean(abs(CADD),na.rm=T)) 


####Phylop
grouped_data_phylop <- res_all %>%
  mutate(abs_phylop = abs(phylop)) %>%
  select(MESuSiE_Conservative, abs_phylop)%>%
  group_by(MESuSiE_Conservative) %>%
  summarise(phylop_values = list(abs_phylop))%>%
  pull(phylop_values)

jt_test_phylop <- JonckheereTerpstraTest(grouped_data_phylop,nperm = 100)

grouped_data_GERP <- res_all %>%
  mutate(GERP = abs(GERP)) %>%
  select(MESuSiE_Conservative, GERP)%>%
  group_by(MESuSiE_Conservative) %>%
  summarise(GERP_values = list(GERP))%>%
  pull(GERP_values)

jt_test_GERP <- JonckheereTerpstraTest(grouped_data_GERP,nperm = 100)


####Phylop
grouped_data_phylop <- res_all %>%
  mutate(abs_phylop = abs(phylop)) %>%
  select(MESuSiE_mlk_Conservative, abs_phylop)%>%
  group_by(MESuSiE_mlk_Conservative) %>%
  summarise(phylop_values = list(abs_phylop))%>%
  pull(phylop_values)

jt_test_phylop <- JonckheereTerpstraTest(grouped_data_phylop,nperm = 100)

grouped_data_GERP <- res_all %>%
  mutate(GERP = abs(GERP)) %>%
  select(MESuSiE_mlk_Conservative, GERP)%>%
  group_by(MESuSiE_mlk_Conservative) %>%
  summarise(GERP_values = list(GERP))%>%
  pull(GERP_values)

jt_test_GERP <- JonckheereTerpstraTest(grouped_data_GERP,nperm = 100)


grouped_data_CADD <- res_all %>%
  mutate(CADD = abs(CADD)) %>%
  select(MESuSiE_mlk_Conservative, CADD)%>%
  group_by(MESuSiE_mlk_Conservative) %>%
  summarise(CADD_values = list(CADD))%>%
  pull(CADD_values)


jt_test_CADD <- JonckheereTerpstraTest(grouped_data_CADD,nperm = 100)

#########Conservation Plot
conserve_data<-test%>%mutate(Signal =MESuSiE_mlk_Conservative, Phylop =abs(phylop) )%>%select(Signal,Phylop,GERP,CADD)

conserve_data<-res_all%>%mutate(Signal =MESuSiE_mlk_Conservative, Phylop =abs(phylop) )%>%select(Signal,Phylop,GERP,CADD)
conserve_data_long<-conserve_data%>%pivot_longer(!(Signal), names_to = "Method", values_to = "Score")
conserve_data_long$Method<-factor(conserve_data_long$Method,levels = c("Phylop","GERP","CADD"))
df <-conserve_data_long %>%group_by(Method,Signal) %>%summarize(ymin = min(Score,na.rm = T), lower = quantile(Score, .25,na.rm = T),  middle = mean(Score,na.rm = T), upper = quantile(Score, .75,na.rm = T),ymax = max(Score,na.rm = T)) 
p_set_phylop = df %>% ggplot(aes(x = factor(Method), fill=Signal)) + geom_boxplot(fatten = 2,aes(ymin = ymin, ymax = ymax, lower = lower, upper = upper, middle = middle), stat = 'identity')+scale_fill_manual(values=c("None"="#DAFFED","Ancestry-specific"="#fffbdb","Shared" = "#7776bc"))
p_set_phylop =p_set_phylop + theme_bw() + xlab("") +ylab("Scores")
p_set_phylop= p_set_phylop  + custom_theme()
p_set_phylop= p_set_phylop+ theme(legend.position="bottom")		

################################################
#
#		MAF Part
#
#
################################################	


maf_data_plot<-data.frame("MAF" = c(maf_shared,maf_ancestry),"Signal"=c(rep("Shared",length(maf_shared)),rep("Ancestry-specific",length(maf_ancestry))))
maf_data_plot$Signal<-factor(maf_data_plot$Signal)
maf_data_plot$Signal<-factor(maf_data_plot$Signal,levels = c(levels(maf_data_plot$Signal),"None"))

p_set = ggplot(data =maf_data_plot,aes(x = factor(Signal), y=MAF,fill=Signal))+geom_boxplot(aes(x = factor(Signal),fill=Signal),outlier.size = 0.01) +scale_fill_manual(values=c("None"="#DAFFED","Ancestry-specific"="#fffbdb","Shared" = "#7776bc"),guide=FALSE)
p_set = p_set + stat_compare_means(size= 5/14*5)
p_set =p_set + theme_bw() + xlab("") +ylab("MAF difference")
p_set= p_set + custom_theme()
p_set= p_set + theme(legend.position="none")	
p_out<-p_set+p_set_phylop+plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect") & theme(legend.position = "bottom",plot.tag = element_text(size = 7,face="bold"))

ggsave(paste0(plot_dir,"MAF_Conservative_mlk_combined.pdf"),p_out,units = 'mm',dpi=500,width = 120,height =120)
