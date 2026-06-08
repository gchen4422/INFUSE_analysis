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

all_sets <- list()
all_CS_enrichment <- list()
all_Signal_enrichment <- list()
all_true_signal_ann <- list()


for(num_causal in 1:3){
  for(h2_num in 1:2){
    
    wrk_dir   <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_", num_causal, "/")
    out_dir<- paste0(wrk_dir, "out/sparse_pca_modified_llk/")
    load(paste0(out_dir,"res_simulation_h2_",h2_num,"_pca_bic.RData"))
    res_all$Trait = "simulation"
    colnames(res_all)[c(18,23)] =c("MESuSiE_res_glmnet_annot_all_cs","MESuSiE_res_mlk_annot_all_cs")
    ############################################################################################################
    #
    #   Stats to report
    #
    #############################################################################################################
    
    # Count unique regions
    unique_regions <- res_all %>% summarise(n_unique_regions = n_distinct(Region))
    print(unique_regions)

    snp_summary <- res_all %>%group_by(Region) %>% summarise(NSNP = n()) %>%summarise(summary(NSNP))
    print(snp_summary)

    
    # Calculate position length by region and summarize
    pos_length_summary <- res_all %>%group_by(Region) %>%summarise(max_pos = max(as.numeric(POS)),min_pos = min(as.numeric(POS)),pos_length = as.numeric(max_pos) - as.numeric(min_pos)) %>%summarise(summary(pos_length / (1024 * 1024)))
    print(pos_length_summary)

    
    

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
    #res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_cs==1))
    #res_all%>%summarise(overlap_prop = sum(MESuSiE_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_cs==1))
    #MESuSiE with annot mlk identified set overlaps with SuSiE, Paintor credible set
    #res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_res_mlk_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_mlk_annot_all_cs==1))
    #res_all%>%summarise(overlap_prop = sum(MESuSiE_res_mlk_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_mlk_annot_all_cs==1))
    #MESuSiE with annot glm identified set overlaps with SuSiE, Paintor credible set
    #res_all%>%group_by(Trait)%>%summarise(overlap_prop = sum(MESuSiE_res_glmnet_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_glmnet_annot_all_cs==1))
    #res_all%>%summarise(overlap_prop = sum(MESuSiE_res_glmnet_annot_all_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_res_glmnet_annot_all_cs==1))
    
    # Median |Z-scores| of SNPs within credible set to quantile the strength of detected SNPs
    MESuSiE_median_z<-res_all%>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    #MESuSiE_mlk_10_median_z<-res_all%>%filter(MESuSiE_res_mlk_annot_10_baseline_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    #MESuSiE_glm_10_median_z<-res_all%>%filter(MESuSiE_res_glmnet_annot_10_baseline_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    MESuSiE_mlk_median_z<-res_all%>%filter(MESuSiE_res_mlk_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    MESuSiE_glm_median_z<-res_all%>%filter(MESuSiE_res_glmnet_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    SuSiE_median_z<-res_all%>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    Paintor_median_z<-res_all%>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    Paintor_all_median_z<-res_all%>%filter(Paintor_cs_all==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    #Paintor_20_median_z<-res_all%>%filter(Paintor_cs_sub_20==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    
    
    cat(round(c(MESuSiE_median_z,MESuSiE_mlk_median_z,MESuSiE_glm_median_z,SuSiE_median_z,Paintor_median_z,Paintor_all_median_z),2))
    
    
    ann_col_name<-c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common")
    
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
    
    #res_all = res_all %>% mutate(cs_flag_true = as.integer(!!sym("Signal") != 0))
    #res_all = res_all %>% mutate(cs_flag_true = as.integer(!!sym("Signal") == 3))
    #res_all <- res_all %>%mutate(cs_flag_true = as.integer(!!sym("Signal") %in% c(1, 2)))
    res_all = res_all %>% mutate(cs_flag_true = as.integer(!!sym("Signal") != 0))
    
    true_Signal_ann = calc_fold_enrichment(res_all, "cs_flag_true", ann_col_name)
    MESuSiE_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_cs", ann_col_name)
    #MESuSiE_mlk_10_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_mlk_annot_10_baseline_cs", ann_col_name)
    #MESuSiE_glm_10_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_glmnet_annot_10_baseline_cs", ann_col_name)
    MESuSiE_mlk_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_mlk_annot_all_cs", ann_col_name)
    MESuSiE_glm_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_res_glmnet_annot_all_cs", ann_col_name)
    SuSiE_PIP_ann <- calc_fold_enrichment(res_all, "SuSiE_cs", ann_col_name)
    Paintor_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs", ann_col_name)
    #Paintor_10_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs_sub_baseline", ann_col_name)
    Paintor_all_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs_all", ann_col_name)
    
    
    
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
      
      
    ) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))#%>%
    #dplyr::select(Trait,cells_cultured_fibroblasts_ind_eQTL ,Method )%>%dplyr::rename(eQTL = cells_cultured_fibroblasts_ind_eQTL)
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
    ) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))#%>% 
    #dplyr::rename(Missense = missense ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,eQTL = cells_cultured_fibroblasts_ind_eQTL)
    
    
    # Pivot to long format
    CS_enrichment_long <- CS_enrichment %>%
      pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
      mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all"))) %>%
      mutate(Cat = factor(Cat, levels = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common")))%>%
      mutate(Prop = round(Prop, 2))
    
    CS_enrichment_long_save <- CS_enrichment_long %>%
      mutate(causal = num_causal,
             h2 = h2_num,
             enrichment_type = "CS")
    
    
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
    
    top_N_signal <- ifelse(num_causal == 3, 150, 100)
    
    bg_an<-res_all%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))
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
      
    ) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))#%>% 
    #dplyr::rename(Missense = missense ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,eQTL = cells_cultured_fibroblasts_ind_eQTL)
    # Pivot to long format
    Signal_enrichment_long <- Signal_enrichment %>%
      pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
      mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all"))) %>%
      mutate(Cat = factor(Cat, levels =c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common")))%>%
      mutate(Prop = round(Prop, 2))
    
    
    Signal_enrichment_long_save <- Signal_enrichment_long %>%
      mutate(causal = num_causal,
             h2 = h2_num,
             enrichment_type = "Signal")
    
    
    MESuSiE_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE")%>%pull(Cat)
    Signal_enrichment_long%>%filter(Method == "MESuSiE",Cat %in% MESuSiE_highest_cat)
    

    
    MESuSiE_mlk_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_mlk")%>%pull(Cat)
    Signal_enrichment_long%>%filter(Method == "MESuSiE_mlk",Cat %in% MESuSiE_mlk_highest_cat)
    

    
    
    MESuSiE_glm_highest_cat<-Signal_enrichment_long%>%group_by(Cat)%>%summarise(Method = Method[which.max(Prop)])%>%filter(Method == "MESuSiE_glm")%>%pull(Cat)
    Signal_enrichment_long%>%filter(Method == "MESuSiE_glm",Cat %in% MESuSiE_glm_highest_cat)
    

    
    
    
    #####################################################################
    #
    #     Plot directory and theme
    #
    ####################################################################
    plot_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/plot/Enrichment/"
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
    
    all_sets_info_long_save <- all_sets_info_long %>%
      mutate(causal = num_causal,
             h2 = h2_num)
    
    
    p_set = ggplot(data =all_sets_info_long,aes(x = Trait, y=Count,fill=Method))+geom_boxplot(aes(x = Trait,fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray")+scale_fill_manual(values=c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot"="#d17f45"),guide=FALSE)
    p_set =p_set + theme_bw() + xlab("") +ylab("Set Size")+coord_cartesian(ylim=c(0,100))
    p_set= p_set+custom_theme()
    
    ################################################
    #
    #		Z-score Part
    #
    ###############################################		
    MESuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))
    MESuSiE_mlk_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_res_mlk_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    MESuSiE_glm_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_res_glmnet_annot_all_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    SuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(SuSiE_cs==1)%>%summarise(zmax =median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    Paintor_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    #Paintor_10_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs_sub_baseline==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
    Paintor_all_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs_all==1)%>%summarise(zmax = median(pmax(abs(zscore_1),abs(zscore_2 ))))%>%pull(zmax)
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
    #ggsave(paste0(plot_dir,"set_size_eqtl.pdf"),p_out,width=140,height =150,dpi=500,units='mm')
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
    ggsave(paste0(plot_dir,"Causal_num_",num_causal,"_h2_",h2_num,"_Enrichment_Comparison_Set_CS.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')
    
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
    #ggsave(paste0(plot_dir,"Causal_num_",num_causal,"_h2_",h2_num,"_Enrichment_Comparison_Set_Top.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')
    
    all_sets[[paste0("causal", num_causal, "_h2", h2_num)]] <- all_sets_info_long_save
    all_CS_enrichment[[paste0("causal", num_causal, "_h2", h2_num)]] <- CS_enrichment_long_save
    all_Signal_enrichment[[paste0("causal", num_causal, "_h2", h2_num)]] <- Signal_enrichment_long_save
    key <- paste0("causal", num_causal, "_h2", h2_num)
    true_Signal_ann = true_Signal_ann %>% mutate(causal = num_causal, h2 = h2_num)
    all_true_signal_ann[[key]] <- true_Signal_ann
    
    }
}


###########Summary plot for average credible set size###########
combined_all_sets <- bind_rows(all_sets)
all_sets_info_long_Combined <-combined_all_sets%>%group_by(Method,Region,Trait)%>%summarise(Count= round(mean(Count),2))
all_sets_info_long_Combined$Method<-factor(all_sets_info_long_Combined$Method,levels =c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot"))

p_set = ggplot(data =all_sets_info_long_Combined,aes(x = Trait, y=Count,fill=Method))+geom_boxplot(aes(x = Trait,fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray")+scale_fill_manual(values=c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot"="#d17f45"),guide=FALSE)
p_set =p_set + theme_bw() + xlab("") +ylab("Set Size")+coord_cartesian(ylim=c(0,100))
p_set= p_set+custom_theme()


###########Summary plot for average enrichment of snps in CS###########

combined_CS_enrichment <- bind_rows(all_CS_enrichment)
CS_enrichment_Combined<-combined_CS_enrichment%>%group_by(Method,Cat)%>%summarise(Prop= round(mean(Prop),2))
CS_enrichment_Combined$Cat<-factor(	 CS_enrichment_Combined$Cat,levels = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common"))
CS_enrichment_Combined$Method<-factor(CS_enrichment_Combined$Method,levels =c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all"))

p_signal <- ggplot(data = CS_enrichment_Combined,aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c("MESuSiE"="#023e8a", "MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 5/14*5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("") + ylab("Fold Enrichment Credible Set") +ylim(0,round(max(CS_enrichment_Combined$Prop))+1) +
  theme_bw() + custom_theme()
p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 6,face="bold"))
ggsave(paste0(plot_dir,"Summary_Enrichment_Comparison_Set_CS.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')



###########Summary plot for average enrichment of top snps###########
combined_Signal_enrichment <- bind_rows(all_Signal_enrichment)
Signal_enrichment_Combined<-combined_Signal_enrichment%>%group_by(Method,Cat)%>%summarise(Prop= round(mean(Prop),2))
Signal_enrichment_Combined$Cat<-factor(	 Signal_enrichment_Combined$Cat,levels = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common"))
Signal_enrichment_Combined$Method<-factor(Signal_enrichment_Combined$Method,levels =c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all"))

p_signal <- ggplot(data = Signal_enrichment_Combined,aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 5/14*5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("") + ylab("Fold Enrichment Top Signal") +ylim(0,round(max(Signal_enrichment_Combined$Prop))+1) +
  theme_bw() + custom_theme()
p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 7,face="bold"))
ggsave(paste0(plot_dir,"Summary_Enrichment_Comparison_Set_Top.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')



###########Summary plot for true enrichment of causal snps###########
true_Signal_ann_df <- bind_rows(all_true_signal_ann)
True_enrichment_long <- true_Signal_ann_df %>%
  pivot_longer(
    cols      = -c(Trait, causal, h2),
    names_to  = "Cat",
    values_to = "Prop"
  ) %>%
  mutate(
    Cat  = factor(
      Cat,
      levels = c(
        "H3K27ac_Hnisz_common",
        "DHS_Trynka_common",
        "H3K4me3_Trynka_common",
        "Conserved_LindbladToh_common",
        "non_synonymous_common"
      )
    ),
    Prop = round(Prop, 2)
  )


True_enrichment_mean = True_enrichment_long%>%group_by(Cat)%>%summarise(Prop= round(mean(Prop),2))
True_enrichment_mean
save(True_enrichment_mean,True_enrichment_long,file = paste0(plot_dir,"/True_enrichment_w2_either.Rdata"))
rm(list = setdiff(ls(), "True_enrichment_long"))