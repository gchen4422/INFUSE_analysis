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
library(purrr)
library(susieR)
library(ggrepel)
library(stringr)


res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
load(paste0(res_dir,"res_all_pca_bic_modified_llk.RData"))
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
res_all<-res_all%>%mutate(utr_comb_flanking = ifelse((UTR_3_UCSC.flanking.500_common+UTR_5_UCSC.flanking.500_common)>0,1,0))


res_all = res_all %>% mutate(MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_WB,MESuSiE_PIP_BB),
                             MESuSiE_PIP_Unique_all_mlk = pmax(MESuSiE_PIP_WB_all_mlk,MESuSiE_PIP_BB_all_mlk),
                             MESuSiE_PIP_Unique_all_glmnet = pmax(MESuSiE_PIP_WB_all_glmnet,MESuSiE_PIP_BB_all_glmnet))
  
custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 7.5),
    axis.text.y = element_text(size = 7.5),
    axis.title.x = element_text(size = 10, face="bold"),
    axis.title.y = element_text(size = 10, face="bold"),
    strip.text.x     = element_text(size = 12, face="bold"),   # ← bump this
    strip.text.y     = element_text(size = 12, face="bold"),   # ← and/or this
    strip.background = element_blank(),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.tag = element_text(size = 12),
    axis.line = element_line(color = "black")
  )
}

# Define annotation names and thresholds
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common", "CADD")



thresholds <- seq(0.5, 0.9, by = 0.05)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

### Compute enrichment for ancestry-shared variants
enrichment_results_shared <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% filter(MESuSiE_PIP_Shared_all_mlk > thresh) %>% nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  
  mes_cadd_mean <- mes_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  mlk_cadd_mean <- mlk_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  ## GLM method
  glm_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared_all_glmnet)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared_all_glmnet)
  
  glm_bin_signal <- glm_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  glm_bin_enrich <- glm_bin_signal / bg_bin
  glm_cadd_mean <- glm_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  glm_gtex_mean <- glm_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  glm_long <- bind_rows(
    glm_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    glm_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    glm_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "GLM", threshold = thresh)
  
  ## Paintor method
  paintor_subset <- res_all %>% 
    arrange(desc(Paintor_PIP)) %>% 
    top_n(n = top_N_signal, wt = Paintor_PIP)
  
  paintor_bin_signal <- paintor_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  paintor_bin_enrich <- paintor_bin_signal / bg_bin
  paintor_cadd_mean <- paintor_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  paintor_gtex_mean <- paintor_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  paintor_long <- bind_rows(
    paintor_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    paintor_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    paintor_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "Paintor", threshold = thresh)
  
  ## Paintor_all method
  paintor_all_subset <- res_all %>% 
    arrange(desc(Paintor_PIP_all)) %>% 
    top_n(n = top_N_signal, wt = Paintor_PIP_all)
  
  paintor_all_bin_signal <- paintor_all_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  paintor_all_bin_enrich <- paintor_all_bin_signal / bg_bin
  paintor_all_cadd_mean <- paintor_all_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  paintor_all_gtex_mean <- paintor_all_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  paintor_all_long <- bind_rows(
    paintor_all_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    paintor_all_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    paintor_all_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "Paintor_all", threshold = thresh)
  
  bind_rows(mes_long, mlk_long, glm_long, paintor_long, paintor_all_long) %>%
    mutate(variant = "Shared")
})

### Compute enrichment for ancestry-specific (unique) variants (only for MESuSiE, GLM, MLK)
enrichment_results_unique <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% filter(MESuSiE_PIP_Unique_all_mlk > thresh) %>% nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  mes_cadd_mean <- mes_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  mlk_cadd_mean <- mlk_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  ## GLM method
  glm_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique_all_glmnet)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique_all_glmnet)
  
  glm_bin_signal <- glm_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  glm_bin_enrich <- glm_bin_signal / bg_bin
  glm_cadd_mean <- glm_subset %>% 
    summarise(CADD = mean(CADD, na.rm = TRUE))
  glm_gtex_mean <- glm_subset %>% 
    summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  glm_long <- bind_rows(
    glm_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    glm_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    glm_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "GLM", threshold = thresh)
  
  bind_rows(mes_long, mlk_long, glm_long) %>% 
    mutate(variant = "Unique")
}) %>%
  # Ensure Unique results are only for MESuSiE, GLM, and MLK
  filter(method %in% c("MESuSiE", "GLM", "MLK"))

### Combine shared and unique results
enrichment_results_combined <- bind_rows(enrichment_results_shared, enrichment_results_unique) %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor", "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common","CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "CADD" = "CADD"))

# Set factor levels for method and variant (so that linetypes are applied correctly)
enrichment_results_combined$method <- factor(enrichment_results_combined$method, 
                                             levels = c("MESuSiE", "Paintor", "Paintor_all", "GLM", "MLK"))
enrichment_results_combined$variant <- factor(enrichment_results_combined$variant, levels = c("Shared", "Unique"))

### Create the combined plot with linetypes
combined_plot <- ggplot(enrichment_results_combined, aes(x = threshold, y = enrichment, color = method, linetype = variant)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Enrichment Results Across PIP Cutoffs for BMI",
    x = "Threshold",
    y = "Enrichment (Signal / Background)"
  ) +
  scale_color_manual(values = c(
    "MESuSiE"     = "#0072B2",   # Blue
    "Paintor"     = "#F0E442",   # Yellow
    "Paintor_all" = "#CC79A7",   # Purple
    "GLM"         = "#E69F00",   # Orange-Red
    "MLK"         = "#009E73"    # Green
  )) +
  scale_linetype_manual(values = c("Shared" = "solid", "Unique" = "dashed")) +
  theme_bw()

combined_plot

### Save the combined plot
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/"
ggsave(paste0(plot_dir, "enrichment_results_plot_bmi_tidy_mesusiemlk_cadd_Shared_Unique.png"),
       combined_plot, width = 225, height = 150, dpi = 500, units = 'mm')
































##########################################################


####Combined Shared vs. Unique plot PIP 0.5 to 0.9


##########################################################




library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)     # for map_df
library(patchwork) # for combining plots

## Preprocess res_all: compute unique PIP columns
res_all <- res_all %>% 
  mutate(MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_WB, MESuSiE_PIP_BB),
         MESuSiE_PIP_Unique_all_mlk = pmax(MESuSiE_PIP_WB_all_mlk, MESuSiE_PIP_BB_all_mlk),
         MESuSiE_PIP_Unique_all_glmnet = pmax(MESuSiE_PIP_WB_all_glmnet, MESuSiE_PIP_BB_all_glmnet))

# Define annotation names and thresholds
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common","CADD")

thresholds <- seq(0.5, 0.9, by = 0.05)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

### Compute enrichment for ancestry-shared variants (only MESuSiE and MLK)
enrichment_results_shared <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Shared_all_mlk > thresh) %>% 
    nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  
  mes_cadd_mean <- mes_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  
  mlk_cadd_mean <- mlk_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  bind_rows(mes_long, mlk_long) %>% 
    mutate(variant = "Shared",
           top_N_signal = top_N_signal)
})

### Compute enrichment for ancestry-unique variants (only MESuSiE and MLK)
enrichment_results_unique <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Unique_all_mlk > thresh) %>% 
    nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  
  mes_cadd_mean <- mes_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  
  mlk_cadd_mean <- mlk_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  bind_rows(mes_long, mlk_long) %>% 
    mutate(variant = "Unique",
           top_N_signal = top_N_signal)
}) 

### Combine Shared and Unique results and recode annotations
enrichment_results_combined <- bind_rows(enrichment_results_shared, enrichment_results_unique) %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor",
                                        "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common",
                                        "CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "CADD" = "CADD"))

# Update factor levels for method and variant (only MESuSiE and MLK)
enrichment_results_combined$method <- factor(enrichment_results_combined$method, 
                                             levels = c("MESuSiE", "MLK"))
enrichment_results_combined$variant <- factor(enrichment_results_combined$variant, 
                                              levels = c("Shared", "Unique"))

### Define x-axis labeling for selected thresholds only: 0.5, 0.6, 0.7, 0.8, and 0.9.
selected_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
top_signals_shared <- enrichment_results_combined %>%
  filter(threshold %in% selected_thresholds) %>%
  group_by(variant,threshold) %>%
  summarise(top_signal = unique(top_N_signal)) %>%
  arrange(threshold) %>% filter(variant == "Shared")

top_signals_unique <- enrichment_results_combined %>%
  filter(threshold %in% selected_thresholds) %>%
  group_by(variant,threshold) %>%
  summarise(top_signal = unique(top_N_signal)) %>%
  arrange(threshold) %>% filter(variant == "Unique")

x_breaks <- selected_thresholds
x_labels <- paste0(x_breaks, " \n(n=", top_signals_shared$top_signal, ")","\n(n=", top_signals_unique$top_signal, ")")

### Split the data into two panels:
panelA <- enrichment_results_combined %>% 
  filter(!(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP","brain_ind_eQTL"))) %>% dplyr::rename(PIP = variant)
panelB <- enrichment_results_combined %>% 
  filter(annotation %in% c("CADD", "brain_ind_eQTL")) %>% dplyr::rename(PIP = variant)

### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations
plotA <- ggplot(panelA, aes(x = threshold, y = enrichment, color = method, linetype = PIP)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel A: Enrichment Results (Signal/Background) for BMI",
       x = "PIP Threshold",
       y = "Enrichment (Signal / Background)",
       color    = "Method",   # legend title for color
       linetype = "Variant"   # legend title for linetype
  ) +
  scale_color_manual(name   = "Method",values = c(
    "MESuSiE" = "#009E73",   # Green
    "MLK"     = "#0072B2"    # Blue
  )) +
  scale_linetype_manual(name   = "Variant",values = c("Shared" = "solid", "Unique" = "dashed")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw()+custom_theme()+  theme(plot.title = element_blank(),legend.position  = "bottom")


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_enrichment_plot_bmi_updated_tidy_mesusiemlk_cadd_Shared_Unique.pdf"),
       plotA, width = 6, height = 3.9, dpi = 500, units = 'in')


### Create Panel B plot: Average Scores for continuous annotations (CADD and GTEx_eQTL_MaxCPP)
plotB <- ggplot(panelB, aes(x = threshold, y = enrichment, color = method, linetype = PIP)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel B: Average Scores for BMI",
       x = "PIP Threshold",
       y = "Average Score") +
  scale_color_manual(values = c(
    "MESuSiE" = "#009E73",   # Green
    "MLK"     = "#0072B2"    # Blue
  )) +
  scale_linetype_manual(values = c("Shared" = "solid", "Unique" = "dashed")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw()+custom_theme()+  theme(plot.title = element_blank(),legend.position  = "bottom")

### Combine the two panels vertically using patchwork
combined_plot <- plotA / plotB

# Display the combined plot
combined_plot


# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_enrichment_and_score_plot_bmi_updated_tidy_mesusiemlk_cadd_Shared_Unique.png"),
       combined_plot, width = 275, height = 300, dpi = 500, units = 'mm')

plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_eqtl_cadd_annotations_plot_bmi_mesusiemlk_Shared_Unique.pdf"),
       plotB, width = 6, height = 3.9, dpi = 500, units = 'in')




















##########################################################


####Combined Shared vs. Unique plot at PIP 0.5 and 0.9


##########################################################




library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)     # for map_df
library(patchwork) # for combining plots

## Preprocess res_all: compute unique PIP columns
res_all <- res_all %>% 
  mutate(MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_WB, MESuSiE_PIP_BB),
         MESuSiE_PIP_Unique_all_mlk = pmax(MESuSiE_PIP_WB_all_mlk, MESuSiE_PIP_BB_all_mlk),
         MESuSiE_PIP_Unique_all_glmnet = pmax(MESuSiE_PIP_WB_all_glmnet, MESuSiE_PIP_BB_all_glmnet))

# Define annotation names and thresholds
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common", "CADD")

thresholds <- c(0.5, 0.9)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

### Compute enrichment for ancestry-shared variants (only MESuSiE and MLK)
enrichment_results_shared <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Shared_all_mlk > thresh) %>% 
    nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  
  mes_cadd_mean <- mes_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Shared_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Shared_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  
  mlk_cadd_mean <- mlk_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  bind_rows(mes_long, mlk_long) %>% 
    mutate(variant = "Shared",
           top_N_signal = top_N_signal)
})

### Compute enrichment for ancestry-unique variants (only MESuSiE and MLK)
enrichment_results_unique <- map_df(thresholds, function(thresh) {
  
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Unique_all_mlk > thresh) %>% 
    nrow()
  
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n())))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique)
  
  mes_bin_signal <- mes_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mes_bin_enrich <- mes_bin_signal / bg_bin
  
  mes_cadd_mean <- mes_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mes_gtex_mean <- mes_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mes_long <- bind_rows(
    mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MESuSiE", threshold = thresh)
  
  ## MLK method
  mlk_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Unique_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Unique_all_mlk)
  
  mlk_bin_signal <- mlk_subset %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
  mlk_bin_enrich <- mlk_bin_signal / bg_bin
  
  mlk_cadd_mean <- mlk_subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
  mlk_gtex_mean <- mlk_subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
  
  mlk_long <- bind_rows(
    mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
    mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
  ) %>% mutate(method = "MLK", threshold = thresh)
  
  bind_rows(mes_long, mlk_long) %>% 
    mutate(variant = "Unique",
           top_N_signal = top_N_signal)
}) 

### Combine Shared and Unique results and recode annotations
enrichment_results_combined <- bind_rows(enrichment_results_shared, enrichment_results_unique) %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor",
                                        "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common",
                                        "brain_ind_eQTL",
                                        "CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "brain_ind_eQTL" = "brain_ind_eQTL",
                             "CADD" = "CADD"))

# Update factor levels for method and variant (only MESuSiE and MLK)
enrichment_results_combined$method <- factor(
  recode(enrichment_results_combined$method, "MLK" = "INFUSE"),
  levels = c("INFUSE","MESuSiE")
)

enrichment_results_combined$variant <- factor(enrichment_results_combined$variant, 
                                              levels = c("Shared", "Unique"))

selected_thresholds <- c(0.5)
top_signals_shared <- enrichment_results_combined %>%
  filter(threshold %in% selected_thresholds) %>%
  group_by(variant,threshold) %>%
  summarise(top_signal = unique(top_N_signal)) %>%
  arrange(threshold) %>% filter(variant == "Shared")

top_signals_unique <- enrichment_results_combined %>%
  filter(threshold %in% selected_thresholds) %>%
  group_by(variant,threshold) %>%
  summarise(top_signal = unique(top_N_signal)) %>%
  arrange(threshold) %>% filter(variant == "Unique")

x_breaks <- selected_thresholds
x_labels <- paste0(x_breaks, " \n(n=", top_signals_shared$top_signal, ")","\n(n=", top_signals_unique$top_signal, ")")

### Split the data into two panels:
panelA <- enrichment_results_combined %>% 
  filter(!(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP"))) %>% dplyr::rename(PIP = variant) %>% filter(threshold == 0.5)

panelA_upper <- enrichment_results_combined %>% 
  filter(!(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP"))) %>% dplyr::rename(PIP = variant) %>% filter(threshold == 0.9)


panelB <- enrichment_results_combined %>% 
  filter(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP")) %>% dplyr::rename(PIP = variant) %>% filter(threshold == 0.5)

### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.5
plotA <- ggplot(panelA, 
                aes(x = threshold, 
                    y = enrichment,
                    linetype = PIP,
                    fill = method 
                    )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at PIP = 0.5 for BMI",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE" = "#009E73",
    "INFUSE"   = "#0072B2"
  )) +
  scale_linetype_manual(values = c(
    "Shared" = "solid",
    "Unique" = "dashed"
  )) +
  theme_bw() +
  custom_theme() +
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
    )


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster_ashg_2025/"
ggsave(paste0(plot_dir, "Figure_5_combined_enrichment_plot_bp_updated_tidy_mesusiemlk_cadd_Shared_Unique.png"),
       plotA, width = 10.7, height = 5.92, dpi = 500, units = 'in')


### Create Panel B plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.9
plotA_upper <- ggplot(panelA_upper, 
                aes(x = threshold, 
                    y = enrichment,
                    linetype = PIP,
                    fill = method 
                )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at PIP = 0.9 for BMI",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE" = "#009E73",
    "INFUSE"   = "#0072B2"
  )) +
  scale_linetype_manual(values = c(
    "Shared" = "solid",
    "Unique" = "dashed"
  )) +
  theme_bw() +
  custom_theme() +
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
  )


### Combine the two panels vertically using patchwork
combined_plot <- (plotA / plotA_upper)+plot_annotation(tag_levels = "a") +  custom_theme()+plot_layout(guides = 'collect') &
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
  )


# Display the combined plot
combined_plot


# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"

ggsave(paste0(plot_dir, "combined_enrichment_and_score_plot_updated_tidy_mesusiemlk_cadd_Shared_Unique_figure5.pdf"),
       combined_plot, width = 7, height = 6.5, dpi = 500, units = 'in')
