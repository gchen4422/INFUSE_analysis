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
library(purrr)


res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"
load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated_evo2_rerun_test.RData")) ### use this as most recent BP results
res_all = res_all_new
rm(res_all_new)
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
res_all<-res_all%>%mutate(utr_comb_flanking = ifelse((UTR_3_UCSC.flanking.500_common+UTR_5_UCSC.flanking.500_common)>0,1,0))
#load("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster/BMI_setsize_enrich.Rdata")




library(dplyr)

# 1. Calculate the PIP gain for every SNP
res_comparison <- res_all %>%
  filter(Trait %in% c("SBP", "DBP", "BMI")) %>%
  mutate(PIP_Gain = PIP_Either_Ensemble_infuse - MESuSiE_PIP_Either_all_mlk)

# 2. Find the top SNP per Region with the highest gain
top_gains_per_region <- res_comparison %>%
  group_by(Region) %>%
  slice_max(order_by = PIP_Gain, n = 1) %>%
  ungroup() %>%
  # Filter for SNPs that actually became "Strong" (e.g., Ensemble PIP > 0.5)
  filter(PIP_Either_Ensemble_infuse > 0.5) %>%
  arrange(desc(PIP_Gain))

# 3. Display the Top 5 "Winner" Loci
winner_loci <- head(top_gains_per_region, 5)
print(winner_loci %>% select(Trait, Region, SNP, MESuSiE_PIP_Either_all_mlk, PIP_Either_Ensemble_infuse, PIP_Gain))


# Check if the winners are heart/artery eQTLs
winner_loci %>% 
  select(SNP, brain_ind_eQTL, PIP_Gain,evo2_delta_score)








ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", "heart_ind_eQTL","artery_ind_eQTL",
                  "GTEx_eQTL_MaxCPP_common", "CADD")




# ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
#                    "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
#                    "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL", "CADD")

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x     = element_text(size = 8, face="bold"),   # ← bump this
    strip.text.y     = element_text(size = 8, face="bold"),   # ← and/or this
    strip.background = element_blank(),
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    plot.title = element_text(size=8, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

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



# 1. Configuration
thresholds <- seq(0.5, 0.9, by = 0.05)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

# Define a mapping of method names to their corresponding PIP columns in res_all
method_map <- list(
  "INFUSE-"      = "MESuSiE_PIP_Either",
  "INFUSE"          = "MESuSiE_PIP_Either_all_mlk",
  #"GLM"          = "MESuSiE_PIP_Either_all_glmnet",
  #"Paintor"      = "Paintor_PIP",
  #"Paintor_all"  = "Paintor_PIP_all",
  "Evo2_Only"    = "PIP_Either_Evo2_Only_infuse", # New
  "Ensemble"     = "PIP_Either_Ensemble_infuse"   # New
)

# 2. Functional Enrichment Calculation
enrichment_results <- map_df(thresholds, function(thresh) {
  
  # Determine top_N based on the baseline MLK model at this threshold
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Either_all_mlk > thresh) %>% 
    nrow()
  
  if(top_N_signal == 0) return(NULL) # Skip if no SNPs meet threshold
  
  # Compute background frequency for binary annotations once per threshold
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))
  
  # Helper function to process each method to avoid repetitive code blocks
  process_method <- function(method_name, pip_col) {
    subset_df <- res_all %>% 
      arrange(desc(!!sym(pip_col))) %>% 
      slice_head(n = top_N_signal)
    
    # Binary Enrichment
    bin_enrich <- subset_df %>% 
      summarise(across(all_of(binary_ann), ~ (sum(.x, na.rm = TRUE) / n()))) / bg_bin
    
    # Continuous Means (CADD & GTEx)
    cont_means <- subset_df %>% 
      summarise(
        CADD = mean(CADD, na.rm = TRUE),
        GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE)
      )
    
    # Pivot to long format
    bind_rows(
      bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
      cont_means %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
    ) %>% 
      mutate(method = method_name, threshold = thresh)
  }
  
  # Apply the helper function across all methods in our map
  map2_df(names(method_map), method_map, process_method)
})

# 3. Post-Processing: Factor Levels and Recoding
annotation_levels <- c(
  "non_synonymous_common", "synonymous", "promotor", "utr_comb",
  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common","heart_ind_eQTL", "artery_ind_eQTL", 
  "GTEx_eQTL_MaxCPP_common", "CADD"
)

enrichment_results <- enrichment_results %>%
  mutate(
    annotation = factor(annotation, levels = annotation_levels),
    annotation = recode(annotation,
                        "non_synonymous_common"        = "Non-synonymous",
                        "synonymous"                   = "Synonymous",
                        "promotor"                     = "Promotor",
                        "utr_comb"                     = "UTR",
                        "Conserved_LindbladToh_common" = "Conserved",
                        "H3K4me3_Trynka_common"        = "H3K4me3",
                        "heart_ind_eQTL" = "Heart ind. eQTL",
                        "artery_ind_eQTL" = "Artery ind. eQTL",
                        "GTEx_eQTL_MaxCPP_common"      = "GTEx_eQTL_MaxCPP",
                        "CADD"                         = "CADD"
    ),
    method = factor(method, levels = c("INFUSE-", "INFUSE", "Evo2_Only", "Ensemble"))
  )

# Create the plot with updated, colorblind-friendly colors
enrichment_plot <- ggplot(
  enrichment_results,
  aes(x = threshold, y = enrichment, color = method, linetype = method)
) +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  geom_point(
    size = 2,
    alpha = 0.9,
    position = position_dodge(width = 0.01)
  ) +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Functional Annotation Patterns Across PIP Thresholds (BPs)",
    x = "PIP threshold",
    y = "Enrichment or mean score",
    color = "Method",
    linetype = "Method"
  ) +
  scale_color_manual(values = c(
    "INFUSE-"   = "#009E73",
    "INFUSE"    = "#56B4E9",
    "Evo2_Only" = "#CC79A7",
    "Ensemble"  = "#999999"
  )) +
  theme_bw() + custom_theme()+
  theme(
    legend.position = "bottom"
  )
# Display plot
print(enrichment_plot)

enrichment_plot
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/"
ggsave(paste0(plot_dir,"enrichment_results_plot_bp_updated_tidy_mesusiemlk_cadd_Either.pdf"),enrichment_plot,width=225,height =150,dpi=500,units='mm')





rank_thresholds <- c(10, 20, 50, 100, 200)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))
global_bg_bin <- res_all %>% 
  summarise(across(all_of(binary_ann), ~ mean(.x, na.rm = TRUE)))

# 2. Functional Enrichment Calculation
enrichment_ranks <- map_df(rank_thresholds, function(top_N) {
  
  process_method <- function(method_name, pip_col) {
    # Sort and slice
    subset_df <- res_all %>% 
      arrange(desc(!!sym(pip_col))) %>% 
      slice_head(n = top_N)
    
    # Binary Enrichment (Signal frequency / Global background frequency)
    bin_enrich <- subset_df %>% 
      summarise(across(all_of(binary_ann), ~ mean(.x, na.rm = TRUE))) / global_bg_bin
    
    # Continuous Means (Actual means, not enrichment ratios)
    cont_means <- subset_df %>% 
      summarise(
        CADD = mean(CADD, na.rm = TRUE),
        GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE)
      )
    
    bind_rows(
      bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment"),
      cont_means %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment")
    ) %>% 
      mutate(method = method_name, rank_cutoff = top_N)
  }
  
  map2_df(names(method_map), method_map, process_method)
})

enrichment_results <- enrichment_ranks %>%
  mutate(
    # Ensure rank_cutoff is treated as a factor so the x-axis isn't stretched linearly
    rank_cutoff = factor(rank_cutoff, levels = c(10, 20, 50, 100, 200)),
    annotation = factor(annotation, levels = annotation_levels),
    annotation = recode(annotation,
                        "non_synonymous_common"        = "Non-synonymous",
                        "synonymous"                   = "Synonymous",
                        "promotor"                     = "Promotor",
                        "utr_comb"                     = "UTR",
                        "Conserved_LindbladToh_common" = "Conserved",
                        "H3K4me3_Trynka_common"        = "H3K4me3",
                        "GTEx_eQTL_MaxCPP_common"      = "GTEx_eQTL_MaxCPP",
                        "CADD"                         = "CADD"
    ),
    method = factor(method, levels = c("MESuSiE", "Paintor", "Paintor_all", "GLM", "MLK", "Evo2_Only", "Ensemble"))
  )

# 2. Create the Rank-Based Plot
enrichment_plot <- ggplot(enrichment_results, aes(x = rank_cutoff, y = enrichment, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  # Use a log scale for Y if some annotations (like Non-synonymous) have huge enrichment
  # scale_y_log10() + 
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Functional Enrichment by SNP Rank (BP Traits)",
    subtitle = "Comparing Top N SNPs across different fine-mapping methods",
    x = "Top N SNPs (Ranked by PIP)",
    y = "Enrichment (Signal / Background)",
    color = "Method"
  ) +
  scale_color_manual(values = c(
    "MESuSiE"     = "#56B4E9",   
    "Paintor"     = "#999999",   
    "Paintor_all" = "#CC79A7",   
    "GLM"         = "#E69F00",   
    "MLK"         = "#009E73",   
    "Evo2_Only"   = "#000000",   
    "Ensemble"    = "#D55E00"    
  )) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate labels for readability
  )

# Display plot
print(enrichment_plot)









######################################


#Barplot at PIP threshold 0.5


######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Define annotation column names
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL","artery_ind_eQTL","CADD")

# Set the threshold to 0.5
thresh <- 0.5

# Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

# Define top_N_signal based on MESuSiE_PIP_Either_all_mlk exceeding the threshold
top_N_signal <- res_all %>% 
  filter(MESuSiE_PIP_Either_all_mlk > thresh) %>% 
  nrow()

# Compute background frequency for binary annotations (assumes n() - top_N_signal > 0)
bg_bin <- res_all %>% 
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))

## ---- Compute metrics for each method (excluding GLM) ----

# MESuSiE method
mes_subset <- res_all %>% 
  arrange(desc(MESuSiE_PIP_Either)) %>% 
  top_n(n = top_N_signal, wt = MESuSiE_PIP_Either)

mes_bin_signal <- mes_subset %>% 
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
mes_bin_enrich <- mes_bin_signal / bg_bin

mes_cadd_mean <- mes_subset %>% 
  summarise(CADD = mean(CADD, na.rm = TRUE))
mes_gtex_mean <- mes_subset %>% 
  summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))

mes_long <- bind_rows(
  mes_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  mes_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  mes_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score")
) %>% mutate(method = "MESuSiE", threshold = thresh)

# MLK method
mlk_subset <- res_all %>% 
  arrange(desc(MESuSiE_PIP_Either_all_mlk)) %>% 
  top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)

mlk_bin_signal <- mlk_subset %>% 
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
mlk_bin_enrich <- mlk_bin_signal / bg_bin

mlk_cadd_mean <- mlk_subset %>% 
  summarise(CADD = mean(CADD, na.rm = TRUE))
mlk_gtex_mean <- mlk_subset %>% 
  summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))

mlk_long <- bind_rows(
  mlk_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  mlk_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  mlk_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score")
) %>% mutate(method = "MLK", threshold = thresh)

# Paintor method
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
  paintor_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  paintor_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  paintor_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score")
) %>% mutate(method = "Paintor", threshold = thresh)

# Paintor_all method
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
  paintor_all_bin_enrich %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  paintor_all_cadd_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score"),
  paintor_all_gtex_mean %>% pivot_longer(cols = everything(), names_to = "annotation", values_to = "score")
) %>% mutate(method = "Paintor_all", threshold = thresh)

# Combine methods (GLM removed)
enrichment_results <- bind_rows(mes_long, mlk_long, paintor_long, paintor_all_long)

# Reorder and recode the 'annotation' factor
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor",
                                        "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL",
                                        "CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "heart_ind_eQTL" = "heart_ind_eQTL",
                             "artery_ind_eQTL"= "artery_ind_eQTL",
                             "CADD" = "CADD")) %>%
  mutate(method = recode(method,
                             "MESuSiE" = "MESuSiE",
                             "MLK" = "INFUSE",
                             "Paintor" = "Paintor",
                             "Paintor_all" = "Paintor+annot"))

enrichment_results$method <- factor(enrichment_results$method, 
                                    levels = c("Paintor", "Paintor+annot","MESuSiE","INFUSE"))

# Split the data into two panels:
# Panel A: All annotations except "CADD" and "GTEx_eQTL_MaxCPP"
enrichment_panelA <- enrichment_results %>% 
  filter(!(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP","heart_ind_eQTL","artery_ind_eQTL")))

# Panel B: Only "CADD" and "GTEx_eQTL_MaxCPP"
enrichment_panelB <- enrichment_results %>% 
  filter(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP"))

# Define fill colors
fill_colors <- c(
  "MESuSiE"      = "#84c3b7",
  "Paintor"      = "#f2b56e",
  "Paintor+annot"  = "#f57c6e",
  "INFUSE"          = "#71b8ed"
)

# Create Panel A plot: Enrichment (Signal / Background)
plotA <- ggplot(enrichment_panelA, aes(x = method, y = score, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Functional Enrichment at PIP = 0.5 for BP",
    x = "Method",
    y = "Enrichment (Signal / Background)"
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_bw()+custom_theme()+ theme(legend.position = "bottom")

# Create Panel B plot: Average Score
plotB <- ggplot(enrichment_panelB, aes(x = method, y = score, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Panel B: Average Scores at Threshold = 0.5 for BP",
    x = "Method",
    y = "Average Score"
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_bw()

# Combine the two panels vertically using patchwork
combined_plot <- plotA / plotB

# Display the combined plot
combined_plot

# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
system(paste0("mkdir -p ",plot_dir))
ggsave(paste0(plot_dir, "combined_enrichment_and_score_plot_bp_updated_tidy_mesusiemlk_cadd_Either_bar_noGLM.png"),
       combined_plot, width = 300, height = 300, dpi = 500, units = 'mm')



######################################


#For PIP threshold from 0.5 to 0.9


######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)     # for map_df
library(patchwork) # for combining plots

# Define annotation column names
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL","artery_ind_eQTL","CADD")

# Define thresholds from 0.5 to 0.9 (by 0.05)
thresholds <- seq(0.5, 0.9, by = 0.05)

# Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

# Loop over thresholds and compute enrichment results for each method
enrichment_results <- map_df(thresholds, function(thresh) {
  
  # Compute top_N_signal based on MESuSiE_PIP_Either_all_mlk exceeding the current threshold
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Either_all_mlk > thresh) %>% 
    nrow()
  
  # Compute background frequency for binary annotations
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Either)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Either)
  
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
    arrange(desc(MESuSiE_PIP_Either_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)
  
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
  
  # Combine methods for this threshold (GLM removed) and add the top signal number as a column  
  bind_rows(mes_long, mlk_long, paintor_long, paintor_all_long) %>% 
    mutate(top_N_signal = top_N_signal)
})

# Reorder and recode the annotation factor
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor",
                                        "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL",
                                        "CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "heart_ind_eQTL" = "heart_ind_eQTL",
                             "artery_ind_eQTL"= "artery_ind_eQTL",
                             "CADD" = "CADD"))

# Update method factor (without GLM)
enrichment_results$method <- factor(enrichment_results$method, 
                                    levels = c( "Paintor", "Paintor_all", "MESuSiE","MLK"))

# Split the data into two panels:
# Panel A: All annotations except CADD and GTEx_eQTL_MaxCPP
panelA <- enrichment_results %>% 
  filter((annotation %in% c("heart_ind_eQTL", "artery_ind_eQTL")))
# Panel B: Only CADD and GTEx_eQTL_MaxCPP (continuous scores)
panelB <- enrichment_results %>% 
  filter(annotation %in% c("CADD"))

panelC <- enrichment_results %>% 
  filter((!annotation %in% c("heart_ind_eQTL","artery_ind_eQTL","CADD","GTEx_eQTL_MaxCPP")))

# Define the thresholds you want for x-axis labeling
selected_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)

# Extract top_N_signal for only these thresholds
top_signals <- enrichment_results %>%
  filter(threshold %in% selected_thresholds) %>%
  group_by(threshold) %>%
  summarise(top_signal = unique(top_N_signal)) %>%
  arrange(threshold)

# Define the x-axis breaks and custom labels
x_breaks <- selected_thresholds
x_labels <- paste0(x_breaks, " \n(n=", top_signals$top_signal, ")")

# Define a color palette for methods
color_palette <- c(  "MESuSiE"      = "#84c3b7",
                     "Paintor"      = "#f2b56e",
                     "Paintor_all"  = "#f57c6e",
                     "MLK"          = "#71b8ed")

# Create Panel A plot: Enrichment (Signal / Background)
panelA = panelA %>% mutate(annotation = recode(annotation,
                                               "heart_ind_eQTL" = "Heart ind. eQTL"  ,
                                               "artery_ind_eQTL"= "Artery ind. eQTL"))

plotA <- ggplot(panelA, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    x = "PIP Threshold",
    y = "Enrichment (Signal / Background)",
    title = "Functional Enrichment at PIP = 0.5 for BP",
    color    = "Method"
  ) +
  scale_color_manual(name = "Method",values = color_palette) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() +
  custom_theme() +
  theme(plot.title = element_blank())

# Create Panel B plot: Average Scores
plotB <- ggplot(panelB, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    x = "PIP Threshold",
    y = "Average Score",
    color    = "Method"
  ) +
  scale_color_manual(name = "Method",values = color_palette) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() +
  custom_theme() +
  theme(plot.title = element_blank())
# Combine the two panels vertically using patchwork
combined_plot <- (
  plotA | plotB
) +
  plot_layout(
    guides  = "collect",    # <-- collect both legends into one
    widths  = c(2, 1)
  ) &
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal"
  )

plotC <- ggplot(panelC, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    x = "PIP Threshold",
    y = "Average Score",
    color    = "Method"
  ) +
  scale_color_manual(name = "Method",values = color_palette) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() +
  custom_theme() +
  theme(
    plot.title        = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal"
  )

# Display the combined plot
combined_plot
# Save the combined plot to the designated directory
#plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster/"
ggsave(paste0(plot_dir, "combined_eqtl_enrichment_and_score_plot_bp_tidy_mesusiemlk_cadd_Either_noGLM.pdf"),
       combined_plot, width = 6, height = 4, dpi = 500, units = 'in')

ggsave(paste0(plot_dir, "combined_6_annotations_plot_bp_mesusiemlk_Either.pdf"),
       plotC, width = 5, height = 4, dpi = 500, units = 'in')





###########################


# Size and Z-score plot


###########################



color_palette <- c(  "MESuSiE"      = "#84c3b7",
                     "Paintor"      = "#f2b56e",
                     "Paintor_all"  = "#f57c6e",
                     "MLK"          = "#71b8ed")



### Median set size by Trait (across all locus)
all_sets_info <- data.frame(
  res_all %>% 
    group_by(Trait, Region) %>% 
    summarise(across(c("Paintor_cs", 
                       "Paintor_cs_all",
                       "MESuSiE_cs", 
                       "MESuSiE_res_mlk_annot_all_cs"), 
                     ~ sum(.x, na.rm = TRUE)))
) %>% 
  filter(MESuSiE_cs != 0, 
         MESuSiE_res_mlk_annot_all_cs != 0, 
         Paintor_cs != 0, 
         Paintor_cs_all != 0)  # Remove loci with zero set sizes for any method

# Reshape to long format for plotting
all_sets_info_long <- all_sets_info %>% 
  pivot_longer(cols = -c(Trait, Region), names_to = "Method", values_to = "Count")

# Define factor levels and recode the method names
all_sets_info_long$Method <- factor(all_sets_info_long$Method, 
                                    levels = c("Paintor_cs", "Paintor_cs_all","MESuSiE_cs", "MESuSiE_res_mlk_annot_all_cs"))

levels(all_sets_info_long$Method) <- c("Paintor", "Paintor+annot","MESuSiE","INFUSE")

# Define the new color palette
color_palette <- c(
  "MESuSiE"      = "#84c3b7",
  "Paintor"      = "#f2b56e",
  "Paintor+annot"  = "#f57c6e",
  "INFUSE"          = "#71b8ed"
)

# Create boxplot for set size
p_set <- ggplot(data = all_sets_info_long, aes(x = Trait, y = log2(Count+1), fill = Method)) +
  geom_boxplot(outlier.size = 0.1, fatten = 0.5, color = "darkgray") +
  scale_fill_manual(values = color_palette, guide = FALSE) +
  theme_bw() + 
  labs(
    title = "Credible Set Size",
    x     = "Trait",
    y     = "log2(Set Size + 1)"
  )+
  #coord_cartesian(ylim = c(0, 300)) +
  custom_theme()


# For each method, compute the median of the maximum absolute Z-score (from either zscore_WB or zscore_BB)
MESuSiE_cs_Z <- res_all %>% 
  group_by(Trait) %>% 
  filter(MESuSiE_cs == 1) %>% 
  summarise(zmax = median(pmax(abs(zscore_WB), abs(zscore_BB))))

MESuSiE_mlk_cs_Z <- res_all %>% 
  group_by(Trait) %>% 
  filter(MESuSiE_res_mlk_annot_all_cs == 1) %>% 
  summarise(zmax = median(pmax(abs(zscore_WB), abs(zscore_BB)))) %>% 
  pull(zmax)

Paintor_cs_Z <- res_all %>% 
  group_by(Trait) %>% 
  filter(Paintor_cs == 1) %>% 
  summarise(zmax = median(pmax(abs(zscore_WB), abs(zscore_BB)))) %>% 
  pull(zmax)

Paintor_all_cs_Z <- res_all %>% 
  group_by(Trait) %>% 
  filter(Paintor_cs_all == 1) %>% 
  summarise(zmax = median(pmax(abs(zscore_WB), abs(zscore_BB)))) %>% 
  pull(zmax)

# Combine into a single data frame.
# Note: MESuSiE_cs_Z already contains the Trait and its corresponding zmax.
set_size_z_info <- data.frame(
  Trait         = MESuSiE_cs_Z$Trait,
  MESuSiE       = MESuSiE_cs_Z$zmax,
  INFUSE         = MESuSiE_mlk_cs_Z,       # or MESuSiE_mlk_cs_Z$zmax if that’s a data.frame
  Paintor       = Paintor_cs_Z,           # same idea here
  `Paintor+annot` = Paintor_all_cs_Z,     # backticks let you use the plus
  check.names   = FALSE                   # ← important to keep the + in your colname
)
# Reshape the data into long format for plotting
set_size_z_info_long <- set_size_z_info %>% 
  pivot_longer(cols = -Trait, names_to = "Method", values_to = "Z") %>% 
  mutate(Method = factor(Method, levels = c( "Paintor", "Paintor+annot","MESuSiE", "INFUSE")))

# Create barplot for median |Z| values
p_z <- ggplot(data = set_size_z_info_long, aes(x = Trait, y = Z, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = color_palette) +
  geom_text(aes(label = round(Z, 2)), 
            position = position_dodge(width = 1), 
            vjust = -0.5, 
            size = 5 * 5/10) +
  theme_bw() +
  labs(
    title = "Median |Z|",
    x     = "Trait",
    y     = "Median |Z|"
  ) +
  ylim(0, max(round(set_size_z_info_long$Z, 2) + 1)) +
  custom_theme()

# (Optional) To display the plots:
print(p_set)
print(p_z)

p_out<-(p_set+ labs(title = NULL))/(p_z+labs(title = NULL))+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom')
p_out

plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "set_size_bp_updated.png"),
       p_out, width = 200, height = 300, dpi = 500, units = 'mm')


###plotA is barplot 0.5 with first 6 annotations

drop_guides <- guides(
  color = "none",
  shape = "none",
  plot.title  = "none",
  linetype = "none"
)

p_total <- 
  (
    p_out | 
      (plotA + labs(title = NULL))    # remove plotA’s title here
  ) +
  plot_layout(
    guides  = "collect",   # collects both legends into one
    heights = c(0.3, 1),
    widths  = c(1, 3)
  ) +drop_guides+
  plot_annotation(
    title      = "Multi-ancestry fine-mapping on All of Us BP data",
    tag_levels = "a"       # will tag panels “A” and “B”
  ) &
  custom_theme() &         # apply your custom theme everywhere
  theme(
    plot.tag           = element_text(size = 12),
    legend.position    = "none",
    legend.justification = "center",
    legend.box         = "horizontal",
    legend.box.just    = "center"
  )

print(p_total)

p_total <- wrap_elements(full = p_total)


combined_plot = (p_total_bmi/p_total)&custom_theme()
combined_plot

plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster_ashg_2025/"
save(combined_plot, file = paste0(plot_dir,"Figure_4_combined_bmi_bp.Rdata"))


ggsave(paste0(plot_dir, "set_size_enrichment_mlk_50_bp_modified_log2.pdf"),
       p_total, width = 8.2, height = 3.9, dpi = 500, units = 'in')

ggsave(paste0(plot_dir, "set_size_enrichment_mlk_50_combined_modified_log2.pdf"),
       combined_plot, width = 8.2, height = 6.8, dpi = 500, units = 'in')

ggsave(paste0(plot_dir, "Figure_4_set_size_enrichment_mlk_50_combined_modified_log2.png"),
       combined_plot, width = 14.07, height = 13.31, dpi = 500, units = 'in')





###########################


# Seperate panel for eQTL


###########################





library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)     # for map_df
library(patchwork) # for combining plots

# Define annotation column names
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", "utr_comb",
                  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
                  "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL","artery_ind_eQTL","CADD")

# Define thresholds from 0.5 to 0.9
thresholds <- seq(0.5, 0.9, by = 0.05)

# Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

enrichment_results <- map_df(thresholds, function(thresh) {
  
  # Define top_N_signal based on MESuSiE_PIP_Either_all_mlk exceeding the current threshold
  top_N_signal <- res_all %>% 
    filter(MESuSiE_PIP_Either_all_mlk > thresh) %>% 
    nrow()
  
  # Compute background frequency for binary annotations (assumes n() - top_N_signal > 0)
  bg_bin <- res_all %>% 
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))
  
  ## MESuSiE method
  mes_subset <- res_all %>% 
    arrange(desc(MESuSiE_PIP_Either)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Either)
  
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
    arrange(desc(MESuSiE_PIP_Either_all_mlk)) %>% 
    top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)
  
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
  
  # Combine methods for this threshold (GLM removed)
  bind_rows(mes_long, mlk_long, paintor_long, paintor_all_long)
})

# Reorder and recode the annotation factor
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor",
                                        "utr_comb",
                                        "Conserved_LindbladToh_common",
                                        "H3K4me3_Trynka_common",
                                        "GTEx_eQTL_MaxCPP_common",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL",
                                        "CADD"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "Conserved_LindbladToh_common" = "Conserved",
                             "H3K4me3_Trynka_common" = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common" = "GTEx_eQTL_MaxCPP",
                             "heart_ind_eQTL" = "heart_ind_eQTL",
                             "artery_ind_eQTL"= "artery_ind_eQTL",
                             "CADD" = "CADD"))

# Update method factor (without GLM)
enrichment_results$method <- factor(enrichment_results$method, 
                                    levels = c("MESuSiE", "Paintor", "Paintor_all", "MLK"))

# Split the data into three panels:
# Panel A: All annotations except "CADD", "GTEx_eQTL_MaxCPP", and excluding those for Panel C ("heart_ind_eQTL", "artery_ind_eQTL")
panelA <- enrichment_results %>% 
  filter(!(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP", "heart_ind_eQTL","artery_ind_eQTL")))
# Panel B: Only the continuous annotations "CADD" and "GTEx_eQTL_MaxCPP"
panelB <- enrichment_results %>% 
  filter(annotation %in% c("CADD", "GTEx_eQTL_MaxCPP"))
# Panel C: Only the annotations "heart_ind_eQTL" and "H3K4me3"
panelC <- enrichment_results %>% 
  filter(annotation %in% c("heart_ind_eQTL", "artery_ind_eQTL"))

# Define a color palette for methods
color_palette <- c("MESuSiE" = "#0072B2", "Paintor" = "#F0E442", 
                   "Paintor_all" = "#CC79A7", "MLK" = "#009E73")

# Create Panel A plot: Enrichment (Signal / Background)
plotA <- ggplot(panelA, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel A: Enrichment Results (Signal / Background)",
       x = "MESuSiE_PIP_Either_all_mlk Threshold",
       y = "Enrichment (Signal / Background)") +
  scale_color_manual(values = color_palette) +
  theme_bw()

# Create Panel B plot: Average Scores
plotB <- ggplot(panelB, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel B: Average Scores",
       x = "MESuSiE_PIP_Either_all_mlk Threshold",
       y = "Average Score") +
  scale_color_manual(values = color_palette) +
  theme_bw()

# Create Panel C plot: Separate enrichment for "heart_ind_eQTL" and "artery_ind_eQTL"
plotC <- ggplot(panelC, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel C: Enrichment for heart_ind_eQTL and artery_ind_eQTL",
       x = "MESuSiE_PIP_Either_all_mlk Threshold",
       y = "Enrichment (Signal / Background)") +
  scale_color_manual(values = color_palette) +
  theme_bw()

# Combine the three panels vertically using patchwork
combined_plot <- plotA  / plotC

# Display the combined plot
combined_plot

# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "test_eqtl_bp_updated.png"),
       combined_plot, width = 225, height = 300, dpi = 500, units = 'mm')







##################################################

## Evaluate based on absolute top cutoff 50 to 100

##################################################


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Annotation names
ann_col_name <- c(
  "non_synonymous_common", "synonymous", "promotor", "utr_comb",
  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
  "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL", "artery_ind_eQTL", "CADD"
)

# Which annotations are binary vs. continuous
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

# Define list sizes to evaluate
list_sizes <- seq(50, 100, by = 10)

enrichment_results <- map_df(list_sizes, function(N) {
  
  # How many SNPs to take for each method
  top_N_signal <- N
  
  # Compute background frequencies (binary only)
  bg_bin <- res_all %>%
    summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))
  
  # Helper to do one method:
  do_method <- function(method_name, pip_col) {
    subset <- res_all %>%
      arrange(desc(.data[[pip_col]])) %>%
      top_n(n = top_N_signal, wt = .data[[pip_col]])
    
    bin_signal  <- subset %>% summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
    bin_enrich  <- bin_signal / bg_bin
    
    cadd_mean   <- subset %>% summarise(CADD = mean(CADD, na.rm = TRUE))
    gtex_mean   <- subset %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
    
    bind_rows(
      bin_enrich  %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment"),
      cadd_mean   %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment"),
      gtex_mean   %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment")
    ) %>%
      mutate(method    = method_name,
             list_size = top_N_signal)
  }
  
  # Apply to each method
  bind_rows(
    do_method("MESuSiE",              "MESuSiE_PIP_Either"),
    do_method("MLK",                  "MESuSiE_PIP_Either_all_mlk"),
    do_method("GLM",                  "MESuSiE_PIP_Either_all_glmnet"),
    do_method("Paintor",              "Paintor_PIP"),
    do_method("Paintor_all",          "Paintor_PIP_all")
  )
})

# Tidy up factor levels and labels
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation,
                             levels = ann_col_name),
         annotation = recode(annotation,
                             non_synonymous_common        = "Non-synonymous",
                             synonymous                   = "Synonymous",
                             promotor                     = "Promotor",
                             utr_comb                     = "UTR",
                             Conserved_LindbladToh_common = "Conserved",
                             H3K4me3_Trynka_common        = "H3K4me3",
                             GTEx_eQTL_MaxCPP_common      = "GTEx_eQTL_MaxCPP",
                             heart_ind_eQTL               = "Heart eQTL",
                             artery_ind_eQTL              = "Artery eQTL",
                             CADD                         = "CADD"),
         method = factor(method, levels = c("MESuSiE","Paintor","Paintor_all","GLM","MLK"))
  )

# Plot: enrichment vs. top‑N list size
enrichment_plot <- ggplot(enrichment_results, aes(x = list_size, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Enrichment of Functional Annotations Across Top‑N SNP Lists",
    x     = "Number of top‑ranked SNPs (N)",
    y     = "Enrichment (Signal / Background)"
  ) +
  scale_color_manual(values = c(
    MESuSiE     = "#0072B2",
    Paintor     = "#F0E442",
    Paintor_all = "#CC79A7",
    GLM         = "#E69F00",
    MLK         = "#009E73"
  )) +
  theme_bw()

# Save
ggsave("enrichment_by_list_size.png", enrichment_plot,
       width = 225, height = 150, units = "mm", dpi = 500)





##################################################

## Evaluate based on pip cutoff 0.5 to 0.9 with
## respective to each method, no same SNPs cutoff

##################################################





library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# 1. Annotation names
ann_col_name <- c(
  "non_synonymous_common", "synonymous", "promotor", "utr_comb",
  "Conserved_LindbladToh_common", "H3K4me3_Trynka_common", 
  "GTEx_eQTL_MaxCPP_common", "heart_ind_eQTL", "artery_ind_eQTL", "CADD"
)

# 2. Thresholds
thresholds <- seq(0.5, 0.9, by = 0.05)

# 3. Which are binary?
binary_ann <- setdiff(ann_col_name, c("CADD", "GTEx_eQTL_MaxCPP_common"))

# 4. Map methods → their PIP columns
method_pips <- c(
  MESuSiE      = "MESuSiE_PIP_Either",
  MLK          = "MESuSiE_PIP_Either_all_mlk",
  GLM          = "MESuSiE_PIP_Either_all_glmnet",
  Paintor      = "Paintor_PIP",
  Paintor_all  = "Paintor_PIP_all"
)

# 5. Compute enrichment by filtering on each method’s PIP > threshold
enrichment_results <- map_df(thresholds, function(thresh) {
  
  # for each method, get its subset, background, and compute ratios
  map_dfr(names(method_pips), function(meth) {
    pip_col <- method_pips[meth]
    
    # signal set: all SNPs with PIP > thresh for this method
    sig <- res_all %>% filter(.data[[pip_col]] > thresh)
    n_sig <- nrow(sig)
    print(n_sig)
    n_total <- nrow(res_all)
    n_bg <- n_total - n_sig
    
    # background frequencies for binary annotations
    bg_bin <- res_all %>%
      summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n_bg))
    
    # signal frequencies for binary annotations
    sig_bin <- sig %>%
      summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n_sig))
    
    # enrichment = signal / background
    bin_enrich <- sig_bin / bg_bin
    
    # continuous means in the signal
    cadd_mean <- sig %>% summarise(CADD = mean(CADD, na.rm = TRUE))
    gtex_mean <- sig %>% summarise(GTEx_eQTL_MaxCPP_common = mean(GTEx_eQTL_MaxCPP_common, na.rm = TRUE))
    
    # bind into long format
    bind_rows(
      bin_enrich   %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment"),
      cadd_mean    %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment"),
      gtex_mean    %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment")
    ) %>%
      mutate(
        method    = meth,
        threshold = thresh
      )
  })
})

# 6. (unchanged) Factor levels, recoding and plotting
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, levels = ann_col_name)) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common"         = "Non‑synonymous",
                             "synonymous"                    = "Synonymous",
                             "promotor"                      = "Promotor",
                             "utr_comb"                      = "UTR",
                             "Conserved_LindbladToh_common"  = "Conserved",
                             "H3K4me3_Trynka_common"         = "H3K4me3",
                             "GTEx_eQTL_MaxCPP_common"       = "GTEx_eQTL_MaxCPP",
                             "heart_ind_eQTL"                = "heart_ind_eQTL",
                             "artery_ind_eQTL"               = "artery_ind_eQTL",
                             "CADD"                          = "CADD"
  )) %>%
  mutate(method = factor(method, levels = c("MESuSiE","Paintor","Paintor_all","GLM","MLK")))

enrichment_plot <- ggplot(enrichment_results, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Enrichment Results Across PIP Cutoffs by Method",
    x     = "PIP Threshold",
    y     = "Enrichment (Signal / Background)"
  ) +
  scale_color_manual(values = c(
    MESuSiE      = "#0072B2",
    Paintor      = "#F0E442",
    Paintor_all  = "#CC79A7",
    GLM          = "#E69F00",
    MLK          = "#009E73"
  )) +
  theme_bw()

# save
ggsave(
  filename = paste0(plot_dir,"enrichment_results_by_method.png"),
  plot     = enrichment_plot,
  width    = 225, height = 150,
  units    = "mm", dpi = 500
)
