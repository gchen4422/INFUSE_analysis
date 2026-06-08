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
library(purrr)  # for map_df

all_sets <- list()
all_CS_enrichment <- list()
all_Signal_enrichment <- list()

for(num_causal in 1:3){
  for(h2_num in 1:2){
    
    wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/causal_num_", num_causal, "/")
    out_dir <- paste0(wrk_dir, "out/")
    load(paste0(out_dir, "res_simulation_h2_", h2_num, ".RData"))
    res_all$Trait <- "simulation"
    colnames(res_all)[c(18,23)] <- c("MESuSiE_res_glmnet_annot_all_cs", "MESuSiE_res_mlk_annot_all_cs")
    
    ############################################################################################################
    #   Stats to report
    ############################################################################################################
    
    # Count unique regions
    unique_regions <- res_all %>% summarise(n_unique_regions = n_distinct(Region))
    print(unique_regions)
    
    # Summarize the number of SNPs by region
    snp_summary <- res_all %>% 
      group_by(Region) %>% 
      summarise(NSNP = n()) %>%
      summarise(summary(NSNP))
    print(snp_summary)
    
    # Calculate position length by region and summarize
    pos_length_summary <- res_all %>%
      group_by(Region) %>%
      summarise(max_pos = max(as.numeric(POS)),
                min_pos = min(as.numeric(POS)),
                pos_length = as.numeric(max_pos) - as.numeric(min_pos)) %>%
      summarise(summary(pos_length / (1024 * 1024)))
    print(pos_length_summary)
    
    # 95% credible set size calculations (example code)
    region_cs <- res_all %>% 
      group_by(Region) %>% 
      summarise(across(c("MESuSiE_cs", "SuSiE_cs", "Paintor_cs", "Paintor_cs_all",
                         "MESuSiE_res_mlk_annot_all_cs", "MESuSiE_res_glmnet_annot_all_cs"), 
                       ~ sum(.x, na.rm = TRUE)))
    region_cs %>% summarise(sum(MESuSiE_cs != 0),
                            sum(SuSiE_cs != 0),
                            sum(Paintor_cs != 0),
                            sum(Paintor_cs_all != 0),
                            sum(MESuSiE_res_mlk_annot_all_cs != 0),
                            sum(MESuSiE_res_glmnet_annot_all_cs != 0))
    res_all %>%
      group_by(Region) %>%
      summarise(
        across(c("MESuSiE_cs", "SuSiE_cs", "Paintor_cs", "Paintor_cs_all",
                 "MESuSiE_res_mlk_annot_all_cs", "MESuSiE_res_glmnet_annot_all_cs"), 
               ~ sum(.x, na.rm = TRUE), .names = "sum_{col}")
      ) %>%
      summarise(
        MESuSiE_cs_median = median(sum_MESuSiE_cs[sum_MESuSiE_cs != 0], na.rm = TRUE),
        SuSiE_cs_median = median(sum_SuSiE_cs[sum_SuSiE_cs != 0], na.rm = TRUE),
        Paintor_cs_median = median(sum_Paintor_cs[sum_Paintor_cs != 0], na.rm = TRUE),
        Paintor_cs_all_median = median(sum_Paintor_cs_all[sum_Paintor_cs_all != 0], na.rm = TRUE),
        MESuSiE_res_mlk_annot_all_cs_median = median(sum_MESuSiE_res_mlk_annot_all_cs[sum_MESuSiE_res_mlk_annot_all_cs != 0], na.rm = TRUE),
        MESuSiE_res_glmnet_annot_all_cs_median = median(sum_MESuSiE_res_glmnet_annot_all_cs[sum_MESuSiE_res_glmnet_annot_all_cs != 0], na.rm = TRUE)
      )
    
    # Median |Z-scores| of SNPs within credible set
    MESuSiE_median_z <- res_all %>% filter(MESuSiE_cs == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    MESuSiE_mlk_median_z <- res_all %>% filter(MESuSiE_res_mlk_annot_all_cs == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    MESuSiE_glm_median_z <- res_all %>% filter(MESuSiE_res_glmnet_annot_all_cs == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    SuSiE_median_z <- res_all %>% filter(SuSiE_cs == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    Paintor_median_z <- res_all %>% filter(Paintor_cs == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    Paintor_all_median_z <- res_all %>% filter(Paintor_cs_all == 1) %>% 
      summarise(zmax = median(pmax(abs(zscore_1), abs(zscore_2)))) %>% pull(zmax)
    
    cat(round(c(MESuSiE_median_z, MESuSiE_mlk_median_z, MESuSiE_glm_median_z,
                SuSiE_median_z, Paintor_median_z, Paintor_all_median_z), 2))
    
    
    ann_col_name <- c("H3K27ac_Hnisz_common","DHS_Trynka_common",
                      "H3K4me3_Trynka_common","Conserved_LindbladToh_common",
                      "non_synonymous_common")
    
    # Define thresholds from 0.5 to 0.9 (adjust step as needed)
    thresholds <- seq(0.5, 0.9, by = 0.05)
    
    # Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
    binary_ann <- c("H3K27ac_Hnisz_common","DHS_Trynka_common",
                    "H3K4me3_Trynka_common","Conserved_LindbladToh_common",
                    "non_synonymous_common")
    
    enrichment_results <- map_df(thresholds, function(thresh) {
      
      # Define top_N_signal based on MESuSiE_PIP_Either_all_mlk exceeding the current threshold
      top_N_signal <- res_all %>% filter(MESuSiE_PIP_Either_all_mlk > thresh) %>% nrow()
      
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
      
      mes_long <- mes_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "MESuSiE", threshold = thresh)
      
      ## MLK method
      mlk_subset <- res_all %>% 
        arrange(desc(MESuSiE_PIP_Either_all_mlk)) %>% 
        top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)
      
      mlk_bin_signal <- mlk_subset %>% 
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
      mlk_bin_enrich <- mlk_bin_signal / bg_bin
      
      mlk_long <- mlk_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "MLK", threshold = thresh)
      
      ## GLM method
      glm_subset <- res_all %>% 
        arrange(desc(MESuSiE_PIP_Either_all_glmnet)) %>% 
        top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_glmnet)
      
      glm_bin_signal <- glm_subset %>% 
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
      glm_bin_enrich <- glm_bin_signal / bg_bin
      
      glm_long <- glm_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "GLM", threshold = thresh)
      
      ## susie 
      susie_subset <- res_all %>% 
        arrange(desc(SuSiE_PIP)) %>% 
        top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_glmnet)
      
      susie_bin_signal <- susie_subset %>% 
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
      susie_bin_enrich <- susie_bin_signal / bg_bin
      
      susie_long <- susie_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "SuSiE", threshold = thresh)
      
      
      ## Paintor method
      paintor_subset <- res_all %>% 
        arrange(desc(Paintor_PIP)) %>% 
        top_n(n = top_N_signal, wt = Paintor_PIP)
      
      paintor_bin_signal <- paintor_subset %>% 
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
      paintor_bin_enrich <- paintor_bin_signal / bg_bin
      
      paintor_long <- paintor_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "Paintor", threshold = thresh)
      
      ## Paintor_all method
      paintor_all_subset <- res_all %>% 
        arrange(desc(Paintor_PIP_all)) %>% 
        top_n(n = top_N_signal, wt = Paintor_PIP_all)
      
      paintor_all_bin_signal <- paintor_all_subset %>% 
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE) / n()))
      paintor_all_bin_enrich <- paintor_all_bin_signal / bg_bin
      
      paintor_all_long <- paintor_all_bin_enrich %>% 
        pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>% 
        mutate(method = "Paintor_all", threshold = thresh)
      
      # Combine all methods for this threshold
      bind_rows(mes_long, mlk_long, glm_long,susie_long, paintor_long, paintor_all_long)
    })
    
    # Reorder and recode the 'annotation' factor
    enrichment_results <- enrichment_results %>%
      mutate(annotation = factor(annotation, 
                                 levels = c("non_synonymous_common",
                                            "Conserved_LindbladToh_common", 
                                            "DHS_Trynka_common",
                                            "H3K4me3_Trynka_common",
                                            "H3K27ac_Hnisz_common"))) %>%
      mutate(annotation = recode(annotation,
                                 "non_synonymous_common" = "Non-synonymous",
                                 "Conserved_LindbladToh_common" = "Conserved",
                                 "H3K4me3_Trynka_common" = "H3K4me3",
                                 "DHS_Trynka_common" = "DHS",
                                 "H3K27ac_Hnisz_common" = "H3K27ac"))
    
    enrichment_results$method <- factor(enrichment_results$method, 
                                        levels = c("MESuSiE", "Paintor", "Paintor_all","SuSiE","GLM", "MLK"))
    
    # Create and print the per-iteration enrichment plot (optional)
    enrichment_plot <- ggplot(enrichment_results, aes(x = threshold, y = enrichment, color = method)) +
      geom_line() +
      facet_wrap(~ annotation, scales = "free_y") +
      labs(
        title = paste("Enrichment Results for causal =", num_causal, "h2 =", h2_num),
        x = "MESuSiE_PIP_Either_all_mlk Threshold",
        y = "Enrichment (Signal / Background)"
      ) +
      scale_color_manual(values = c(
        "MESuSiE"     = "#0072B2",   # Blue
        "Paintor"     = "#F0E442",   # Yellow
        "Paintor_all" = "#CC79A7",   # Purple
        "GLM"         = "#E69F00",   # Orange-Red
        "MLK"         = "#009E73"    # Green
      )) +
      theme_bw()
    print(enrichment_plot)
    
    # Append this iteration's enrichment results to the list for later averaging
    all_Signal_enrichment[[paste0("num_causal_", num_causal, "_h2_", h2_num)]] <- enrichment_results
  }
}

#########################
# Aggregated Plot: Average Enrichment Across All Simulation Runs
#########################

# Combine all enrichment results across simulations
combined_enrichment <- bind_rows(all_Signal_enrichment)

# Compute the average enrichment grouped by method, annotation, and threshold
avg_enrichment <- combined_enrichment %>%
  group_by(method, annotation, threshold) %>%
  summarise(avg_enrichment = mean(enrichment, na.rm = TRUE)) %>%
  ungroup()

# Reorder the annotation factor if needed
avg_enrichment$annotation <- factor(avg_enrichment$annotation, 
                                    levels = c("Non-synonymous", "Conserved", "DHS", "H3K4me3", "H3K27ac"))

# Create the aggregated enrichment plot
avg_enrichment_plot <- ggplot(avg_enrichment, aes(x = threshold, y = avg_enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Average Enrichment Across All num_causal and h2_num Combinations",
    x = "MESuSiE_PIP_Either_all_mlk Threshold",
    y = "Average Enrichment (Signal / Background)"
  ) +
  scale_color_manual(values = c(
    "MESuSiE"     = "#0072B2",
    "Paintor"     = "#F0E442",
    "Paintor_all" = "#CC79A7",
    "GLM"         = "#E69F00",
    "MLK"         = "#009E73"
  )) +
  theme_bw()



plot_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/plot/Enrichment/"
ggsave(paste0(plot_dir,"Summary_Enrichment_Comparison_Set_Top_Either.pdf"),avg_enrichment_plot,dpi=500,width=180,height =120,units = 'mm')

##########################

###Bar plot at PIP = 0.5

##########################

Signal_enrichment_Combined = avg_enrichment %>% filter(threshold == 0.5) %>% filter(method != "GLM")

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
    strip.background = element_blank(),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12, face="bold"),
    plot.title = element_text(size=7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

Signal_enrichment_Combined$annotation<-factor(	 Signal_enrichment_Combined$annotation,levels = c("H3K27ac","DHS","H3K4me3","Conserved","Non-synonymous"))
Signal_enrichment_Combined$method<-factor(Signal_enrichment_Combined$method,levels =c("MESuSiE","MLK", "SuSiE","Paintor","Paintor_all"))


bar_width <- 0.8
pd        <- position_dodge(width = bar_width)

p_signal <- ggplot(
  data = Signal_enrichment_Combined,
  aes(x = annotation, y = avg_enrichment, fill = method)
) +
  # use explicit dodge width
  geom_col(
    position = pd,
    width    = bar_width
  ) +
  scale_fill_manual(values = c(
    "MESuSiE"      = "#84c3b7",
    "MLK"          = "#71b8ed",
    "SuSiE"        = "#fbe79e",
    "Paintor"      = "#f2b56e",
    "Paintor_all"  = "#f57c6e"
  )) +
  # put labels at top‐center of each bar
  geom_text(
    aes(label = round(avg_enrichment, 2)),
    position = pd,
    vjust    = -0.3,    # slightly above bar
    hjust    = 0.5,     # center text
    size     = 5*5/10
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = NULL, y = "Fold Enrichment Top Signal") +
  ylim(0, round(max(Signal_enrichment_Combined$avg_enrichment)) + 1) +
  theme_bw()

p_signal

load(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/plot/summarize_low_dim_w2_baseline.Rdata"))

# 1) Completely remove any legend from p_signal


p_roc_either = p_roc_either+
  # this drops *any* fill, colour, size, etc. legends
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none")

p_size_box$data <- p_size_box$data %>%
  mutate(
    Method = fct_recode(Method,
      INFUSE         = "MESuSiE_fun",
      `Paintor+annot` = "Paintor_fun"
    ),
    Method = fct_relevel(Method,
       "SuSiE", "Paintor", "Paintor+annot","MESuSiE","INFUSE"
    )
  )

p_size_box = p_size_box +
  scale_fill_manual(values = c(
    "MESuSiE"       = "#84c3b7",
    "INFUSE"         = "#71b8ed",
    "SuSiE"         = "#fbe79e",
    "Paintor"       = "#f2b56e",
    "Paintor+annot" = "#f57c6e"
  ))+
  # this drops *any* fill, colour, size, etc. legends
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none")

p_power_bar$data <- p_power_bar$data %>%
  mutate(
    Method = fct_recode(Method,
                        INFUSE         = "MESuSiE_fun",
                        `Paintor+annot` = "Paintor_fun"
    ),
    Method = fct_relevel(Method,
                         "SuSiE", "Paintor", "Paintor+annot","MESuSiE","INFUSE"
    ))


p_power_bar$layers <- Filter(
  function(l) !inherits(l$geom, "GeomText"),
  p_power_bar$layers
)


p_power_bar = p_power_bar +
  scale_fill_manual(values = c(
    "MESuSiE"       = "#84c3b7",
    "INFUSE"         = "#71b8ed",
    "SuSiE"         = "#fbe79e",
    "Paintor"       = "#f2b56e",
    "Paintor+annot" = "#f57c6e"
  ))+
  # put labels at top‐center of each bar
  geom_text(
    aes(label = round(p_power_bar$data$Power_name, 2)),
    position = pd,
    vjust    = -0.3,    # slightly above bar
    hjust    = 0.5,     # center text
    size     = 5*5/10
  ) +
  
  # this drops *any* fill, colour, size, etc. legends
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none")




p_signal_nolegend <- p_signal +
  # this drops *any* fill, colour, size, etc. legends
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none")


p_FDR_Power$data <- p_FDR_Power$data %>%
  mutate(
    Method = fct_recode(Method,
                        INFUSE         = "MESuSiE_fun",
                        `Paintor+annot` = "Paintor_fun"
    ),
    Method = fct_relevel(Method,
                         "SuSiE", "Paintor", "Paintor+annot","MESuSiE","INFUSE"
    ))


p_FDR_Power$layers <- Filter(
  function(l) !inherits(l$geom, "GeomText"),
  p_FDR_Power$layers
)


pd <- position_dodge(width = 0.9)

p_FDR_Power = p_FDR_Power +
  scale_fill_manual(values = c(
    "MESuSiE"       = "#84c3b7",
    "INFUSE"         = "#71b8ed",
    "SuSiE"         = "#fbe79e",
    "Paintor"       = "#f2b56e",
    "Paintor+annot" = "#f57c6e"
  ))+
  # put labels at top‐center of each bar
  geom_text(
    aes(label = round(p_FDR_Power$data$Power, 2)),
    position = pd,
    vjust    = -0.3,    # slightly above bar
    hjust    = 0.5,     # center text
    size     = 5*5/10
  )

size_power<-p_size_box/p_power_bar
roc_fdr_power<-p_roc_either/p_FDR_Power+theme(legend.position = 'bottom')

p_total <- (
  size_power | roc_fdr_power
) /
  p_signal_nolegend +
  plot_annotation(tag_levels = 'a') +
  plot_layout(
    guides = "collect",
    # make top row 1 unit tall, bottom row only 0.3 units tall:
    heights = c(1, 0.3),
    # keep the two top panels equal width
    widths  = c(1, 1)
  ) & 
  theme(
    plot.tag        = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )+ custom_theme()
p_total
ggsave(paste0(plot_dir,"summarize_enrichment_low_dim_w2_baseline.pdf"),p_total,height=9.4, width=9, units = "in",dpi=500)



save(p_total,file = paste0(plot_dir,"summarize_low_dim_w2_baseline.Rdata"))



####################################################

## Mean square error plot of enrichment at PIP = 0.5

####################################################


load("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/plot/Enrichment//True_enrichment_w2_either.Rdata")


custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
    strip.background = element_blank(),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12, face="bold"),
    plot.title = element_text(size=12, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

Signal_enrichment_5 = combined_enrichment%>% filter(threshold == 0.5) %>% filter(method != "GLM")


# 1) define a lookup from your Cat → annotation names
cat_map <- c(
  H3K27ac_Hnisz_common        = "H3K27ac",
  DHS_Trynka_common           = "DHS",
  H3K4me3_Trynka_common       = "H3K4me3",
  Conserved_LindbladToh_common = "Conserved",
  non_synonymous_common       = "Non-synonymous"
)

# add an 'annotation' column to the true‐values df
true2 <- True_enrichment_long %>%
  mutate(annotation = cat_map[as.character(Cat)])

# 2) join to the estimated enrichments
merged <- Signal_enrichment_5 %>%
  inner_join(true2, by = c("annotation", "causal", "h2"))


mse_by_annotation <- merged %>%
  group_by(method, causal, h2, annotation) %>%
  summarise(
    mse = mean((enrichment - Prop)^2, na.rm = TRUE),
    .groups = "drop"
  )

mse_by_annotation_mean = mse_by_annotation%>%group_by(method,annotation)%>%summarise(mse= round(mean(mse),2))


bar_width <- 0.8
pd        <- position_dodge(width = bar_width)


mse_by_annotation <- mse_by_annotation %>%
  mutate(
    method = recode(
      method,
      MLK         = "INFUSE",
      Paintor_all = "Paintor+annot"
    )
  ) %>%
  # now force the order you like:
  mutate(
    method = factor(
      method,
      levels = c("SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE")
    )
  )

# same manual fill colors as p_signal
method_colors <- c(
  "MESuSiE"     = "#84c3b7",
  "INFUSE"         = "#71b8ed",
  "SuSiE"       = "#fbe79e",
  "Paintor"     = "#f2b56e",
  "Paintor+annot" = "#f57c6e"
)

p_mse_box <- ggplot(mse_by_annotation,
                    aes(x = method, y = mse, fill = method)) +
  # side‐by‐side boxplots
  geom_boxplot(position = position_dodge(width = 0.7),
               width = 0.7,
               outlier.size = 0.5) +
  # custom colors
  scale_fill_manual(values = method_colors) +
  # zero‐reference line
  geom_hline(yintercept = 0, linetype = "dashed") +
  # one facet per annotation, all in a single row
  facet_wrap(~ annotation,
             scales = "free_y",
             nrow   = 1,
             ncol   = 5) +
  # axis labels and title
  labs(
    x     = NULL,
    y     = "MSE of Enrichment",
    title = "MSE Distribution by Annotation\n(d=10)"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(hjust = 0.5),
    strip.text       = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_mse_box + custom_theme()
save(p_mse_box,file = "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/plot/MSE_low_dim_w2.Rdata")

