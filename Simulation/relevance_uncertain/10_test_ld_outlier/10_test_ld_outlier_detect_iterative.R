library(Rcpp)
library(MESuSiE)
library(susieR)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/modified_kriging_rss.R")
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(patchwork)


#args <- as.numeric(commandArgs(TRUE))

results <- expand.grid(
  LD_BLOCK   = 1:100,
  h2_num     = 1:2,
  num_causal = 1:3
) %>%
  mutate(
    flag_eu_initial    = integer(n()),
    flag_bb_initial    = integer(n()),
    flag_joint_initial = integer(n()),
    tp_combined = integer(n()),
    fp_combined = integer(n()),
    fn_combined = integer(n()),
    joint_tp    = integer(n()),
    joint_fp    = integer(n()),
    joint_fn    = integer(n()),
    # ← NEW: initialize FDR columns
    fdr_combined = numeric(n()),
    fdr_joint    = numeric(n())
  )

for (i in 1:nrow(results)) {
  
  
  LD_BLOCK = results$LD_BLOCK[i]
  h2_num = results$h2_num[i]
  num_causal = results$num_causal[i]
  
  
  cat(sprintf(
    "Processing LD_BLOCK=%d | h2_num=%d | num_causal=%d\n",
    LD_BLOCK, h2_num, num_causal
  ))
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_",num_causal,"/")
  data_dir<-paste0(wrk_dir,"summary_data/w2_startover/")
  
  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
  LD_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
  LD_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  
  zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
  summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
  summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )
  #summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 
  
  EU_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
  colnames(EU_cov)<-zfile$RSID
  BB_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
  colnames(BB_cov)<-zfile$RSID
  #R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)
  
  
  # ———— INITIAL DIAGNOSTICS ————
  
  EU_diagnostic <- kriging_rss(summary_stat_1$Z, EU_cov, n = median(summary_stat_1$N))
  BB_diagnostic <- kriging_rss(summary_stat_2$Z, BB_cov, n = median(summary_stat_2$N))
  joint_diagnostic <- kriging_rss_joint(summary_stat_1$Z, summary_stat_2$Z, EU_cov,BB_cov,median(summary_stat_1$N), median(summary_stat_2$N))
  
  
  results$flag_eu_initial[i]    <- length(BB_diagnostic$plot$plot_env$idx)
  results$flag_bb_initial[i]    <- length(EU_diagnostic$plot$plot_env$idx)
  idx_flag_joint = ifelse(joint_diagnostic$plot$plot_env$idx_flag_joint > nrow(summary_stat_1),joint_diagnostic$plot$plot_env$idx_flag_joint - nrow(summary_stat_1),joint_diagnostic$plot$plot_env$idx_flag_joint)
  results$flag_joint_initial[i] <- length(unique(c(joint_diagnostic$plot$plot_env$idx_flag1,joint_diagnostic$plot$plot_env$idx_flag2,idx_flag_joint)))
  
  
  
  # ——— MANUAL FLIPS (with tryCatch guards) ——————————————————
  # 1) Identify SNPs with |Z| in each ancestry
  idx1_all <- which(abs(summary_stat_1$Z) > 0)
  idx2_all <- which(abs(summary_stat_2$Z) > 0)
  
  # 2) Pick n_common SNPs to flip in both ancestries
  n_common <- 0
  common  <- intersect(idx1_all, idx2_all)
  if(length(common) < n_common) stop("Not enough common significant SNPs")
  idx_common <- sample(common, n_common)
  
  # Initialize flipped vectors
  summary_stat_1$Zf <- summary_stat_1$Z
  summary_stat_2$Zf <- summary_stat_2$Z
  
  # Flip the common ones in both
  summary_stat_1$Zf[idx_common] <- - summary_stat_1$Z[idx_common]
  summary_stat_2$Zf[idx_common] <- - summary_stat_2$Z[idx_common]
  
  
  # 3) Now pick ancestry‐specific flips *excluding* the common ones
  
  # For EUR
  n_eur_only <- 6
  eur_only_pool <- setdiff(idx1_all, idx_common)
  idx1_specific <- tryCatch(
    sample(eur_only_pool, n_eur_only),
    error = function(e) {
      warning("Could only find ", length(eur_only_pool),
              " EUR‐only SNPs >2; flipping all of them")
      eur_only_pool
    }
  )
  summary_stat_1$Zf[idx1_specific] <- - summary_stat_1$Z[idx1_specific]
  
  
  
  # For AFR
  n_afr_only <- 6
  afr_only_pool <- setdiff(idx2_all, idx_common)
  idx2_specific <- tryCatch(
    sample(afr_only_pool, n_afr_only),
    error = function(e) {
      warning("Could only find ", length(afr_only_pool),
              " AFR‐only SNPs >2; flipping all of them")
      afr_only_pool
    }
  )
  summary_stat_2$Zf[idx2_specific] <- - summary_stat_2$Z[idx2_specific]
  
  
  idx_list = list(
    common      = idx_common,
    eur_specific = idx1_specific,
    afr_specific = idx2_specific
  )
  
  eu_df    <- kriging_rss(summary_stat_1$Zf, EU_cov, n = median(summary_stat_1$N))
  bb_df    <- kriging_rss(summary_stat_2$Zf, BB_cov, n = median(summary_stat_2$N))
  joint_df <- kriging_rss_joint(
    summary_stat_1$Zf, summary_stat_2$Zf,
    EU_cov, BB_cov,
    median(summary_stat_1$N), median(summary_stat_2$N)
  )
  
  eu_flagged    <- sort(eu_df$plot$plot_env$idx)
  bb_flagged    <- sort(bb_df$plot$plot_env$idx)
  idx_flag_joint = ifelse(joint_df$plot$plot_env$idx_flag_joint > nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint - nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint)
  joint_flagged <- sort(unique(c(
    joint_df$plot$plot_env$idx_flag1,
    joint_df$plot$plot_env$idx_flag2,
    idx_flag_joint)))
  
  
  
  
  # start with all SNPs
  keep_idx <- seq_along(summary_stat_1$Zf)
  iter     <- 0
  
  repeat {
    iter <- iter + 1
    
    # subset to the current set
    Zf1_sub    <- summary_stat_1$Zf[keep_idx]
    Zf2_sub    <- summary_stat_2$Zf[keep_idx]
    N1_sub     <- summary_stat_1$N [keep_idx]
    N2_sub     <- summary_stat_2$N [keep_idx]
    EU_cov_sub <- EU_cov   [keep_idx, keep_idx]
    BB_cov_sub <- BB_cov   [keep_idx, keep_idx]
    
    eu_df    <- kriging_rss(Zf1_sub, EU_cov_sub, n = median(N1_sub))
    bb_df    <- kriging_rss(Zf2_sub, BB_cov_sub, n = median(N2_sub))
    
    eu_flagged    <- sort(eu_df$plot$plot_env$idx)
    bb_flagged    <- sort(bb_df$plot$plot_env$idx)
    
    # print the iteration number and flagged indices
    cat("Iteration", iter, "- dropping indices eu:", 
        if (length(eu_flagged)) eu_flagged else "none", "\n")
    
    cat("Iteration", iter, "- dropping indices bb:", 
        if (length(bb_flagged)) bb_flagged else "none", "\n")
    
    
    # run your joint test
    joint_df <- kriging_rss_joint(
      Zf1_sub, Zf2_sub,
      EU_cov_sub, BB_cov_sub,
      median(N1_sub), median(N2_sub)
    )
    
    idx_flag_joint = ifelse(joint_df$plot$plot_env$idx_flag_joint > nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint - nrow(summary_stat_1),joint_df$plot$plot_env$idx_flag_joint)
    
    
    joint_flagged <- sort(unique(c(
      joint_df$plot$plot_env$idx_flag1,
      joint_df$plot$plot_env$idx_flag2,
      idx_flag_joint)))
    
    # print the iteration number and flagged indices
    cat("Iteration", iter, "- dropping indices:", 
        if (length(joint_flagged)) joint_flagged else "none", "\n")
    
    # stop if nothing new to drop
    if (length(joint_flagged) == 0) break
    
    # otherwise remove them and loop again
    keep_idx <- keep_idx[-joint_flagged]
  }
  
  
  # ground truth sets
  
  # 1) EUR “truth” = all SNPs flipped in EUR (common + EUR‐only)
  truth_eu <- sort(c(idx_list$common, idx_list$eur_specific))
  
  # 2) AFR “truth” = all SNPs flipped in AFR (common + AFR‐only)
  truth_bb <- sort(c(idx_list$common, idx_list$afr_specific))
  
  # 3) Joint “truth” = any SNP flipped in either ancestry
  truth_joint <- sort(union(truth_eu, truth_bb))
  
  
  
  
  
  # right after you define truth_* and flagged_*:
  n_snps_eur   <- nrow(summary_stat_1)      # total EUR SNPs
  n_snps_afr   <- nrow(summary_stat_2)      # total AFR SNPs
  n_true_eu    <- length(truth_eu)          
  n_true_bb    <- length(truth_bb)          
  n_true_joint <- length(truth_joint)       
  
  
  
  tp_c_raw <- length(
    unique(c(
      intersect(eu_flagged, truth_eu),
      intersect(bb_flagged, truth_bb)
    ))
  )
  
  fp_c_raw <- length(
    unique(c(
      setdiff(eu_flagged, truth_eu),
      setdiff(bb_flagged, truth_bb)
    ))
  )
  tp_j_raw <- length(intersect(joint_flagged, truth_joint))
  fp_j_raw <- length(setdiff(joint_flagged, truth_joint))
  
  
  
  # now store rates instead of counts:
  results$tp_combined[i] <- (
    length(unique(c(intersect(eu_flagged, truth_eu), intersect(bb_flagged, truth_bb))))
  ) / (n_true_joint)
  
  results$fn_combined[i] <- (
    length(unique(c(setdiff(truth_eu, eu_flagged), setdiff(truth_bb, bb_flagged))))
  ) / (n_true_joint)
  
  results$fp_combined[i] <- (
    length(unique(c(setdiff(eu_flagged, truth_eu), setdiff(bb_flagged, truth_bb))))
  ) / ((n_snps_eur + n_snps_afr - n_true_joint))
  
  results$fdr_combined[i] <- if ((tp_c_raw + fp_c_raw) > 0) {
    fp_c_raw / (tp_c_raw + fp_c_raw)
  } else {
    NA_real_
  }
  
  
  
  results$joint_tp[i] <- length(intersect(joint_flagged, truth_joint)) / n_true_joint
  results$joint_fn[i] <- length(setdiff(truth_joint,      joint_flagged)) / n_true_joint
  results$joint_fp[i] <- length(setdiff(joint_flagged,   truth_joint)) / (n_snps_eur + n_snps_afr - n_true_joint)
  results$fdr_joint[i] <- if ((tp_j_raw + fp_j_raw) > 0) {
    fp_j_raw / (tp_j_raw + fp_j_raw)
  } else {
    NA_real_
  }
  
}




# summary_metrics <- results %>%
#   summarise(
#     combined_sensitivity       = mean(tp_combined),
#     combined_false_negative    = mean(fn_combined),
#     combined_false_positive    = mean(fp_combined),
#     joint_sensitivity          = mean(joint_tp),
#     joint_false_negative       = mean(joint_fn),
#     joint_false_positive       = mean(joint_fp)
#   )


save(results, file = "10_test_ld_results.Rdata")

#estimate_s_rss(summary_stat_1$Z_flipped, EU_cov, n = median(summary_stat_1$N))
#estimate_s_rss(summary_stat_2$Z_flipped, BB_cov, n = median(summary_stat_2$N))
#estimate_s_rss_joint(summary_stat_1$Z_flipped, summary_stat_2$Z_flipped, EU_cov,BB_cov,median(summary_stat_1$N), median(summary_stat_2$N))
#estimate_s_rss_two(summary_stat_1$Z_flipped, summary_stat_2$Z_flipped, EU_cov,BB_cov,median(summary_stat_1$N), median(summary_stat_2$N))
#estimate_s_rss_two_penalty(summary_stat_1$Z_flipped, summary_stat_2$Z_flipped, EU_cov,BB_cov,median(summary_stat_1$N), median(summary_stat_2$N),rho = 1,alpha = 500)


#test <- kriging_rss_two(summary_stat_1$Z_flipped, summary_stat_2$Z_flipped, EU_cov,BB_cov,median(summary_stat_1$N), median(summary_stat_2$N),rho = 1,alpha = 5)
#test$plot$plot_env$idx1
#test$plot$plot_env$idx2



load("10_test_ld_results.Rdata")


results <- results %>%
  mutate(
    n_snp = integer(n()),
    type_i_error_combined    = numeric(n()),
    type_i_error_joint    = numeric(n())
  )






for (i in 1:nrow(results)) {
  
  
  LD_BLOCK = results$LD_BLOCK[i]
  h2_num = results$h2_num[i]
  num_causal = results$num_causal[i]
  
  
  cat(sprintf(
    "Processing LD_BLOCK=%d | h2_num=%d | num_causal=%d\n",
    LD_BLOCK, h2_num, num_causal
  ))
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_",num_causal,"/")
  data_dir<-paste0(wrk_dir,"summary_data/w2_startover/")
  
  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
  LD_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
  LD_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  
  zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
  results$n_snp[i] = nrow(zfile)
  
}


results = results %>% mutate(type_i_error_combined = (flag_eu_initial + flag_bb_initial)/n_snp,
                             type_i_error_joint = (flag_joint_initial)/n_snp,)


summary_metrics <- results %>%
  group_by(h2_num, num_causal) %>%
  summarise(
    combined_sensitivity    = mean(tp_combined),
    combined_false_negative = mean(fn_combined),
    combined_false_positive = mean(fp_combined),
    combined_fdr = mean(fdr_combined,na.rm = T),
    combined_type_i_error = mean(type_i_error_combined),
    joint_sensitivity       = mean(joint_tp),
    joint_false_negative    = mean(joint_fn),
    joint_false_positive    = mean(joint_fp),
    joint_fdr = mean(fdr_joint),
    joint_type_i_error = mean(type_i_error_joint),
    .groups = "drop"
  )


sens_long <- summary_metrics %>%
  dplyr::select(h2_num, num_causal, combined_sensitivity, joint_sensitivity) %>%
  pivot_longer(
    cols      = c(combined_sensitivity, joint_sensitivity),
    names_to  = "metric",
    values_to = "value"
  )


# False discovery rate (FDR)
fdr_long <- summary_metrics %>%
  dplyr::select(h2_num, num_causal, combined_fdr, joint_fdr) %>%
  pivot_longer(
    cols      = c(combined_fdr, joint_fdr),
    names_to  = "metric",
    values_to = "value"
  )

# False positives
fp_long <- summary_metrics %>%
  dplyr::select(h2_num, num_causal, combined_false_positive, joint_false_positive) %>%
  pivot_longer(
    cols      = c(combined_false_positive, joint_false_positive),
    names_to  = "metric",
    values_to = "value"
  )




# False negatives
fn_long <- summary_metrics %>%
  dplyr::select(h2_num, num_causal, combined_false_negative, joint_false_negative) %>%
  pivot_longer(
    cols      = c(combined_false_negative, joint_false_negative),
    names_to  = "metric",
    values_to = "value"
  )



# Type I error
type_i_error_long <- summary_metrics %>%
  dplyr::select(h2_num, num_causal, combined_type_i_error, joint_type_i_error) %>%
  pivot_longer(
    cols      = c(combined_type_i_error, joint_type_i_error),
    names_to  = "metric",
    values_to = "value"
  )

# 2) Build each plot

p_sens <- ggplot(sens_long, aes(
  x    = factor(num_causal),
  y    = value,
  fill = metric
)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  facet_wrap(~ h2_num, nrow = 1, labeller = label_both) +
  scale_fill_brewer(
    palette = "Set1",
    labels  = c("Combined Sensitivity", "Joint Sensitivity")
  ) +
  labs(
    x     = "Number of causal variants",
    y     = "Sensitivity",
    title = "Sensitivity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )


p_sens <- p_sens +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1)
  ) +
  labs(y = "Sensitivity")  # if you like percent labels


p_fdr <- ggplot(fdr_long, aes(
  x    = factor(num_causal),
  y    = value,
  fill = metric
)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  # dashed line at FDR = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
  facet_wrap(~ h2_num, nrow = 1, labeller = label_both) +
  scale_fill_brewer(
    palette = "Set3",
    labels  = c("Combined FDR", "Joint FDR")
  ) +
  labs(
    x     = "Number of causal variants",
    y     = "False discovery rate",
    title = "FDR"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  ) +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_y_continuous(
    breaks = seq(0, 0.1, by = 0.01),
    labels = seq(0, 0.1, by = 0.01)
  )




p_fp <- ggplot(fp_long, aes(
  x    = factor(num_causal),
  y    = value,
  fill = metric
)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  facet_wrap(~ h2_num, nrow = 1, labeller = label_both) +
  scale_fill_brewer(
    palette = "Set2",
    labels  = c("Combined FP", "Joint FP")
  ) +
  labs(
    x     = "Number of causal variants",
    y     = "False positives",
    title = "False Positives"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  ) +
  coord_cartesian(ylim = c(0, 0.0001)) +
  scale_y_continuous(
    breaks = seq(0, 0.0001, by = 0.00001),
    labels = seq(0, 0.0001, by = 0.00001)
  )




p_typei <- ggplot(type_i_error_long, aes(
  x    = factor(num_causal),
  y    = value,
  fill = metric
)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  facet_wrap(~ h2_num, nrow = 1, labeller = label_both) +
  scale_fill_brewer(
    palette = "Set2",
    labels  = c("Combined Type I Error", "Joint Type I Error")
  ) +
  labs(
    x     = "Number of causal variants",
    y     = "Type I error",
    title = "Type I error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  ) +
  coord_cartesian(ylim = c(0, 0.0001)) +
  scale_y_continuous(
    breaks = seq(0, 0.0001, by = 0.00001),
    labels = seq(0, 0.0001, by = 0.00001)
  )





p_fn <- ggplot(fn_long, aes(
  x    = factor(num_causal),
  y    = value,
  fill = metric
)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  facet_wrap(~ h2_num, nrow = 1, labeller = label_both) +
  scale_fill_brewer(
    palette = "Set3",
    labels  = c("Combined FN", "Joint FN")
  ) +
  labs(
    x     = "Number of causal variants",
    y     = "False negatives",
    title = "False Negatives"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )+
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1)
  )

# 3) Combine with patchwork

combined_fig <- (p_typei | p_sens | p_fdr) +
  plot_annotation(
    title = "Performance Metrics by Scenario",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Print it
print(combined_fig)
plot_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/plot/"

ggsave(
  filename = file.path(plot_dir, "performance_metrics_by_scenario_fdr.pdf"),
  plot     = combined_fig,      # explicitly pass your patchwork object
  width    = 16,                # inches
  height   = 5,                 # inches
  units    = "in", 
  device   = cairo_pdf
)

# 2) Save as high‐res TIFF
ggsave(
  filename    = file.path(plot_dir, "performance_metrics_by_scenario_fdr.tiff"),
  plot        = combined_fig,
  width       = 16,             # inches
  height      = 5,              # inches
  units       = "in",
  dpi         = 600,            # high resolution
  device      = "tiff",
  compression = "lzw"           # lossless
)


