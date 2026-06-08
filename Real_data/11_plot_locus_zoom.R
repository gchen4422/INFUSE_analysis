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
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/"
Gene_List<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/Gencode_genelist/Gencode_GRCh38_Genes_UniqueList2025.txt",header=T)

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



##############################################################
#
#
#       Real Data Example Plotter
#
#
################################################################

gwas_plot_fun <- function(data_plot, xlab_name, ylab_name, yintercept) {
  
  p_manhattan = ggplot() + geom_point(data = data_plot%>%filter(Lead_SNP==0), aes(x = POS, y = PIP, color = r2), size = 1)
  p_manhattan = p_manhattan + geom_point(data = data_plot%>%filter(Lead_SNP==1), aes(x = POS, y = PIP), size = 1.5, color = "red") +
    geom_text_repel(data = data_plot%>%filter(Lead_SNP==1), mapping = aes(x = POS, y = PIP, label = SNP), vjust = 1.2, size = 7/14*5, show.legend =FALSE) 
  p_manhattan = p_manhattan +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = 'grey50',
      name = expression(R^2)
    )
  p_manhattan = p_manhattan +
    geom_hline(
      yintercept = yintercept,
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) 
  p_manhattan = p_manhattan +
    geom_vline(
      xintercept = data_plot%>%filter(lead_SNP==1)%>%pull(POS),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) 
  p_manhattan = p_manhattan + xlim(min(data_plot$POS),max(data_plot$POS))
  p_manhattan = p_manhattan + expand_limits(x = round(max(data_plot$POS)/1e6)*1e6)
  if(max(data_plot$POS>1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e6, " MB"))
  }
  if(max(data_plot$POS<1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e3, " KB"))
  }
  p_manhattan = p_manhattan + xlab(xlab_name) +ylab(ylab_name)
  p_manhattan = p_manhattan + guides(fill = guide_legend(title = as.expression(bquote(R^2))))
  p_manhattan = p_manhattan + theme_bw()+custom_theme()
  return(p_manhattan)
}

###Function used for PIP plot	
finemap_plot_fun<-function(data_plot,xlab_name,ylab_name,yintercept){
  pos_rng <- range(data_plot$POS, na.rm = TRUE)
  
  p_manhattan = ggplot() + geom_point(data = data_plot, aes(x = POS, y = PIP, color = r2,shape = cat))+scale_shape_manual(name="Category",drop=FALSE,values=c(20,24,25,23,22))
  p_manhattan = p_manhattan + geom_text_repel(data =data_plot%>%filter(Lead_SNP==1), mapping=aes(x=POS, y=PIP, label=SNP),vjust=1.2, size= 7/14*5,show.legend = FALSE)
  p_manhattan = p_manhattan + theme_bw()+scale_color_stepsn(
    colors = c("navy", "lightskyblue", "green", "orange", "red"),
    breaks = seq(0.2, 0.8, by = 0.2),
    limits = c(0, 1),
    show.limits = TRUE,
    na.value = 'grey50',
    name = expression(R^2)
  )
  p_manhattan = p_manhattan + geom_hline(
    yintercept =yintercept,
    linetype = "dashed",
    color = "grey50",
    size = 0.5
  ) + geom_vline(
    xintercept = data_plot%>%filter(lead_SNP==1)%>%pull(POS),
    linetype = "dashed",
    color = "grey50",
    size = 0.5
  ) 
  p_manhattan = p_manhattan + xlim(min(data_plot$POS),max(data_plot$POS))
  p_manhattan = p_manhattan + expand_limits(x = round(max(data_plot$POS)/1e6)*1e6)
  if(max(data_plot$POS>1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e6, " MB"),expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng)
  }
  if(max(data_plot$POS<1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e3, " KB"),expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng)
  }
  p_manhattan= p_manhattan+xlab(xlab_name)+ylab(ylab_name)
  p_manhattan= p_manhattan+guides(fill=guide_legend(title=as.expression(bquote(R^2))))
  p_manhattan = p_manhattan + theme_bw()+custom_theme()
  return(p_manhattan)
}			
# Function used for gene plot 
gene_range_plot_fun<-function(gene_list_data,plot.range){
  p<-ggplot(data = gene_list_data) +
    geom_linerange(aes(x = Gene, ymin = Start, ymax = End))+ylim(plot.range)+ expand_limits(y = round(max(plot.range[2])/1e6)*1e6)+scale_y_continuous(limits = plot.range,labels = function(y) paste0(y / 1e6, " MB"),expand = expansion(mult = c(0.05, 0.05)),)+coord_flip()+
    geom_text(aes(x = Gene, y = Start, label = Gene), hjust = "right", size = 5/14*5) + ylab(paste0("chr",unique(gsub("chr","",gene_list_data$Chrom))))+ xlab("Gene") + 
    theme_bw() +  theme(
      axis.text.x = element_text(size = 7.5),
      axis.text.y = element_blank(),  
      axis.ticks.y =  element_blank(), 
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
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black")
    )
  return(p)
}

create_annotation_comparison <- function(candidate_region, 
                                       INFUSE_pip_col = "MESuSiE_PIP_Shared_all_mlk",
                                       mesusie_pip_col = "MESuSiE_PIP_Shared",
                                       custom_annotations = NULL) {
  
  # Default key annotations if not provided
  if (is.null(custom_annotations)) {
    key_annotations <- c(
      # Conserved elements
      "Conserved_LindbladToh_common",
      
      # Chromatin accessibility & marks
      "DHS_Trynka_common", 
      "H3K27ac_Hnisz_common", 
      "H3K4me3_Trynka_common",
      "Enhancer_Andersson_common",
      
      # Transcription factor binding
      "TFBS_ENCODE_common",
      
      # Regulatory regions
      "promotor",
      
      # Coding impact (for comparison)
      "Coding_UCSC_common",
      "Intron_UCSC_common",
      "missense",
      "synonymous"
    )
  } else {
    key_annotations <- custom_annotations
  }
  
  # Get top SNPs based on PIP columns
  top_INFUSE_snp <- candidate_region %>% 
    slice_max(.data[[INFUSE_pip_col]], n = 1) %>%
    pull(SNP)
  
  top_mesusie_snp <- candidate_region %>% 
    slice_max(.data[[mesusie_pip_col]], n = 1) %>%
    pull(SNP)
  
  # Check if annotations exist in the data
  available_annotations <- intersect(key_annotations, colnames(candidate_region))
  if (length(available_annotations) == 0) {
    stop("None of the specified annotations are found in the candidate_region data.")
  }
  
  if (length(available_annotations) < length(key_annotations)) {
    missing_annotations <- setdiff(key_annotations, available_annotations)
    warning(paste("Missing annotations:", paste(missing_annotations, collapse = ", ")))
  }
  
  # Create raw annotation comparison
  annotation_comparison_raw <- candidate_region %>%
    filter(SNP %in% c(top_INFUSE_snp, top_mesusie_snp)) %>%
    select(SNP, all_of(available_annotations))
  
  # Create cleaned annotation comparison
  annotation_comparison <- annotation_comparison_raw %>%
    # Convert to long format first
    pivot_longer(cols = -SNP, names_to = "Annotation", values_to = "Present") %>%
    # Clean annotation names
    mutate(
      Annotation = str_replace_all(Annotation, "_", " "),
      Annotation = str_replace_all(Annotation, "common", ""),
      Annotation = str_trim(Annotation),
      Present = as.numeric(Present),  # Convert to 0/1 for heatmap
      SNP_Label = case_when(
        SNP == top_INFUSE_snp ~ paste0(SNP, " (INFUSE Top)"),
        SNP == top_mesusie_snp ~ paste0(SNP, " (MESuSiE Top)"),
        TRUE ~ SNP
      )
    )
  
  # Return both raw and processed data along with top SNPs
  return(list(
    annotation_comparison = annotation_comparison,
    annotation_comparison_raw = annotation_comparison_raw,
    top_INFUSE_snp = top_INFUSE_snp,
    top_mesusie_snp = top_mesusie_snp,
    available_annotations = available_annotations,
    INFUSE_pip = candidate_region %>% filter(SNP == top_INFUSE_snp) %>% pull(.data[[INFUSE_pip_col]]),
    mesusie_pip = candidate_region %>% filter(SNP == top_mesusie_snp) %>% pull(.data[[mesusie_pip_col]])
  ))
}


####################################
#
# For shared causal variants for BMI
#
####################################

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData"))

# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_shared_col <- "MESuSiE_PIP_Shared_all_mlk"
mesusie_col      <- "MESuSiE_PIP_Shared"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_shared <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_shared_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_shared_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_shared %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_shared %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_shared_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_shared)
head(novel_regions)
head(novel_summary)






####################################
#
# For WB causal variants for BMI
#
####################################

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData"))


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_WB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_WB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_wb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_wb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_wb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_wb)
head(novel_regions)
head(novel_summary)



####################################
#
# For BB causal variants for BMI
#
####################################

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData"))


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_BB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_BB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_bb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_bb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_bb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_bb)
head(novel_regions)
head(novel_summary)


####################################
#
# For shared causal variants for DBP
#
####################################


load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "DBP")

# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_shared_col <- "MESuSiE_PIP_Shared_all_mlk"
mesusie_col      <- "MESuSiE_PIP_Shared"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_shared <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_shared_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_shared_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_shared %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_shared %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_shared_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_shared)
head(novel_regions)
head(novel_summary)




####################################
#
# For WB causal variants for DBP
#
####################################

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "DBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_WB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_WB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_wb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_wb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_wb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_wb)
head(novel_regions)
head(novel_summary)



####################################
#
# For BB causal variants for DBP
#
####################################

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "DBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_BB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_BB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_bb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_bb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_bb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_bb)
head(novel_regions)
head(novel_summary)



####################################
#
# For shared causal variants for SBP
#
####################################


load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "SBP")

# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_shared_col <- "MESuSiE_PIP_Shared_all_mlk"
mesusie_col      <- "MESuSiE_PIP_Shared"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_shared <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_shared_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_shared_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_shared %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_shared %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_shared_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_shared)
head(novel_regions)
head(novel_summary)


####################################
#
# For WB causal variants for SBP
#
####################################

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "SBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_WB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_WB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_wb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_wb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_wb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_wb)
head(novel_regions)
head(novel_summary)



####################################
#
# For BB causal variants for SBP
#
####################################

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "SBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
INFUSE_col <- "MESuSiE_PIP_BB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_BB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_bb <- res_all %>%
  filter(
    # called by INFUSE as shared
    .data[[INFUSE_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(INFUSE_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_bb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_bb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[INFUSE_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_bb)
head(novel_regions)
head(novel_summary)



###################################################################
#
#
#          BMI  Upstream of NGER and RPL31P12 eQTL shared signal
#
#
##################################################################


LD_BLOCK = 5
trait = "BMI"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs34361149 is the eQTL of RPL31P12 gene, highlighted in the paper
lead_SNP = "rs34361149"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                    MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                    TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs34361149")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs34361149")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_BB)

###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "UKBB", "-log10(P-value)", -log10(5e-8))
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "GLGC", "-log10(P-value)", -log10(5e-8))

###Finemap Plot
EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_EUR_SuSiE/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_AFR_SuSiE/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",region ,"_",trait_name,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)





###################################################################
#
#
#            BMI Misssense and sQTL of RPGRIP1L shared signal
#
#
##################################################################


LD_BLOCK = 168
trait = "BMI"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs34361149 is the eQTL of RPL31P12 gene, highlighted in the paper
lead_SNP = "rs61743997"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                    TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs61743997")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs61743997")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_BB)

###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "UKBB", "-log10(P-value)", -log10(5e-8))
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "GLGC", "-log10(P-value)", -log10(5e-8))

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",region ,"_",trait_name,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)












###################################################################
#
#
#             BMI Intron and eQTL of TFAP2B shared signal!!! https://pubmed.ncbi.nlm.nih.gov/22952648/
#             https://onlinelibrary.wiley.com/doi/10.1038/oby.2009.83
#
##################################################################

load(paste0(res_dir,"res_bmi_pca_bic_modified_llk.RData"))


drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)

LD_BLOCK = 81
trait = "BMI"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs34361149 is the eQTL of RPL31P12 gene, highlighted in the paper
lead_SNP = "rs3798519"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                          TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs3798519")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs3798519")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_WB)

annotation_results <- create_annotation_comparison(candidate_region)



###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5) + drop_guides + theme(legend.position = "none")
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)+theme(legend.position = "bottom")

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000& Gene !="EFHC1")
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together

combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",region ,"_",trait_name,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)













# Alternative version with even cleaner styling
p_heatmap_clean <- annotation_comparison %>%
  mutate(
    # Clean annotation names further
    Annotation = case_when(
      Annotation == "Conserved LindbladToh" ~ "Conserved",
      Annotation == "DHS Trynka" ~ "DHS",
      Annotation == "H3K27ac Hnisz" ~ "H3K27ac",
      Annotation == "H3K4me3 Trynka" ~ "H3K4me3",
      Annotation == "TFBS ENCODE" ~ "TFBS",
      TRUE ~ Annotation
    )
  ) %>%
  ggplot(aes(x = reorder(Annotation, -Present), y = SNP_Label, fill = factor(Present))) +
  geom_tile(color = "white", size = 1.5, width = 0.85, height = 0.85) +
  scale_fill_manual(
    values = c("0" = "grey95", "1" = "steelblue"),  # Slight grey for better contrast
    name = "",
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.position = "none",  # Remove legend for cleaner look
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_fixed(ratio = 0.8, expand = FALSE)



# Or use the cleaner version
p_heatmap_clean

###################################################################
#
#
#             BMI Intron and eQTL of BDNF WB signal https://pmc.ncbi.nlm.nih.gov/articles/PMC6121551/#abstract1
#             https://www.nature.com/articles/s41588-023-01392-0#Sec2
#       https://pubmed.ncbi.nlm.nih.gov/18160634/, https://pubmed.ncbi.nlm.nih.gov/18753648/
#
##################################################################



drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)

LD_BLOCK = 126
trait = "BMI"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs34361149 is the eQTL of RPL31P12 gene, highlighted in the paper
lead_SNP = "rs11030119"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                          TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs11030119")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs11030119")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_WB)


annotation_results <- create_annotation_comparison(candidate_region,"MESuSiE_PIP_Either_all_mlk")


###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5) + drop_guides + theme(legend.position = "none")
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)+theme(legend.position = "bottom")

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000& Gene !="EFHC1")
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together

combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)



###################################################################
#
#
#             BMI Intron and eQTL of INPP4B BB signal  https://www.nature.com/articles/s42003-021-01940-6
#
#
##################################################################


LD_BLOCK = 60
trait = "BMI"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs34361149 is the eQTL of RPL31P12 gene, highlighted in the paper
lead_SNP = "rs13140802"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                          TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs13140802")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs13140802")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_BB)

###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "UKBB", "-log10(P-value)", -log10(5e-8))
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "GLGC", "-log10(P-value)", -log10(5e-8))

###Finemap Plot
EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",region ,"_",trait_name,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)





###################################################################
#
#
#             DBP Missense and eQTL of SH2B3 shared signal, https://www.ahajournals.org/doi/10.1161/HYPERTENSIONAHA.118.12094
#             https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.121.320625 , https://www.snpedia.com/index.php/Rs3184504
#
##################################################################

drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)


LD_BLOCK = 28
trait = "DBP"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))


#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs3184504 is the eQTL of SH2B3 gene, highlighted in the paper
lead_SNP = "rs3184504"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                          TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs3184504")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs3184504")%>%dplyr::select(MESuSiE_PIP_Shared_all_mlk,MESuSiE_PIP_WB_all_mlk,MESuSiE_PIP_Shared_all_mlk,MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_BB)


annotation_results <- create_annotation_comparison(candidate_region)


###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)


annotation_results <- create_annotation_comparison(candidate_region,"MESuSiE_PIP_Either_all_mlk")


annotation_colors <- c(
  "Coding_UCSC_common"   = "#0072B2", # A nice blue
  "H3K27ac_Hnisz_common" = "#FFB6C1", # A light pink
  "Promoter_UCSC_common" = "#009E73", # A bluish green
  "default"              = "gray30"  # A dark gray for other annotations
)

plot_annotation_highlight <- function(data, 
                                      annotation_col, 
                                      highlight_snp_id,
                                      snp_col = "SNP",
                                      pos_col = "POS",
                                      chr_col = "CHR",
                                      plot_range = NULL,
                                      custom_theme = NULL) {
  
  # --- 1. PREPARE VARIABLES ---
  
  # Select the color for the '1' values from the palette
  track_color <- annotation_colors[[annotation_col]]
  if (is.null(track_color)) {
    track_color <- annotation_colors[["default"]]
  }
  
  plot_data <- data[data[[annotation_col]] == 1, ]


  # Convert annotation column to a factor for discrete coloring
  plot_data[[annotation_col]] <- as.factor(plot_data[[annotation_col]])

  chromosome_label <- paste0("chr", unique(plot_data[[chr_col]])[1])

  if (is.null(plot_range)) {
    pos_rng <- range(data[[pos_col]], na.rm = TRUE)
  } else {
    pos_rng <- plot_range
  }
  
  
  highlight_snp <- data[data[[snp_col]] == highlight_snp_id, ]
  
  # --- 2. BUILD THE PLOT ---
  
  p <- ggplot(data = plot_data, aes(x = .data[[pos_col]], y = 1)) +
    # This layer's color now depends on the value (0 or 1)
    geom_point(
      aes(color = .data[[annotation_col]]), # Color is now an aesthetic
      shape = "|", 
      size = 5, 
      alpha = 0.7,
      show.legend = FALSE # Hide the color legend
    ) +
    
    # Manually set the colors for '0' and '1'
    scale_color_manual(values = c("1" = track_color)) +
    
    # Layer to Highlight the specified SNP
    geom_point(
      data = highlight_snp,
      color = "red",
      size = 1.6
    ) +
    
    # Set the genomic position axis
    scale_x_continuous(
      limits = pos_rng,
      labels = function(x) paste0(x / 1e6, " MB")
    ) +
    
    # NEW: Force the y-axis to only show labels for 0 and 1
    scale_y_discrete(breaks = c(0, 1)) +
    
    labs(
      y = annotation_col, 
      x = chromosome_label
    ) +
    theme_bw()
  
  if (!is.null(custom_theme) && is.function(custom_theme)) {
    p <- p + custom_theme()
  }
  
  return(p)
}

p_coding = plot_annotation_highlight(candidate_region, "Coding_UCSC_common", "rs3184504",custom_theme = custom_theme() ) + custom_theme() +    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x =   element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

p_H3K27ac= plot_annotation_highlight(candidate_region, "H3K27ac_Hnisz_common", "rs3184504",custom_theme = custom_theme() ) + custom_theme()+    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x =   element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )


combined_plot<-(p_EUR/p_coding/p_H3K27ac/p_INFUSE/p2+plot_layout(heights = c(1,0.1,0.1,1,1))|p_AFR/p_coding/p_H3K27ac/p_MESuSiE/p2+plot_layout(heights = c(1,0.1,0.1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")

ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.pdf"),combined_plot,width = 7.42, height = 210,units="mm",dpi=300)

plot_dir = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster_ashg_2025/"
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.png"),combined_plot,width = 11.02, height = 7.42,dpi=600)
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.pdf"),combined_plot,width = 11.02, height = 7.42,dpi=600)


###################################################################
#
#
#             SBP Missense and eQTL of FGF5 shared signal, https://pmc.ncbi.nlm.nih.gov/articles/PMC5911643/
#             https://academic.oup.com/edrv/article/40/3/825/5257798
#
##################################################################

drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)


LD_BLOCK = 50
trait = "SBP"

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")

#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK-38,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK-38,".ld")))
rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs1458038 is the eQTL of SH2B3 gene, highlighted in the paper
lead_SNP = "rs1458038"
lead_SNP_index<-which(candidate_region$SNP==lead_SNP)
candidate_region<-candidate_region%>%mutate(r2_EUR  = unname(unlist((EU_cov[,lead_SNP_index])^2)) ,r2_AFR = unname(unlist((BB_cov[,lead_SNP_index])^2)),POS = as.numeric(POS))
####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(susie_EU>0.5&susie_BB>0.5 ~ 3,
                                                                  susie_EU>0.5&susie_BB<0.5 ~ 1,
                                                                  susie_EU<0.5&susie_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            Paintor_cat = case_when(Paintor_PIP>0.5~4,
                                                                    TRUE~0),
                                            Paintor_annot_cat = case_when(Paintor_PIP_all>0.5~4,
                                                                          TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            INFUSE_cat = case_when(MESuSiE_PIP_WB_all_mlk>0.5~1,
                                                                  MESuSiE_PIP_BB_all_mlk>0.5~2,
                                                                  MESuSiE_PIP_Shared_all_mlk>0.5~3,
                                                                  TRUE~0))

###Stat of the region
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs1458038")%>%summarise(2*pnorm(-abs(zscore_WB)),2*pnorm(-abs(zscore_BB)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(zscore_BB)) )%>%filter(Pvalue_BB<5e-8)
candidate_region%>%filter(SNP=="rs1458038")%>%dplyr::select(MESuSiE_PIP_Shared,MESuSiE_PIP_WB,MESuSiE_PIP_BB,SuSiE_Shared,susie_EU,susie_BB,Paintor_PIP)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_WB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(zscore_WB))))

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(zscore_BB))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise((2*pnorm(-abs(zscore_BB))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%pull(zscore_BB)

###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(zscore_WB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(zscore_BB))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_annot_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_all,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_annot_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
INFUSE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either_all_mlk,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(INFUSE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_Paintor_annot<-finemap_plot_fun(Paintor_annot_plot_data, "Paintor+annot", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_INFUSE<-finemap_plot_fun(INFUSE_plot_data, "INFUSE", "PIP", 0.5)

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_Paintor/p_MESuSiE/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_Paintor_annot/p_INFUSE/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)