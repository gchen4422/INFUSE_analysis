# Load necessary package
library(dplyr)


df = fread("BaselineLF2.2_description.csv")


# (Optional) Rename columns for ease of use
df <- df %>%
  dplyr::rename(Annotation_Name = `Annotation Name`,
         Detailed_Description = `Detailed Description`)


# Assign groups based on keywords in Category or Annotation_Name
df <- df %>%
  mutate(Group = case_when(
    # Group 1: Coding and Gene Structure (coding exons, introns, variant classifications)
    grepl("Coding", Category) | grepl("synonymous", Annotation_Name) | 
      grepl("non_synonymous", Annotation_Name) | grepl("Intron", Category) ~ 1,
    
    # Group 2: Promoter and UTR Annotations (promoters, UTRs, and related Villar/ExAC annotations)
    grepl("Promoter", Category) | grepl("UTR", Category) | 
      grepl("Human_Promoter", Annotation_Name) | grepl("Ancient_Sequence_Age_Human_Promoter", Annotation_Name) |
      grepl("ExAC", Annotation_Name) ~ 2,
    
    # Group 3: Enhancer Annotations (enhancers, weak enhancers, super enhancers, Villar enhancers)
    grepl("Enhancer", Category) | grepl("Human_Enhancer", Annotation_Name) | 
      grepl("Ancient_Sequence_Age_Human_Enhancer", Annotation_Name) ~ 3,
    
    # Group 4: Conservation Metrics (conserved regions and conservation scores)
    grepl("Conserved", Category) | grepl("GERP", Annotation_Name)|grepl("Conserved", Annotation_Name)  ~ 4,
    
    # Group 10: QTL Annotations (QTL MaxCPP metrics)
    grepl("QTL", Category) ~ 10,
    
    # Group 5: Histone Modification Signals (various histone marks)
    grepl("H3K27ac", Annotation_Name) | grepl("H3K4me1", Annotation_Name) | 
      grepl("H3K4me3", Annotation_Name) | grepl("H3K9ac", Annotation_Name) ~ 5,
    
    # Group 6: Chromatin Accessibility and TF Binding (DNase, DGF, CTCF, TFBS, TSS, transcribed regions)
    grepl("CTCF", Category) | grepl("Digital Genomic Footprinting", Category) | 
      grepl("DNase", Category) | grepl("FetalDHS", Annotation_Name) | 
      grepl("TFBS", Category) | grepl("TSS", Annotation_Name) | grepl("Transcr", Category) ~ 6,
    
    # Group 7: Chromatin States and Repression (repressed regions and bivalent flanking regions)
    grepl("Repressed", Category) | grepl("BivFlnk", Annotation_Name) ~ 7,
    
    # Group 8: Allele Frequency and LD Metrics (MAF bins and MAF-adjusted metrics)
    grepl("MAFbin", Category) | grepl("MAF_Adj", Annotation_Name)|grepl("MAFbin", Annotation_Name) ~ 8,
    
    # Group 9: Population Genetics Metrics (recombination rate, nucleotide diversity, background selection, CpG)
    grepl("Recomb_Rate", Annotation_Name) | grepl("Nucleotide_Diversity", Annotation_Name) | 
      grepl("Backgrd_Selection", Annotation_Name) | grepl("CpG_Content", Annotation_Name) ~ 9,
    
    TRUE ~ NA_real_
  ))

table(df$Group, useNA = "always")

# View the resulting data frame with assigned groups
head(df)

save(df,file = "BaselineLF2.2_description_grouped.Rdata")

load("BaselineLF2.2_description_grouped.Rdata")

