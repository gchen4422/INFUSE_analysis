#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Publish matched-size results consistently across canonical summary interfaces.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
root<-normalizePath(Sys.getenv("AOU_UNION_ROOT", file.path(.infuse_work, "INFUSE_variant_availability/real_data_AoU_union")),mustWork=TRUE)
s<-file.path(root,"summary");m<-file.path(s,"matched_size")
archive<-file.path(root,"archive/own_cutoff_43044259/summary")
dir.create(archive,recursive=TRUE,showWarnings=FALSE)
z<-fread(file.path(m,"functional_enrichment.tsv"));sz<-fread(file.path(m,"selection_sizes.tsv"))
stopifnot(nrow(z)==180L*uniqueN(z$annotation),all(sz[,.(k=uniqueN(N)),by=.(group,target,INFUSE_threshold)]$k==1L))
for(n in c("functional_enrichment.tsv","summary_data.rds","SUMMARY_README.txt")){
  f<-file.path(s,n);b<-file.path(archive,n)
  if(file.exists(f)&&!file.exists(b))stopifnot(file.copy(f,b,overwrite=FALSE))
}
r<-readRDS(file.path(s,"summary_data.rds"))
r$tables$functional_enrichment<-z
r$evaluation_annotations<-fread(file.path(s,"evaluation_annotation_rows.tsv.gz"))
added<-data.table(annotation=c("Conserved_LindbladToh_common","H3K4me3_Trynka_common"),
 definition=c("Conserved sequence indicator from allele-matched baselineLF186; also used in fitting",
 "H3K4me3 indicator from allele-matched baselineLF186; also used in fitting"))
defs<-as.data.table(r$tables$evaluation_feature_definitions)
r$tables$evaluation_feature_definitions<-rbind(defs[!annotation %in% added$annotation],added)
r$provenance$features<-r$tables$evaluation_feature_definitions$annotation
r$provenance$added_enrichment_annotations<-file.path(s,"reference_enrichment_annotation_provenance.tsv")
r$tables$enrichment_selection_sizes<-sz
pooled_dir<-file.path(m,"pooled_BMI_BP")
if(file.exists(file.path(pooled_dir,"functional_enrichment.tsv"))){
 r$tables$pooled_configuration_enrichment<-fread(file.path(pooled_dir,"functional_enrichment.tsv"))
 r$tables$pooled_configuration_selection_sizes<-fread(file.path(pooled_dir,"selection_sizes.tsv"))
 r$provenance$pooled_configuration_enrichment<-file.path(pooled_dir,"functional_enrichment.tsv")
}
r$provenance$enrichment<-"Matched-size top N derived from INFUSE target PIP > cutoff; each method's own ranking; annotation-specific complete-case selected/background frequency ratio"
r$provenance$enrichment_groups<-"BMI separately; pooled DBP/SBP records before ranking; cross-trait repeats retained"
r$provenance$enrichment_universe<-"Each method's finite target-score records; no artificial zero for unanalysed SuSiEx SNPs"
r$provenance$enrichment_ties<-"descending target score; ascending trait,locus,SNP"
r$provenance$enrichment_membership<-c(ranked_universe=file.path(m,"ranked_universe.tsv.gz"),selected_records=file.path(m,"selected_records.tsv.gz"))
r$provenance$enrichment_updated<-as.character(Sys.time())
# Old common-universe and own-cutoff threshold counts remain valid separate diagnostics.
r$provenance$universe<-"PIP overlap diagnostics use common finite Either scores; matched-size enrichment uses method-specific finite scores (see enrichment_universe)."
stage<-file.path(s,"summary_data_matched_staging.rds")
saveRDS(r,stage);stopifnot(file.rename(stage,file.path(s,"summary_data.rds")))
fwrite(z,file.path(s,"functional_enrichment.tsv"),sep="\t",na="NA")
fwrite(sz,file.path(s,"enrichment_selection_sizes.tsv"),sep="\t")
fwrite(r$tables$evaluation_feature_definitions,file.path(s,"evaluation_feature_definitions.tsv"),sep="\t")
