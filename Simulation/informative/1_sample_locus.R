wrk_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/"
library(data.table)
library(dplyr)
#candidate_locus_clump<-fread("/scratch/negishi/chen4422/hapnest/Simulation_update/define_region/merged_regions_clumping.bed")
#colnames(candidate_locus_clump) = c("CHR","MinBP","MaxBP")


candidate_locus<-fread("/scratch/negishi/chen4422/hapnest/Simulation_update/define_region/merged_regions_1mb.bed")[,c(1:3)]
candidate_locus = candidate_locus[!duplicated(candidate_locus)]
colnames(candidate_locus) = c("CHR","MinBP","MaxBP")
candidate_locus_subset<-candidate_locus%>%group_by(CHR)%>%summarise(NUM_SNP = n())
set.seed(0114)
sample_index<-sample(seq_len(nrow(candidate_locus)),100,replace = F)
#summary(candidate_locus_subset$NUM_SNP[sample_index])
selected_region<-candidate_locus[sample_index]

#candidate_locus_select_region<-candidate_locus%>%filter(Region%in%selected_region)%>%arrange(Region)
#candidate_locus_select_region<-candidate_locus_select_region%>%group_by(Trait,Region)%>%mutate(sim_region = cur_group_id())
selected_region = selected_region %>% arrange(CHR,MinBP)
#write.table(candidate_locus_select_region,paste0(wrk_dir,"simulation_locus_info.txt"),col.names=T,row.names=F,quote=F)
write.table(selected_region,paste0(wrk_dir,"simulation_locus_info.txt"),col.names=T,row.names=F,quote=F)
