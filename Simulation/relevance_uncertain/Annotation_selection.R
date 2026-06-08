###First round is just selecting the top annotations we get from polyfun, 
###The mutual correlation between annotations are less than 0.3 (calculate across whole genome not per loci)
###The improvement is not significant with these annotations
load("../../PAINTOR_top_annotation.Rdata")
# Combine all elements into a single vector
all_elements <- unlist(top_annotations_list)

# Count occurrences of each unique element
element_counts <- sort(table(all_elements), decreasing = T)

annot_name  = names(head(element_counts,16))
annot_name


###Second round, we still select the top annotations, but with only binary annotations

test_annot = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_100.annot")
idx = unname(apply(test_annot,2,function(x) {all(x %in% 0:1)}))
test_annot_binary = colnames(test_annot[,idx,with = F])

ori_name = c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common")
annot_name_updated  = head(test_annot_binary[which(names(head(element_counts,50)) %in% test_annot_binary)],20)
annot_name_updated = head(unique(c(ori_name, annot_name_updated)),20)
annot_name_updated

save(annot_name_updated, file = "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/select_annot.Rdata")


