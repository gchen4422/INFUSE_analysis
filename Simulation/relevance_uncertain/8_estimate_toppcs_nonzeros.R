library(Rcpp)
library(dplyr)
library(data.table)
library(elasticnet)

set.seed(123)

idx_sample = sample(1:100,9)
idx = c(2,idx_sample)
idx


for (LD_BLOCK in idx){
  h2_num = sample(1:2,1)
  num_causal = sample(1:3,1)
  
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/causal_num_",num_causal,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  system(paste0("mkdir -p ",data_dir))
  result_dir<-paste0(wrk_dir,"result/")
  system(paste0("mkdir -p ",result_dir))
  
  
  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
  LD_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
  LD_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")
  
  
  annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/all_annotes/annot/")  # we take the full set candidate of annotations
  
  
  zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
  summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
  summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )
  summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 
  
  EU_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
  colnames(EU_cov)<-zfile$RSID
  BB_cov<-as.matrix(fread(paste0(LD_dir_afr,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
  colnames(BB_cov)<-zfile$RSID
  R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)
  
  library(Rcpp)
  library(nloptr)
  library(glmnet)
  sourceCpp("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6.cpp")
  source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_utility.R")
  source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/meSuSie_R6_pca.R")
  
  bayes_fac=4
  ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))
  annot_file_subset<-fread(paste0(annot_dir,"/Region_",LD_BLOCK,".annot")) #%>% select(all_of(annot_name_updated))
  
  MESuSiE_res_mlk_annot_subset<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,
                                             ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =F,max_iter =100,
                                             annot = as.matrix(annot_file_subset), 
                                             annot_method = "mlk")
  
  
  
  
  
  
  input_name = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/test_pca_input/test_cs_input_pca",l_index,".Rdata")
  load(input_name)
  
  sds <- apply(annotation, 2, sd)
  annotation_nonconst <- as.matrix(annotation)[, sds > 0]
  annotation_nonconst_scaled <- scale(annotation_nonconst)
  
  
  # Save the scaling attributes for later back-transformation
  scaling_means <- attr(annotation_nonconst_scaled, "scaled:center")
  scaling_scales <- attr(annotation_nonconst_scaled, "scaled:scale")
  # Define grid parameters for k and nonzero loadings per component
  k_values <- 5:10        # Number of components
  nonzero_values <- 5:10  # Number of nonzero loadings per component
  
  # Create a data frame to store grid search results
  results <- expand.grid(k = k_values, p = nonzero_values)
  results$BIC <- NA
  
  
  for(i in 1:nrow(results)) {
    k <- results$k[i]
    nonzero <- results$p[i]
    
    # Run sparse PCA with current grid parameters.
    spca_result <- spca(annotation_nonconst_scaled, 
                        K = k, 
                        type = "predictor", 
                        sparse = "varnum", 
                        para = rep(nonzero, k))
    
    # Compute the scores and reconstruct the data.
    annotation_spca <- annotation_nonconst_scaled %*% spca_result$loadings
    #reconstruction <- annotation_spca %*% t(spca_result$loadings)
    
    # Calculate the Residual Sum of Squares (RSS) for PCA (not used in this BIC)
    # RSS <- sum((annotation_nonconst_scaled - reconstruction)^2)
    
    # Next, compute beta estimates
    inter_weight <- optimize_llk(as.matrix(annotation_spca), alpha, rep(0, ncol(annotation_spca)))
    gamma <- inter_weight 
    V_k <- spca_result$loadings
    beta_std <- V_k %*% gamma
    beta_orig <- beta_std / scaling_scales
    
    n_beta <- length(beta_orig)
    
    # Compute beta_RSS: sum of squared beta estimates
    beta_RSS <- sum(beta_orig^2)
    
    # Compute the modified BIC based on beta
    results$BIC[i] <- n_beta * log(beta_RSS / n_beta) + log(n_beta) * (k * nonzero)
    
    cat("Processed row", i, "\n")
  }
  
  print(results)
  save(results,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/test_pca_input/","CV_results_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num,"_l4.Rdata"))
  # Identify the optimal parameters (k and nonzero) that minimize the BIC
  optimal_params <- results[which.min(results$BIC), ]
  print(optimal_params)
  
}

files <- list.files(pattern = "^CV_results.*\\.Rdata$")
optimal_results <- data.frame()  # Create an empty data frame to store results

for(file in files) {
  load(file)
  optimal_params <- results[which.min(results$BIC), ]
  optimal_params$file <- file   # Optionally add a column with the file name
  optimal_results <- rbind(optimal_results, optimal_params)
}

optimal_results = optimal_results %>% arrange(k,p)

print(optimal_results)

