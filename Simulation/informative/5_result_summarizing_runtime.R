source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/utility.R")
simulation_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/"
plot_dir<-paste0(simulation_dir,"plot/")
system(paste0("mkdir -p ",plot_dir))

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/res_summary/")
system(paste0("mkdir -p ",res_dir))
#res_out<-"shared_50_more_annotes_mlk_round2_w2_startover_functional_pc_bic_modified_llk.RData"
res_out<-"shared_50_mlk_round2_w2_startover_runtime.RData"


library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)

num_BB = 300000
num_EU = 300000
################################################################
#
#        Result Summarizing
#
#
################################################################
ancestry_all_ROC_data_list<-list()
either_all_ROC_data_list<-list()
shared_all_ROC_data_list<-list()
all_Set_data_list<-list()
pip_res_all_list<-list()
time_all<-c()
data_all<-c()
empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
MESuSiE_ROC_Either_All<-SuSiE_ROC_Either_All<-MESuSiE_fun_ROC_Either_All<-Paintor_ROC_Either_All<-Paintor_fun_ROC_Either_All<-empty_mat





compare_table = data.frame(causal_num = NA,h2=NA,LOCI_num = NA,Category = NA)

for(causal_num in 1:3){
  for(h2 in 1:2){
    wrk_dir<-paste0(simulation_dir,"causal_num_",causal_num,"/")
    data_dir<-paste0(wrk_dir,"summary_data/w2_startover/")
    result_dir<-paste0(wrk_dir,"result/")
    
    setwd(result_dir)
    
    set = which(paste0("Res_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
    paintor_set = which(paste0("CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".mcmc.paintor")%in%list.files())
    fun_paintor_set = which(paste0("CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".fun.mcmc.paintor")%in%list.files())
    set = Reduce(intersect, list(set,paintor_set,fun_paintor_set))
    empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
    MESuSiE_ROC<-SuSiE_ROC<- MESuSiE_fun_ROC<-Paintor_ROC<-Paintor_fun_ROC<-empty_mat
    MESuSiE_ROC_Either<-SuSiE_ROC_Either<-MESuSiE_fun_ROC_Either<-Paintor_ROC_Either<-Paintor_fun_ROC_Either<-empty_mat
    MESuSiE_ROC_EU<-MESuSiE_ROC_BB<-SuSiE_ROC_EU<-Paintor_ROC_EU<-Paintor_fun_ROC_EU<-empty_mat
    SuSiE_ROC_BB<-MESuSiE_fun_ROC_EU<-MESuSiE_fun_ROC_BB<-Paintor_ROC_BB<-Paintor_fun_ROC_BB<-empty_mat
    num_shared_region<-num_EUR_region<-num_AFR_region<-0
    set_res_all<-c()
    pip_res_all<-c()  
    
    for(LOCI_num in set){
      paintor_res_name<-paste0(result_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".mcmc.paintor")
      paintor_res<-fread(paintor_res_name)
      paintor_fun_res_name<-paste0(result_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".fun.mcmc.paintor")
      paintor_fun_res<-fread(paintor_fun_res_name)
      SuSiE_res_name<-paste0(result_dir,"Res_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      load(SuSiE_res_name)
      zfile<-read.table(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2),header=T)%>%dplyr::rename(SNP = RSID)
      
      #Set Size/Power
      
      MESuSiE_SET = compute_set_susie(unique(unlist(MESuSiE_res$cs$cs)),ifelse(zfile$Signal!=0,1,0))		
      SuSiE_SET = compute_set_susie(unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs))),ifelse(zfile$Signal!=0,1,0))	
      MESuSiE_fun_SET = compute_set_susie(unique(unlist(MESuSiE_res_mlk_annot_subset$cs$cs)),ifelse(zfile$Signal!=0,1,0))
      #Case 1: Great or equal power, lower or equal FDR 
      #Case 2: Great or equal power, greater FDR
      #Case 3: Lower power,  lower or equal FDR  
      #Case 4: Lower power,  greater FDR
      if (MESuSiE_fun_SET[1] >= MESuSiE_SET[1] & MESuSiE_fun_SET[2] <= MESuSiE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 1)
      } else if (MESuSiE_fun_SET[1] >= MESuSiE_SET[1] & MESuSiE_fun_SET[2] > MESuSiE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 2)
      } else if (MESuSiE_fun_SET[1] < MESuSiE_SET[1] & MESuSiE_fun_SET[2] <= MESuSiE_SET[2]) {
        locus_summ = c(causal_num, h2, LOCI_num, 3)
      } else { 
        locus_summ = c(causal_num, h2, LOCI_num, 4)
      }
      compare_table = rbind(compare_table,locus_summ)
      
      
      
      Paintor_SET = compute_set_pip(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob,"CAUSAL"=ifelse(zfile$Signal!=0,1,0)))
      Paintor_fun_SET = compute_set_pip(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob,"CAUSAL"=ifelse(zfile$Signal!=0,1,0)))
      set_res_all<-rbind(set_res_all,rbind(c(MESuSiE_SET,"MESuSiE"),c( SuSiE_SET,"SuSiE"),c(MESuSiE_fun_SET,"INFUSE")),c(Paintor_SET,"Paintor"),c(Paintor_fun_SET,"Paintor+annot"))
      
      #ROC for shared/either/ancestry-specific
      if(sum(zfile$Signal==3)>0){
        num_shared_region<-num_shared_region+1
        MESuSiE_ROC = MESuSiE_ROC + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res$pip_config[,3]),ifelse(zfile$Signal==3,1,0))
        SuSiE_ROC = SuSiE_ROC + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = pmin(susie_EU$pip,susie_BB$pip)),ifelse(zfile$Signal==3,1,0))
        MESuSiE_fun_ROC = MESuSiE_fun_ROC + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res_mlk_annot_subset$pip_config[,3]),ifelse(zfile$Signal==3,1,0))
        Paintor_ROC = Paintor_ROC + compute_ROC(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob),ifelse(zfile$Signal==3,1,0))
        Paintor_fun_ROC = Paintor_fun_ROC + compute_ROC(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob),ifelse(zfile$Signal==3,1,0))
      }
      
      MESuSiE_ROC_Either = MESuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res$pip),ifelse(zfile$Signal!=0,1,0))
      SuSiE_ROC_Either = SuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = pmax(susie_EU$pip,susie_BB$pip)),ifelse(zfile$Signal!=0,1,0))
      MESuSiE_fun_ROC_Either = MESuSiE_fun_ROC_Either + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res_mlk_annot_subset$pip),ifelse(zfile$Signal!=0,1,0))
      Paintor_ROC_Either = Paintor_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob),ifelse(zfile$Signal!=0,1,0))
      Paintor_fun_ROC_Either = Paintor_fun_ROC_Either + compute_ROC(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob),ifelse(zfile$Signal!=0,1,0))
      
      MESuSiE_ROC_Either_All = MESuSiE_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res$pip),ifelse(zfile$Signal!=0,1,0))
      SuSiE_ROC_Either_All = SuSiE_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = pmax(susie_EU$pip,susie_BB$pip)),ifelse(zfile$Signal!=0,1,0))
      MESuSiE_fun_ROC_Either_All = MESuSiE_fun_ROC_Either_All + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res_mlk_annot_subset$pip),ifelse(zfile$Signal!=0,1,0))
      Paintor_ROC_Either_All = Paintor_ROC_Either_All + compute_ROC(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob),ifelse(zfile$Signal!=0,1,0))
      Paintor_fun_ROC_Either_All = Paintor_fun_ROC_Either_All + compute_ROC(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob),ifelse(zfile$Signal!=0,1,0))
      
      if(sum(zfile$Signal==1)>0){
        num_EUR_region<-num_EUR_region+1
        MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res$pip_config[,1]),ifelse(zfile$Signal==1,1,0))
        SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC_uni_susie(data.frame("SNP" =zfile$SNP ,"PIP_1" =susie_EU$pip ,"PIP_2" = susie_BB$pip),ifelse(zfile$Signal==1,1,0))
        MESuSiE_fun_ROC_EU = MESuSiE_fun_ROC_EU + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res_mlk_annot_subset$pip_config[,1]),ifelse(zfile$Signal==1,1,0))
        Paintor_ROC_EU = Paintor_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob),ifelse(zfile$Signal==1,1,0))
        Paintor_fun_ROC_EU = Paintor_fun_ROC_EU + compute_ROC(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob),ifelse(zfile$Signal==1,1,0))
      }
      
      if(sum(zfile$Signal==2)>0){
        num_AFR_region<-num_AFR_region+1
        MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res$pip_config[,2]),ifelse(zfile$Signal==2,1,0))
        SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC_uni_susie(data.frame("SNP" =zfile$SNP ,"PIP_2" =susie_EU$pip ,"PIP_1" = susie_BB$pip),ifelse(zfile$Signal==2,1,0))  
        MESuSiE_fun_ROC_BB = MESuSiE_fun_ROC_BB + compute_ROC(data.frame("SNP" =zfile$SNP ,"PIP" = MESuSiE_res_mlk_annot_subset$pip_config[,2]),ifelse(zfile$Signal==2,1,0))
        Paintor_ROC_BB = Paintor_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$RSID ,"PIP" = paintor_res$Posterior_Prob),ifelse(zfile$Signal==2,1,0))
        Paintor_fun_ROC_BB = Paintor_fun_ROC_BB + compute_ROC(data.frame("SNP" =paintor_fun_res$RSID ,"PIP" = paintor_fun_res$Posterior_Prob),ifelse(zfile$Signal==2,1,0))
        
      }
      
      
      time_out_name<-paste0(result_dir,"time_susie_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".txt")
      time_paintor = fread(paste0(result_dir,"time_paintor_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".txt")) %>% pull(V1)
      
      time_vec<-c(t(fread(time_out_name)),time_paintor,nrow(paintor_res))
      time_all<-rbind(time_all,time_vec)
      
      cat(LOCI_num)
    }
    
    #######################################
    #
    #   credible set combined 
    #
    ####################################### 
    set_res_all_data = data.frame("Power" = as.numeric(set_res_all[,1]),"FDR" = as.numeric(set_res_all[,2]),"Size"=as.numeric(set_res_all[,4]),"Method"=set_res_all[,5])
    set_res_all_data$h2 = h2
    set_res_all_data$causal_num = causal_num
    set_res_summary<-set_res_all_data %>% group_by(Method)%>%summarise_at(vars(Power,Size),list(name = mean)) 
    
    set_res_summary$h2 = h2
    set_res_summary$causal_num = causal_num
    all_Set_data_list[[paste0("causal_",causal_num,"_h2_",h2)]]<-set_res_all_data		
    
    #######################################
    #
    #   ROC combined 
    #
    ####################################### 
    # Either signal part
    roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either, 
                          SuSiE_ROC_Either = SuSiE_ROC_Either, 
                          MESuSiE_fun_ROC_Either = MESuSiE_fun_ROC_Either,
                          Paintor_ROC_Either = Paintor_ROC_Either,
                          Paintor_fun_ROC_Either = Paintor_fun_ROC_Either)
    method_names <- c("MESuSiE", "SuSiE", "INFUSE","Paintor","Paintor+annot")
    either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, length(set), h2, causal_num)
    
    # Shared signal part
    roc_data_list <- list(MESuSiE_ROC = MESuSiE_ROC, SuSiE_ROC = SuSiE_ROC, MESuSiE_fun_ROC = MESuSiE_fun_ROC, Paintor_ROC = Paintor_ROC, Paintor_fun_ROC = Paintor_fun_ROC)
    method_names <- c("MESuSiE", "SuSiE", "INFUSE","Paintor","Paintor+annot")
    shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, num_shared_region, h2, causal_num)
    
    # Ancestry specific part
    roc_data_list <- list(MESuSiE_ROC_EU = MESuSiE_ROC_EU,MESuSiE_ROC_BB = MESuSiE_ROC_BB, 
                          SuSiE_ROC_EU = SuSiE_ROC_EU, SuSiE_ROC_BB = SuSiE_ROC_BB,
                          MESuSiE_fun_ROC_EU = MESuSiE_fun_ROC_EU,MESuSiE_fun_ROC_BB = MESuSiE_fun_ROC_BB,
                          Paintor_ROC_EU = Paintor_ROC_EU, Paintor_ROC_BB = Paintor_ROC_BB,
                          Paintor_fun_ROC_EU = Paintor_fun_ROC_EU, Paintor_fun_ROC_BB = Paintor_fun_ROC_BB)
    method_names <- c("MESuSiE_EU", "MESuSiE_BB","SuSiE_EU", "SuSiE_BB","INFUSE_EU", "INFUSE_BB","Paintor_EU","Paintor_BB","Paintor+annot_EU","Paintor+annot_BB")
    ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, length(set), h2, causal_num)
    
    
  }
}



#######################################
#
#  Loci performance summary
#
####################################### 
compare_table = compare_table[-1,]
summary_table <- compare_table %>%
  group_by(Category) %>%
  summarise(
    Count = n()  # Count of rows per categor
  )

summary_table
#######################################
#
#  Set Size and Power
#
####################################### 

all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
all_Set_data_dataframe<-all_Set_data_dataframe%>%mutate(Method = fct_relevel(Method, c("SuSiE","Paintor","Paintor+annot","MESuSiE", "INFUSE")))
#levels(all_Set_data_dataframe$Method)<-c("MeSuSiE","MESuSiE_fun_","SuSiE")
cols <- c("Method","h2","causal_num")
set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
colnames(set_power_summary)[4]<-"Power_name"

set_power_summary<-set_power_summary%>%mutate(Method = fct_relevel(Method, c("SuSiE","Paintor","Paintor+annot","MESuSiE", "INFUSE")))
#levels(set_power_summary$Method)<-c("MESuSiE","SuSiE","MESuSiE_fun_")
set_power_summary
set_size_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars( Size),list(name = mean)))
colnames(set_size_summary)[4]<-" Size_name"
set_size_summary

#######################################
#
#   ROC data 
#
####################################### 
shared_all_ROC_data_dataframe<-ROC_list_to_data(shared_all_ROC_data_list)
either_all_ROC_data_dataframe<-ROC_list_to_data(either_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ROC_list_to_data(ancestry_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ancestry_all_ROC_data_dataframe%>%mutate(Method = fct_recode(Method,
                                                                                              "SuSiE BB" = "SuSiE_BB",
                                                                                              "SuSiE WB" = "SuSiE_EU",
                                                                                              "Paintor BB" = "Paintor_BB",
                                                                                              "Paintor WB" = "Paintor_EU",
                                                                                              "Paintor+annot BB" = "Paintor+annot_BB",
                                                                                              "Paintor+annot WB" = "Paintor+annot_EU",
                                                                                              "MESuSiE BB" = "MESuSiE_BB",
                                                                                              "MESuSiE WB" = "MESuSiE_EU",
                                                                                              "INFUSE BB" = "INFUSE_BB",
                                                                                              "INFUSE WB" = "INFUSE_EU"))

########################################################################
#
#
#FDR Power based on threshold of 0.01,0.05,0.1,0.5
#
#
########################################################################

FDR_Power_shared<-FDR_Power(shared_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "INFUSE","Paintor","Paintor+annot")))
FDR_Power_shared%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))


FDR_Power_either<-FDR_Power(either_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "INFUSE","Paintor","Paintor+annot")))
FDR_Power_either%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

FDR_Power_ancestry<-FDR_Power(ancestry_all_ROC_data_dataframe)%>%
  mutate(Method = fct_relevel(Method, "SuSiE BB" ,"SuSiE WB" ,"MESuSiE BB","MESuSiE WB","INFUSE BB","INFUSE WB","Paintor BB","Paintor WB","Paintor+annot BB","Paintor+annot WB" ))
FDR_Power_ancestry%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))


#######################################
#
#   time data 
#
####################################### 
time_all<-data.frame(time_all)
colnames(time_all)<-c("MESuSiE","SuSiE","INFUSE","Paintor","Paintor+annot","NUM_SNP")


#######################################
#
#   save the result
#
####################################### 
save(either_all_ROC_data_dataframe,shared_all_ROC_data_dataframe,ancestry_all_ROC_data_dataframe,FDR_Power_either,FDR_Power_shared,FDR_Power_ancestry,all_Set_data_dataframe,set_power_summary,time_all,data_all,file = paste0(res_dir,res_out))

################################################################
#
#        Numbers to report
#
#
################################################################


##95% Credible set size/power/fdr
all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))
all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))


#Power with a FDR = 0.05 
FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))



FDR_Power_either%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.1)%>%group_by(Method)%>%summarise(round(mean(Power),2))



#FDR and Power with a PIP threshold of 0.5	
shared_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
ancestry_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))


################################################################
#
#        Plot
#
#
################################################################
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/utility.R")
upper_limit<-round(all_Set_data_dataframe%>%filter(Method == "Paintor",causal_num=="Num~Causal  == 5 ",h2=="~h^2 == 10^-4")%>%summarise(upper = quantile(Size,0.75))%>%pull(upper))
p_size_box<-Set_Size_fun(all_Set_data_dataframe%>%mutate(Size= log2(Size+1)),upper_limit =log2(upper_limit+1))
p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")
p_size_box

p_power_bar<-Set_Power_fun(set_power_summary)
size_power<-p_size_box/p_power_bar+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 7,face="bold"))

ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=180, width=270, units = "mm",dpi=500)
##########################################################
#
#   Either Signal 
#  ROC | FDR Power 
#
##########################################################

###################
#
#ROC
#
###################
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE"))
p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################

p_FDR_Power_Either<-FDR_Power_shared_fun(FDR_Power_either%>%filter(FDR!=0.01)%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE")))+ylim(0,min(ceiling(max(FDR_Power_either$Power)*10)/10+0.1,1))


ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Either.pdf"),ROC_FDR_Power_Either,height=270, width=270, units = "mm",dpi=500)
##########################################################
#
# Shared Signal 
#  ROC | FDR Power 
#
##########################################################
###################
#
#ROC
# 
###################
shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE"))
p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################
p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared%>%filter(FDR!=0.01)%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE")))+ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1))

ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=270, width=270, units = "mm",dpi=500)


##########################################################
#
# Ancestry-specific Signal 
#  ROC | FDR Power | PIP calibration  
#
##########################################################
###################
#
#ROC
#
###################

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>% 
  mutate(Method = as.character(Method))
split_list <- strsplit(ancestry_all_ROC_data_dataframe %>% pull(Method), " +")
ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)	%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

p_ROC_ancestry<-ROC_ancestry_fun(ancestry_all_ROC_data_dataframe)

###################
#
#FDR&Power
#
###################

FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method,"SuSiE","Paintor","Paintor+annot","MESuSiE","INFUSE"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))


p_FDR_Power_ancestry<-FDR_Power_ancestry_fun(FDR_Power_ancestry%>%filter(FDR!=0.01))+ylim(0, min(ceiling(max(FDR_Power_ancestry$Power)*10)/10+0.1,1))
ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(2, 1))
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=270, width=270, units = "mm",dpi=500)









##############################################################################################################################################################




##########################################################
#
#
#             Summarized result
#
#
##########################################################

###Time Plot
#load(paste0(res_dir,"shared_50_baseline.RData"))
load(paste0(res_dir,res_out))

time_all_long <- time_all %>%
  pivot_longer(
    cols = 1:5,
    names_to = "Method",
    values_to = "Time"
  ) %>%
  mutate(
    Method = recode(Method, "INFUSE" = "INFUSE"),
    Method = factor(
      Method,
      levels = c("SuSiE",  "Paintor", "Paintor+annot","MESuSiE", "INFUSE")
    )
  )


runtime_summary <- time_all_long %>%
  group_by(Method) %>%
  summarise(
    n_regions = n(),
    time_mean = mean(Time, na.rm = TRUE),
    time_median = median(Time, na.rm = TRUE),
    time_q25 = quantile(Time, 0.25, na.rm = TRUE),
    time_q75 = quantile(Time, 0.75, na.rm = TRUE),
    time_p95 = quantile(Time, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Method)

runtime_summary_report <- runtime_summary %>%
  mutate(
    across(
      starts_with("time_"),
      ~ round(.x, 4)
    )
  )

runtime_summary_report

write.csv(
  runtime_summary_report,
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/Supplementary_files/supp_figures/runtime/runtime_summary_by_method.csv",
  row.names = FALSE
)

method_colors <- c(
  "SuSiE" = "#fbe79e",
  "Paintor" = "#f2b56e",
  "Paintor+annot" = "#f57c6e",
  "MESuSiE" = "#84c3b7",
  "INFUSE" = "#71b8ed"
)

custom_theme_half <- function() {
  theme(
    axis.text.x = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.5),
    axis.title.x = element_text(size = 7.5, face = "bold"),
    axis.title.y = element_text(size = 7.5, face = "bold"),
    strip.text.x = element_text(size = 5.5),
    strip.text.y = element_text(size = 5.5),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6.5, face = "bold"),
    plot.title = element_text(size = 7.5, hjust = 0.5),
    plot.tag = element_text(size = 7.5, face = "bold")
  )
}

p_time_1 <- ggplot(
  time_all_long %>% filter(Time < 300),
  aes(x = NUM_SNP, y = Time, color = Method)
) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "Number of SNPs",
    y = "Time (minutes)",
    color = NULL
  ) +
  theme_bw() +
  custom_theme() +
  custom_theme_half() +
  theme(
    legend.position = "bottom"
  )



p_time_1_low = p_time_1

save(
  p_time_1_low,
  file = "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/Supplementary_files/supp_figures/runtime/p_time_1_low.RData"
)


p_size_box = ggplot(data=all_Set_data_dataframe, aes(x = Method, y =log2(Size+1), fill=Method)) +
  geom_boxplot(aes(fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray") +
  scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed"), guide = FALSE) +
  ylab("log2(Set Size + 1)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)+ custom_theme() 



p_power_bar = ggplot(data=set_power_summary%>%group_by(Method)%>%summarize(Power_name= mean(Power_name)),aes(x = Method, y =Power_name,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+
  scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed"),guide=FALSE)+
  ylab("Power")+xlab("")
p_power_bar = p_power_bar+geom_text(aes(x=Method,y=Power_name,label=round(Power_name,2)),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)
p_power_bar = p_power_bar +  theme(
  axis.text.x = element_text(size = 7,face="bold"))+ 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)+theme_bw()+custom_theme()



roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either_All, 
                      SuSiE_ROC_Either = SuSiE_ROC_Either_All, 
                      MESuSiE_fun_ROC_Either = MESuSiE_fun_ROC_Either_All,
                      Paintor_ROC_Either = Paintor_ROC_Either_All,
                      Paintor_fun_ROC_Either = Paintor_fun_ROC_Either_All)
method_names <- c("MESuSiE", "SuSiE", "INFUSE","Paintor","Paintor+annot")
either_all_ROC_data_dataframe_mean<- create_data_frame(roc_data_list, method_names, 600, h2, causal_num)

#	either_all_ROC_data_dataframe_mean<-either_all_ROC_data_dataframe%>%filter(h2=="~h^2 == 2%*%10^-4",causal_num == "Num~Causal  == 3 ")
p_roc_either <- ggplot(either_all_ROC_data_dataframe_mean, aes(x = Power, y = 1 - FDR)) + 
  geom_line(aes(color = Method), size = 0.5) +
  scale_color_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed")) + 
  xlab("Recall") + 
  ylab("Precision") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  geom_point(data = either_all_ROC_data_dataframe_mean %>% filter(Cutoff == 0.5), 
             aes(x = Power, y = 1 - FDR, group = Method, color = Method, shape = 'PIP > 0.5'), 
             size = 2.5, show.legend = c(color = FALSE)) +
  scale_shape_manual(name = ' ', values = c('PIP > 0.5' = 21)) + guides(color = guide_legend(order = 1), 
                                                                        fill = guide_legend(order = 1), 
                                                                        shape = guide_legend(order = 2))+custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)





FDR_Power_Combined<-rbind( FDR_Power_either%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Either"),
                           FDR_Power_shared%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Shared"),
                           FDR_Power_ancestry%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.05)%>%mutate(Cat = "Ancestry"))




FDR_Power_Combined$Method <- factor(
  FDR_Power_Combined$Method,
  levels = c(
    "SuSiE",
    "Paintor",
    "Paintor+annot",
    "MESuSiE",
    "INFUSE"
  )
)


FDR_Power_Combined$Cat<-factor(	 FDR_Power_Combined$Cat,levels = c("Either","Shared","Ancestry"))
#FDR_Power_Combined$Method<-factor(FDR_Power_Combined$Method,levels = c("SuSiE","Paintor","Paintor_fun","MESuSiE","MESuSiE_fun"))
p_FDR_Power<-ggplot(  FDR_Power_Combined, aes(Cat, Power, fill = Method))+
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(x=Cat,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed"))+
  ylab("Power")+xlab("FDR = 0.05")+ ylim(0,0.6)+
  theme_bw() + custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+theme(legend.position = 'bottom')


FDR_Power_Combined_1_percent<-rbind( FDR_Power_either%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                                       filter(FDR==0.1)%>%mutate(Cat = "Either"),
                                     FDR_Power_shared%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                                       filter(FDR==0.1)%>%mutate(Cat = "Shared"),
                                     FDR_Power_ancestry%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.1)%>%mutate(Cat = "Ancestry"))



FDR_Power_Combined_1_percent$Method <- factor(
  FDR_Power_Combined_1_percent$Method,
  levels = c(
    "SuSiE",
    "Paintor",
    "Paintor+annot",
    "MESuSiE",
    "INFUSE"
  )
)
FDR_Power_Combined_1_percent$Cat<-factor(	 FDR_Power_Combined_1_percent$Cat,levels = c("Either","Shared","Ancestry"))
p_FDR_Power_1_percent<-ggplot(  FDR_Power_Combined_1_percent, aes(Cat, Power, fill = Method))+
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(x=Cat,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed"))+
  ylab("Power")+xlab("FDR = 0.05")+ ylim(0,0.6)+
  theme_bw() + custom_theme() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+theme(legend.position = 'bottom')




size_power<-p_size_box/p_power_bar
roc_fdr_power<-p_roc_either/p_FDR_Power+theme(legend.position = 'bottom')
p_total <- (size_power | roc_fdr_power) + 
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides = "collect") & 
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

p_total
save(p_size_box,p_power_bar,p_roc_either,p_FDR_Power,p_FDR_Power_1_percent,file = paste0(plot_dir,"summarize_high_dim_w2_baseline.Rdata"))
ggsave(paste0(plot_dir,"summarize.pdf"),p_total,height=200, width=300, units = "mm",dpi=500)

p_total <- (size_power | roc_fdr_power)  +  
  plot_layout(guides = "collect") & 
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

p_total
ggsave(paste0(plot_dir,"summarize_notag.pdf"),p_total,height=140, width=180, units = "mm",dpi=500)
