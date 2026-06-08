################################################
#
#   ROC, Credible set size 
#
#################################################

compute_ROC<-function(pip_data,true_causal){
  do.call(rbind,lapply(seq(1,0.01,-0.01),function(x){
    
    
    TP = sum(pip_data$SNP[pip_data$PIP>x]%in%pip_data$SNP[true_causal!=0])
    FP = sum(pip_data$PIP>x)-TP
    power = TP/sum(true_causal!=0)
    fdr = FP/(FP+TP)
    return(c(x,ifelse(is.na(power),0,power),ifelse(is.na(fdr),0,fdr)))
    
  }))
}

compute_ROC_uni_susie<-function(pip_data,true_causal){
  do.call(rbind,lapply(seq(1,0.01,-0.01),function(x){
    
    TP = sum(pip_data$SNP[pip_data$PIP_1>x&pip_data$PIP_2<x]%in%pip_data$SNP[true_causal!=0])
    FP = sum(pip_data$PIP_1>x&pip_data$PIP_2<x)-TP
    power = TP/sum(true_causal!=0)
    fdr = FP/(FP+TP)
    return(c(x,ifelse(is.na(power),0,power),ifelse(is.na(fdr),0,fdr)))
    
  }))
}


compute_set_pip<-function(pip_data,coverage = 0.95){
  
  pip_data$PIP = pip_data$PIP/sum(pip_data$PIP)
  snps_selected = pip_data[order(pip_data$PIP,decreasing = T),][which(cumsum(pip_data$PIP[order(pip_data$PIP,decreasing = T)]) < coverage),]
  set_length = nrow(snps_selected)
  total_signal = sum(pip_data$CAUSAL)
  identified_signal =sum(snps_selected$CAUSAL) 
  
  if(set_length==0){
    FDR = 0
  }else{
    FDR = (set_length - identified_signal)/set_length
  }
  
  return(Paintor_set = c(identified_signal/total_signal,FDR,identified_signal,set_length))
  
}

compute_set_susie<-function(snps_selected,true_causal,coverage = 0.95){
  
  
  set_length = length(snps_selected)
  total_signal = sum(true_causal)
  identified_signal =sum(true_causal[snps_selected]) 
  
  if(set_length==0){
    FDR = 0
  }else{
    FDR = (set_length - identified_signal)/set_length
  }
  
  if(total_signal == 0){
    power = 0
  }else{
    power = identified_signal/total_signal
  }
  
  
  return(c(power,FDR,identified_signal,set_length))
  
}
#############################################################
#
#
#
#						Plot function
#
#
#
#################################################################
custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),  
    axis.title.x = element_text(size = 13, face="bold"),
    axis.title.y = element_text(size = 13, face="bold"),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    strip.background = element_blank(),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11, face="bold"),
    plot.title = element_text(size = 11, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

##################################################
#
#         Set Size Plot
#
##################################################

Set_Size_fun<-function(all_Set_data_dataframe,upper_limit){
  p_size_box = ggplot(data=all_Set_data_dataframe, aes(x = Method, y =Size, fill=Method)) +
    geom_boxplot(aes(fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray") +
    coord_cartesian(ylim = c(0, upper_limit)) +
    scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed","INFUSE+correct" = "#c984f2"), guide = FALSE) +
    ylab("Set Size") +
    xlab("") +
    facet_grid(h2~causal_num, labeller=label_parsed, scales = "free", space = "free") +
    theme_bw() + custom_theme() +
    theme(axis.text.x = element_text(size = 10,face="bold")) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)
  return(p_size_box)
}

##################################################
#
#         Set Power Plot
#
##################################################
Set_Power_fun<-function(power_summary){
  p_power_bar = ggplot(data=power_summary,aes(x = Method, y =Power_name,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+
    scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed","INFUSE+correct" = "#c984f2"),guide=FALSE)+
    ylab("Power")+xlab("")+facet_grid(h2~causal_num,labeller=label_parsed)+ylim(c(0,1))+theme_bw()
  p_power_bar = p_power_bar+geom_text(aes(x=Method,y=Power_name,label=round(Power_name,2)),position = position_dodge(width = 1),vjust=-0.5,size = 5*10/14)
  p_power_bar = p_power_bar + custom_theme()+  theme(
    axis.text.x = element_text(size = 10,face="bold"))+ 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)
  return(p_power_bar)
}

##################################################
#
#         ROC for either/shared signal detection 
#
##################################################
ROC_shared_fun <- function(all_ROC_data_dataframe) {
  p <- ggplot(all_ROC_data_dataframe, aes(x = Power, y = 1 - FDR)) + 
    geom_line(aes(color = Method), size = 0.5) +
    scale_color_manual(values = c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed","INFUSE+correct" = "#c984f2")) + 
    xlab("Recall") + 
    ylab("Precision") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    facet_grid(h2~causal_num, labeller = label_parsed) + 
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    geom_point(data = all_ROC_data_dataframe %>% filter(Cutoff == 0.5), 
               aes(x = Power, y = 1 - FDR, group = Method, color = Method, shape = 'PIP > 0.5'), 
               size = 2.5, show.legend = c(color = FALSE)) +
    scale_shape_manual(name = ' ', values = c('PIP > 0.5' = 21)) + guides(color = guide_legend(order = 1), 
                                                                          fill = guide_legend(order = 1), 
                                                                          shape = guide_legend(order = 2))+custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  
  # p <- grid.arrange(p, left = textGrob(ylab_title, rot = 90, gp = gpar(fontface = "bold", fontsize=7)))
  
  return(p)
}


###################################################
#
#
#       ROC for ancestry-specific signal detection
#
#
###################################################
ROC_ancestry_fun <- function(all_ROC_data_dataframe) {
  
  p <- ggplot(all_ROC_data_dataframe, aes(x = Power, y = 1 - FDR)) +
    geom_line(
      aes(
        x = Power, 
        y = 1 - FDR,
        linetype = Ancestry,
        color = Method
      ),
      size = 0.8
    ) +
    scale_color_manual(
      values = c(
        "SuSiE" = "#fbe79e",
        "Paintor" = "#f2b56e",
        "Paintor+annot" = "#f57c6e",
        "MESuSiE" = "#84c3b7",
        "INFUSE" = "#71b8ed",
        "INFUSE+correct" = "#c984f2"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "White British" = "dashed", 
        "Black British" = "solid"
      )
    ) +
    xlab("Recall") + 
    ylab("Precision") + 
    theme_bw() + 
    theme(
      legend.position = "bottom"
    ) + 
    facet_grid(h2 ~ causal_num, labeller = label_parsed) + 
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) + 
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    geom_point(
      data = all_ROC_data_dataframe %>% 
        filter(abs(Cutoff - 0.5) < 0.001), 
      aes(
        x = Power, 
        y = 1 - FDR, 
        group = Method, 
        color = Method, 
        shape = "PIP > 0.5"
      ), 
      size = 3, 
      show.legend = c(color = FALSE)
    ) +
    scale_shape_manual(
      name = " ", 
      values = c("PIP > 0.5" = 21)
    ) + 
    guides(
      color = guide_legend(order = 1),  
      linetype = guide_legend(order = 2), 
      shape = guide_legend(order = 3)
    ) + 
    custom_theme() + 
    annotate(
      "rect", 
      xmin = -Inf, 
      xmax = Inf, 
      ymin = -Inf, 
      ymax = -Inf, 
      colour = "black", 
      size = 0.2
    )
  
  return(p)
}



##################################################
#
#         FDR Power for either/shared signal detection 
#
##################################################
FDR_Power_shared_fun<-function(FDR_Power){
  p<-ggplot(FDR_Power, aes(FDR, Power, fill = Method))+
    geom_bar(stat="identity", position = "dodge")+
    geom_text(aes(x=FDR,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 5*6/14)+
    scale_fill_manual(values=c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed","INFUSE+correct" = "#c984f2"))+
    facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+
    ylab("Power")+xlab("FDR")+ ylim(0,0.7)+
    theme_bw() + custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  return(p)
}
###################################################
#
#     FDR Power for ancestry-specific signal detection
#
#####################################################
FDR_Power_ancestry_fun<-function(FDR_Power){
  p =  ggplot(FDR_Power_ancestry, aes(FDR, Power,fill = Method,pattern = Ancestry)) + 
    geom_col_pattern(position = "dodge", pattern_density = .01, pattern_spacing = 0.05,pattern_key_scale_factor = 0.3,pattern_fill ="#DEDBD2",pattern_colour ="ivory3" ) +
    scale_fill_manual(name = "Method",values = c("SuSiE"="#fbe79e","Paintor" = "#f2b56e" ,"Paintor+annot" = "#f57c6e","MESuSiE"="#84c3b7","INFUSE"="#71b8ed","INFUSE+correct" = "#c984f2"))+   
    scale_pattern_manual(values=c("White British"="stripe","Black British"="circle"))+
    facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+
    ylab("Power")+xlab("FDR")+
    theme_bw() + custom_theme() + theme(legend.position = "bottom")+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+
    guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))
  return(p)
}

#######################################################
#
#   PIP calibration for shared signal
#
#########################################################
PIP_calibration_shared_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  return(p)
}

PIP_calibration_shared_byh2_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(vars(h2),vars(Method),labeller=label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  return(p)
}
#######################################################
#
#   PIP calibration for ancestry specific signal
#
#########################################################
PIP_calibration_ancestry_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(Method~causal_num,labeller = label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  return(p)
}
PIP_calibration_ancestry_byh2_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(h2~Method,labeller = label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  return(p)
}
#########################################################################
#
#             PIP scatter plot
#
#########################################################################
PIP_scatter_fun<-function(pip_all_dataframe_long){
  p = ggplot(data=pip_all_dataframe_long,aes(x=MESuSiE,y=PIP)) +geom_point(alpha=0.3)+ geom_point(data=pip_all_dataframe_long[pip_all_dataframe_long$Signal!="None",],aes(x=MESuSiE,y=PIP,color=Signal),size=3)+scale_colour_manual(values=c("#F2C1B6","#ABDB9F", "#6162B0"))+geom_abline(intercept=0,slope=1, col="red")+ylab("PIP")+theme_bw()
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18),
                axis.title.y = element_text( size = 18),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))	
  return(p)
}
#########################################################################
#
#             PIP scatter plot shared only
#
#########################################################################
PIP_scatter_shared_fun<-function(pip_all_dataframe_long){
  p = ggplot(data=pip_all_dataframe_long,aes(x=MESuSiE,y=PIP)) +geom_point(alpha=0.3)+ geom_point(data=pip_all_dataframe_long[pip_all_dataframe_long$Signal!="None",],aes(x=MESuSiE,y=PIP,color=Signal),size=3)+scale_colour_manual(values=c( "#6162B0"))+geom_abline(intercept=0,slope=1, col="red")+ylab("PIP")+theme_bw()
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18),
                axis.title.y = element_text( size = 18),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))	
  return(p)
}



create_data_frame <- function(roc_data_list, method_names, num_region, h2, causal_num) {
  # Create data frames for each method, normalize by the number of replicates, and add method names
  all_data_list <- lapply(1:length(roc_data_list), function(i) {
    data.frame(roc_data_list[[i]] / num_region, Method = method_names[i])
  })
  
  # Combine data frames into one and set column names
  all_data_combine <- do.call(rbind, all_data_list)
  colnames(all_data_combine) <- c("Cutoff", "Power", "FDR", "Method")
  
  # Add h2 and causal_num columns
  all_data_combine$h2 <- h2
  all_data_combine$causal_num <- causal_num
  
  return(all_data_combine)
}

#########################ROC Data ################################
ROC_list_to_data<-function(dat){
  
  dat_frame<-data.frame(Reduce(rbind,dat))
  dat_frame$h2<-factor(dat_frame$h2)
  levels(dat_frame$h2)=c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
  dat_frame$causal_num<-factor(dat_frame$causal_num)
  levels(dat_frame$causal_num)=c(paste0("Num~","Causal "," == 2 "), paste0("Num~","Causal "," == 3 "), paste0("Num~","Causal "," == 5 "))
  dat_frame$Method<-factor(dat_frame$Method)
  return(dat_frame)
}

FDR_Power<-function(dat){
  res<-do.call(rbind,lapply(c(0.01,0.05,0.1,0.5),function(x){
    intermediate_data<-dat %>% group_by( Method,h2,causal_num) %>% filter(abs(FDR - x) == min(abs(FDR - x))) %>%  dplyr::slice(which.max(Power))
    intermediate_data$FDR = x
    return(intermediate_data)
  }))
  res$FDR = factor(res$FDR,levels=c("0.01","0.05","0.1","0.5"))
  res$Power<-round(res$Power,2)
  return(res)
}



###########PIP Calibration Part
create_obs_frq<-function(pip,signal_num,col_names){
  pip_signal = pip%>%mutate(Signal = ifelse(Signal%in%signal_num,1,0))
  pip_signal_long<-melt(pip_signal, id.vars=c("Signal","h2","causal_num"))
  colnames(pip_signal_long)[4:5]<-c("Method","PIP")
  pip_signal_long<-pip_signal_long%>%filter(Method%in%col_names)
  PIP_calibration<-pip_signal_long%>% group_by(Method,causal_num,gr=cut(PIP, breaks= seq(0, 1, by = 0.1)))%>%summarize(mean_pip = mean(PIP),obs_frq = mean(Signal),N_obs= n())
  PIP_calibration$obs_min = PIP_calibration$obs_frq - 1.96*sqrt(PIP_calibration$obs_frq*(1-PIP_calibration$obs_frq)/PIP_calibration$N_obs)
  PIP_calibration$obs_max = PIP_calibration$obs_frq + 1.96*sqrt(PIP_calibration$obs_frq*(1-PIP_calibration$obs_frq)/PIP_calibration$N_obs)
  PIP_calibration$causal_num<-factor(PIP_calibration$causal_num)
  levels(PIP_calibration$causal_num)=c(paste0("Num~","Causal "," == 1 "), paste0("Num~","Causal "," == 3 "), paste0("Num~","Causal "," == 5 "))
  PIP_calibration<-droplevels(PIP_calibration)
  PIP_calibration$Method<-factor(PIP_calibration$Method)
  return(PIP_calibration)
}
create_obs_frq_byh2<-function(pip,signal_num,col_names){
  pip_signal = pip%>%mutate(Signal = ifelse(Signal%in%signal_num,1,0))
  pip_signal_long<-melt(pip_signal, id.vars=c("Signal","h2"))
  colnames(pip_signal_long)<-c("Signal","h2","Method","PIP")
  pip_signal_long<-pip_signal_long%>%filter(Method%in%col_names)
  PIP_calibration<-pip_signal_long%>% group_by(h2,Method,gr=cut(PIP, breaks= seq(0, 1, by = 0.1)))%>%summarize(mean_pip = mean(PIP),obs_frq = mean(Signal),N_obs= n())
  PIP_calibration$obs_min = PIP_calibration$obs_frq - 1.96*sqrt(PIP_calibration$obs_frq*(1-PIP_calibration$obs_frq)/PIP_calibration$N_obs)
  PIP_calibration$obs_max = PIP_calibration$obs_frq + 1.96*sqrt(PIP_calibration$obs_frq*(1-PIP_calibration$obs_frq)/PIP_calibration$N_obs)
  PIP_calibration$Method<-factor(PIP_calibration$Method)
  PIP_calibration$h2<-factor(PIP_calibration$h2)
  levels(PIP_calibration$h2) =c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
  return(PIP_calibration)
}