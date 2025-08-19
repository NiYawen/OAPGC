source("OAPGC_function.R")

#alpha.change
dt_index <- list(c("sgb","shannon"),
                 c("sgb","obs")) 
for (cp in 1:length(dt_index)){
  m=dt_index[[cp]]
  i=m[[1]]
  j=m[[2]]
  dt.input=eval(parse(text=i))
  
  p<- disease.diversity.change(dt=dt.input,sample_map=meta.disease, ID="Sample",index=j,
                               Group="Group",Project="Project",
                               SampleType="SampleType",DiseaseType="disease_type",
                               sample.color=color.disease,title=paste(i,j,sep="."),
                               disease_order=disease_order,
                               type_order=sampleseq[sampleseq %in% unique(meta.disease.re$SampleType)],
                               dir=paste("03.disease/alpha",i,j,sep="."))
  
  ggsave(paste("03.disease/alpha",i,j,"plot.pdf",sep="."),p,width=12,height=10)
}


#R2
res_list = list()
for(proj in unique(sample_map$Project)){
    #proj= unique(sample_map$Project)[1]
    smpf <- sample_map %>% dplyr::select(Sample, Group, Project) %>% filter(Project == proj) %>% unique() %>% .[complete.cases(.),]
    data <- dt[,smpf$Sample] %>% .[, colSums(.) > 0]
    smpf <- filter(smpf,Sample %in% colnames(data))
   
    mydist = vegdist(t(data), method = "bray")
    formula_terms <- paste0("mydist ~ smpf$Group")    
    ##adonis
    ado = adonis2(as.formula(formula_terms),by="margin")  
    adj_r2 = get_adjusted_r2(ado)
    r2 = ado$R2[1]
    p = ado$`Pr(>F)`[1]
    tmp=data.frame(Project=proj, adj_r2=adj_r2, r2=r2,pval=p)
    
    res_list = append(res_list, list(tmp))
  }
result = do.call("rbind",res_list)
proj_map = sample_map %>% dplyr::select(Project,disease_type,SampleType) %>% unique()
 data.plot = merge(result,proj_map, by="Project")
  
##plot
data.plot = data.plot %>%
  mutate(shape = ifelse(pval < 0.01, "**",ifelse(pval < 0.05, "*", NA))) 
data.plot$disease_type <- factor(data.plot$disease_type,levels= disease_order)
data.plot$SampleType <- factor(data.plot$SampleType,levels=type_order)
data.plot$Project <- factor(data.plot$Project,levels=rev(project_order))
  
adonis_plot <- ggplot(data=data.plot, aes(x=adj_r2, y=Project, color=disease_type, fill=disease_type))+
    geom_segment(aes(x=0,xend=adj_r2, y=Project,yend=Project))+
    geom_point(aes(label=shape, x=adj_r2), size=3)+
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", fill=NA),
          text=element_text(size=20))+
    geom_text(aes(label=shape, x=adj_r2+0.06), nudge_y=-.2,size=10, color="black")+ 
    scale_color_manual(values=color.disease)+
    scale_fill_manual(values=color.disease)


#AUC
res_list=list()
for (proj in unique(sample_map$Project)){
    #proj= unique(sample_map$Project)[22]
    smpf <- sample_map %>% dplyr::select(Sample, Group, Project) %>% filter(Project == proj) %>% unique() %>% .[complete.cases(.),]
    data <- dt[,smpf$Sample] %>% .[, colSums(.) > 0]
    smpf <- filter(smpf,Sample %in% colnames(data))
    message("project -> ",proj,"\t\tsample_num -> ",nrow(smpf))
    
    filelist <- zy_format_class_name(rf_dt=data, rf_map=smpf, zy_sample="Sample")
    data <- filelist[["rf_dt"]]
    smpf <- filelist[["rf_map"]]
    
    predict_result<- zy_RF_two_class(rf_dt=data,rf_map=smpf,zy_sample="Sample",
                                     group="Group",
                                     ntree=999, cross_n = 5,
                                     nspecies = NA,
                                     seed=123)
    
    auc = calc_auc(predict_result, pred="Control", true="Group")$table
    auc$project = proj
    tmp = auc
    res_list = append(res_list, list(tmp))
}
  
result.auc <- do.call("rbind",res_list)
  
###plot
plot.data <- merge(result.auc,sample_map,by.x="project",by.y="Project") %>% dplyr::select(c("low","auc","high","project","disease_type","SampleType")) %>% unique()
colnames(plot.data) = c("low","auc","high","Project","disease_type","SampleType")
plot.data$disease_type <- factor(plot.data$disease_type,levels= disease_order)
plot.data$SampleType <- factor(plot.data$SampleType,levels=type_order)
plot.data$Project <- factor(plot.data$Project,levels=rev(project_order))
  
auc_plot <- ggplot(data=plot.data, aes(x=auc, y=Project, color=disease_type, fill=disease_type))+
    geom_bar(stat='identity',width=0.8)+
    geom_errorbar(aes(xmax=high, xmin=low), color='black',width=0)+
    geom_vline(xintercept = 50, linetype="dashed")+
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", fill=NA),
          text=element_text(size=20))+
    scale_color_manual(values=color.disease)+
    scale_fill_manual(values=color.disease)+
    scale_x_continuous(breaks = scales::pretty_breaks())
  
  
