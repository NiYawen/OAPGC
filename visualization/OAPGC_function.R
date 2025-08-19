library(ggplot2)
library(gridExtra)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(ggtext)
library(ggrepel)
library(vegan)
library(ggpubr)
library(reshape2)
library(randomForest)
library(tidyselect)
library(pROC)

zy_pie <- function(dt, value, fill, facet_my=NULL, col=NULL,label=T){
  total_color1 = c(brewer.pal(12,"Set3"), brewer.pal(12,"Paired"))
  
  if(typeof(col) == "NULL"){
    if(length(unique(dt[,fill])) > length(total_color1)){
      message("ERROR!!!\n分类太多，请不要超过默认的24个")
      exit(1237)
    }else{
      col = total_color1[1:length(unique(dt[,fill]))]
    }
  }
  
  if (typeof(facet_my) == "NULL"){
    data = dt[, c(fill,value)]
    colnames(data) = c("fill", "value")
    data = data[order(data$value, decreasing=T), ]
    
    data$fill = factor(data$fill, levels=unique(data$fill))
    
    ss2 = sum(data$value)
    plot_dt <- data %>%
      mutate(
        Perc =  round(value/ss2, digits=4) ) %>% 
      mutate(
        label = paste0(fill, " ,", value,", ", round((Perc)*100, digits = 4), "%")
        ,ypos = cumsum(Perc) - 0.5 * Perc
        ,wght=runif(length(fill))
        ,wght=wght/sum(wght)
        ,wght=round(wght, digits=2)
      )
  }
  
  p <- ggplot(plot_dt, aes(x = "", y = Perc , fill = fill)) +
    geom_bar(stat = "identity", width = 1 , color = "white", show.legend = T)+
    coord_polar("y", start = 0) +
    theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.key = element_rect(fill = 'transparent'), 
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    scale_fill_manual(values=col)+
    guides(fill = guide_legend(title = ""))
  
  if(isTRUE(label)){
    P <- p+
      geom_text_repel(aes(x = 1.5,y=1-ypos, label = label),
                      ,color="black"
                      , nudge_x = 0.3
                      ,hjust=0
                      ,size = 3
                      , segment.color = "gray50",
                      segment.linetype="dashed",force=2,
                      min.segment.length = 0,
                      box.padding = 0.5 ,
                      segment.size = .2) 
  }
  p
  
}


alpha_multi_group = function(dt=NA, sample_map=NA, group="Group", 
                                 color.group=NA,ID="Sample", 
                                 index="shannon", 
                                 sample.color=NA, 
                                 box_width=0.5, 
                                 title="alpha diversity", 
                                 violin = F,dir=NA
){
  ## align dt and group
  inter <- intersect(colnames(dt),sample_map[,ID])
  if (length(inter) == 0) {  
    stop("Error: No overlapping columns found between dt and sample_map.")  
  } 
  dt <- dt[,inter]
  dt = dt[rowSums(dt)!=0,]
  sample_map <- filter(sample_map, sample_map[,ID] %in% inter)
  
  #alpha
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }else{
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
  
  if(!is.na(dir)){
    write.csv(dm,paste(dir,".csv",sep=""))
  }
  
  p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[color.group]]))
  
  if(isTRUE(violin)){
    p <- p+
      geom_violin()+
      geom_boxplot(width=box_width, fill="white",
                   position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.colour = NA)
  }else{
    p <- p+ 
      geom_boxplot(position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
  }
  
  ylabs = structure(c("Obs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                    names=c("obs", "shannon", "simpson","invsimpson"))
  ylab = ylabs[tolower(index)]
  
  
  p <- p+theme_bw()+
    theme(strip.text = element_blank())+
    scale_fill_manual(values=sample.color)+
    labs(title=title, y = ylab, x=NULL)
  
  p
}

align_dt_sample <- function(dt, sample_map, ID=NA){
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(.) !=0)
  list(dt=dt, sample_map=sample_map)
}

group_compositions_others <- function(dt=NA, sample_map=NA, ID=NA, group=NA,   
                                         title=NA, taxo.color = NULL, others.color=NA,
                                         width=0.9, label_order=NA,other_abund=NA,
                                         order_func = NA, order_n = 1, group_level=NA){  
  
  if(typeof(sample_map) != "list"){  
    cat("\n\033[1;31m[ERROR!!!]\033[0m\t \033[1msample_map\033[0m is not a list or data.frame\n\n")  
    return ()  
  }  
  
  x = align_dt_sample(dt, sample_map, ID=ID)  
  dt = x$dt  
  sample_map = x$sample_map  
  
  dt = dt[rowSums(dt) != 0,] 
  
  avg_abundance = rowMeans(dt, na.rm = TRUE)  
  dt.dom = dt[avg_abundance >= other_abund,] 
  dt.other=dt[avg_abundance < other_abund,]
  dt.other["others",]=colSums(dt.other)
  dt.all = rbind(dt.dom,dt.other["others",])
  dt = dt.all[order(rowMeans(dt.all), decreasing=T),]
  
  dl = melt(as.matrix(dt))  
  
  if(order_func %in% c("order", "cluster", "specific")){  
    if (order_func == "order"){  
      label_order = dt[order_n,] %>% t() %>% as.data.frame() %>% arrange_all() %>% rownames()  
    } else if(order_func == "cluster"){  
      otu.dist = vegdist(t(dt), method="bray")  
      hc = hclust(otu.dist)  
      label_order = hc$labels[hc$order]  
    } else if(order_func == "specific"){  
      label_order = label_order  
    }  
  }  
  
  tax_ord = rownames(dt)  
  tax_ord = c(setdiff(tax_ord, "others"), "others") 
  taxo.color <- taxo.color[1:length(tax_ord)]
  taxo.color <- c(setdiff(taxo.color,others.color),others.color)
  color.map <- setNames(taxo.color,tax_ord)
  
  dm = merge(dl, sample_map, by.x='Var2', by.y=ID)  
  dm$Var1 = factor(dm$Var1, level=tax_ord)  
  dm$Var2 = factor(dm$Var2, level=label_order)  
  dm[,group]=factor(dm[,group],level=group_level)
  
  p <- ggplot(dm, aes(x=Var2, y=value, fill=Var1))+  
    geom_bar(stat='identity', width=width)+
    theme_classic2()+
    facet_grid(as.formula(paste(". ~", group)), scale = "free", space = "free_x",switch = "both")+
    ggtitle(label=title)+
    scale_fill_manual(values=taxo.color)+
    scale_y_continuous(expand = c(0,0))+
    labs(x=NULL, y="% Relative Abundance")+
    theme(panel.grid = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_text(angle=90, hjust=1),
          legend.position = "bottom")+
    guides(fill = guide_legend(title = ""))
  
  
  p  
}  


disease.diversity.change=function(dt=NA,sample_map=NA, ID="Sample",index="shannon",
                                  Group="Group",Project="Project",sample.color=NA,title=NA,
                                  SampleType="SampleType",DiseaseType="disease_type",
                                  disease_order=disease_order,type_order=NA,dir=NA){
  
  sample_map <- filter(sample_map, !is.na(!!sym(Group)))
  inter <- intersect(colnames(dt),sample_map[,ID])
  if (length(inter) == 0) {  
    stop("Error: No overlapping columns found between dt and sample_map.")  
  } 
  dt <- dt[,colnames(dt) %in% inter]
  dt = dt[rowSums(dt)!=0,]
  sample_map <- filter(sample_map, sample_map[,ID] %in% inter)
  
  
  #alpha
  alpha=data.frame()
  
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }
  
  if(tolower(index) == "shannon"){
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  alpha <- alpha %>% rownames_to_column(ID)
  dm = merge(alpha,sample_map, by=ID)
  
  #wilcox
  vprojtype <- unique(sample_map[,Project])
  res_list = list()
  for (proj in vprojtype){
    #proj=vprojtype[1]
    print(proj)
    tmpf <- dm %>% filter(!!sym(Project) == proj) %>% unique()
    if (length(unique(tmpf[[Group]])) > 1) {
      d1 = filter(tmpf,!!sym(Group) == "Control") %>% pull(alpha)
      d2 = filter(tmpf,!!sym(Group) != "Control") %>% pull(alpha)
      p = wilcox.test(d1,d2)$p.value
      res <- data.frame(Project=proj, pval=p)
      res_list = append(res_list, list(res))
    }
  }
  res_pval <- do.call("rbind", res_list)
  
  #foldchange
  dm[[Group]][dm[[Group]] != "Control"] <- "Disease"
  
  #text_pos
  data.fc <- dm %>%
    dplyr::group_by(!!sym(Project),!!sym(Group)) %>%
    dplyr::summarise(alpha = mean(alpha)) %>%
    pivot_wider(names_from = !!sym(Group),values_from = alpha) %>% 
    dplyr::mutate(fold_change=(Disease-Control)/Control,
                  text_pos = ifelse(fold_change < 0, fold_change-0.06, fold_change+0.06)) %>%
    merge(res_pval, by.x=Project,by.y="Project") %>% 
    left_join(distinct(dm[,c(Project,DiseaseType,SampleType)]),by=Project)
  
  if(!is.na(dir)){
    write.csv(data.fc,paste(dir,"table.csv",sep="."))
  }
  
  
  ###plot
  data.fc[,DiseaseType] <- factor(data.fc[,DiseaseType],levels= disease_order)
  data.fc[,SampleType] <- factor(data.fc[,SampleType],levels=type_order)
  data.fc[,Project] <- factor(data.fc[,Project],levels=rev(project_order))
  
  data.fc <- data.fc %>%
    mutate(shape = ifelse(pval < 0.01, "**",ifelse(pval < 0.05, "*", NA)))
  data.fc$fold_change <- data.fc$fold_change*100
  data.fc$text_pos <- data.fc$text_pos*100
  
  
  alpha_plot <- ggplot(data=data.fc, aes(x=fold_change, y=!!sym(Project), color=!!sym(DiseaseType), fill=!!sym(DiseaseType)))+
    geom_segment(aes(x=0,xend=fold_change, y=!!sym(Project),yend=!!sym(Project)))+
    geom_point(size=3)+
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", fill=NA),
          text=element_text(size=20))+
    geom_text(aes(label=shape, x=text_pos), nudge_y=-0.2,size=10,color="black")+ 
    scale_color_manual(values=color.disease)+
    scale_fill_manual(values=color.disease)+
    ggtitle(paste(title,"(%increase/decrease)",sep="\n"))+
    ylab("")
  
  alpha_plot
}

zy_format_class_name <- function(rf_dt=NA, rf_map=NA, zy_sample=NA){
  row.names(rf_dt) = make.names(row.names(rf_dt))
  colnames(rf_dt) = make.names(colnames(rf_dt))
  rf_map$zy_RF_temp_ID = make.names(rf_map[,zy_sample])
  rf_map[,zy_sample] = make.names(rf_map[,zy_sample])
  ndt = ncol(rf_dt)
  nmap = nrow(rf_map)
  intersect_id = intersect(rf_map$zy_RF_temp_ID, colnames(rf_dt))
  if(nmap != ndt){
    message("rf_dt -> nsample: ",ndt)
    message("rf_map -> nsample: ",nmap)
    message("intersect -> ", length(intersect_id))
  }
  rf_dt = rf_dt[,intersect_id]
  rf_map = rf_map[match(intersect_id, rf_map[,zy_sample]),]
  return(list(rf_dt=rf_dt, rf_map=rf_map))
}


zy_RF_two_class <- function(rf_dt=NA, rf_map=NA, 
                            zy_sample="zy_RF_temp_ID", group=NA, 
                            ntree=999, cross_n = 10,
                            nspecies = NA,
                            seed=123){
  if(! is.na(nspecies)){
    rf_dt = rf_dt[1:nspecies,] %>% filter(rowSums(.) !=0)
  }
  set.seed(seed)
  gs = rf_map %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(value=n()) %>%
    as.data.frame()
  
  g1 <- rf_map %>%
    dplyr::filter(across({{group}})==gs[1,1]) %>%
    dplyr::mutate(rf_temp_cross_n=rep(sample(1:cross_n), gs[1,2]/cross_n+1)[1:gs[1,2]])
  
  g2 <- rf_map %>%
    dplyr::filter(across({{group}})==gs[2,1]) %>%
    dplyr::mutate(rf_temp_cross_n = rep(sample(1:cross_n), gs[2,2]/cross_n+1)[1:gs[2,2]])
  
  rf_map = rbind(g1,g2)
  rf_map[,group] = as.factor(rf_map[,group])
  
  predict_result = list() 
  
  for(i in 1:cross_n){
    cat("\rcorss: ", i, " / ", cross_n)
    test_sample = rf_map[rf_map$rf_temp_cross_n == i,]
    test_dt = rf_dt[,test_sample[,zy_sample]]
    
    train_sample = rf_map[rf_map$rf_temp_cross_n != i,]
    train_dt = rf_dt[,train_sample[,zy_sample]]
    
    fit = randomForest(train_sample[,group]~.,data=t(train_dt), ntree=ntree,importance=TRUE, proximity=TRUE)
    pred = as.data.frame(predict(fit, t(test_dt), type='prob'))
    predict_result = append(predict_result, list(pred)) 
  }
  predict_result <- do.call("rbind", predict_result)
  predict_result = merge(rf_map[,c(zy_sample, group)], predict_result, by.x=zy_sample, by.y="row.names")
  predict_result
}

