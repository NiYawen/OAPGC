source("OAPGC_function.R")

#alpha.bar
i_plots <- list()  
desired_combinations <- list(  
  c("sgb", "shannon"),  
  c("sgb", "obs"),  
  c("arg", "obs"),
)  
for (combination in desired_combinations) { 
  #combination= desired_combinations[[1]]
  i <- combination[1]  
  j <- combination[2]  
  
  print(paste("Processing:", i, j))  
  
  dt <- eval(parse(text = i))  
  
  color.map <- color.site
   
  p1 <- alpha_multi_group(dt = dt, sample_map = meta.f, group = "project_number", 
                              color.group="SampleType", ID = "Used_name", 
                              sample.color = color.map, index = j,
                              dir=paste("01.alpha.diversity",i,j,sep = "_")) +  
    labs(title = paste(j, "of", i, sep = " ")) +   
    theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust = 1),
          text=element_text(size=12))
  
  ggsave(paste("01.habitat/01.alpha.diversity/alpha", i, j, "plot.pdf", sep = "_"), p1, height = 8, width = 16)  
  
  i_plots[[length(i_plots) + 1]] <- p1  
}  


#bar.composition
v2 <- c("genus","phylum")
for (i in v2){
  dt = eval(parse(text=i))
  print(i)
  dt1 <- dt %>% mutate(dt=rownames(.))
  dt_long <- melt(dt1,id.vars="dt",variable.name="Used_name",value.name="abundance")
  dfm <- left_join(dt_long,meta.f,by="Used_name")
  projabund <- dfm %>% dplyr::group_by(dt,project_number) %>% dplyr::summarise(mean_abundance = mean(abundance,na.rm=T),.groups="drop")
  dfmean2 <- projabund %>% 
    pivot_wider(names_from = project_number,values_from = mean_abundance) %>%
    column_to_rownames("dt") %>% as.data.frame()
  p2 <- group_compositions_others(dt=dfmean2, sample_map=sample.map, 
                                     ID="project_number",group="SampleType",
                                     taxo.color=color32,others.color="#45bcec",other_abund=0.1,
                                     order_func = "specific", label_order = genus.order, 
                                     group_level=sampleseq)
  p2
  ggsave(paste("02.composition/plot_",i,"_stackbar.pdf"),p2 ,width=17, height=8)
  
}

#scatter plot
sgbs = read.table("./sgb.list", sep="\t")

all_res= rbind()
for(i in sgbs$V1){
  # load(paste(i,"/adonis_patristic.RData", sep=""))
  load(paste(i,"/adonis_ani.RData", sep=""))
  all_res = rbind(all_res, res)
}

x = "R2"



########## all
tax = dtf$taxonomy
tax = do.call('rbind', strsplit(dtf$taxonomy,";"))
colnames(tax) = c("Domain","Phylum","Class","Order","Family","Genus","Species")
pdt = cbind(dtf, tax)
pdt$x = pdt[,x]
pdt = pdt %>%  
  mutate(labs = paste(Species, "(", SGB,")", sep=""),
         sig = ifelse(pvalue<0.01,"+",ifelse(pvalue<0.05,"*",NA)))
ord = pdt %>% 
  arrange(desc(Phylum), desc(vars), x) %>%
  mutate(labs = paste(Species, "(", SGB, ")", sep="" ))

pdt$y = factor(pdt$labs, levels=unique(ord$labs))


max_axis = max(pdt$R2)
min_axis = max(min(pdt$R2), 0)
limt_axis = c(min_axis, max_axis)

mm = dcast(pdt, y + Phylum + total + pvalue ~ vars, value.var="R2") %>%
  mutate(xx = ifelse(is.na(Continent), 0, Continent),
         yy = ifelse(is.na(SampleType), 0, SampleType)) %>%
  mutate(dist = abs(xx - yy)/sqrt(2) )

q3 = quantile(mm$dist, 0.90, na.rm=T)
mm <- mm %>%
  mutate(labs = ifelse(dist > q3, as.character(y), NA) )

p3 <- ggplot(mm, aes(x=xx, y=yy,fill=Phylum))+
  geom_point(aes(size=total),shape=21, color="black", alpha=0.8)+
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label=labs), segment.color = "gray50",force=2, min.segment.length = 0, box.padding = 0.5 )+
  scale_x_continuous(trans="sqrt", limits = limt_axis)+
  scale_y_continuous(trans="sqrt", limits = limt_axis)+
  scale_fill_manual(values=colors.sgb.phylum)+
  scale_size_continuous(breaks = c(100,200,400,600,800,1000,1200,1400))+
  coord_fixed()+
  labs(x="Continent's Adjusted R² (controlled for Sampletype)", y="Sampletype's Adjusted R² (controlled for Continent)", title="ANI distance")+
  theme_bw()




