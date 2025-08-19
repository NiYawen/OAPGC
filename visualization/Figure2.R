source("OAPGC_function.R")

#sunburst
my <- read.table("sunplot.lineage.tsv", sep="\t", header=T, stringsAsFactors = F, check.names=F)

my = my[,-1]
my$plot_label = my$species

my <- my %>% 
  arrange(phylum, class, order, family, genus, species)


fake_circle<-c()
for (i in 1:nrow(my)){
  fake_circle<-append(fake_circle,rep(my$species[i],round(my$value[i]))) } # 这边需要改

edges<- data.frame(rbind(
  cbind(rep('origin',length(unique(my$phylum))),unique(as.character(my$phylum))),
  as.matrix(my[!duplicated(my[c('phylum','class')]),1:2]),
  as.matrix(my[!duplicated(my[c('class','order')]),2:3]),
  as.matrix(my[!duplicated(my[c('order','family')]),3:4]),
  as.matrix(my[!duplicated(my[c('family','genus')]),4:5]),
  as.matrix(my[!duplicated(my[c('genus','species')]),5:6])
)
)

colnames(edges)<-c('from','to')
vertices0<-data.frame(name=unique(c(as.character(edges$from), as.character(edges$to))))
my_leaf <- my[,c("species", 'value','plot_label')] 


vertices<-left_join(vertices0,my_leaf,by=c('name'='species'))

my_color<- data.frame(rbind(
  as.matrix(my[!duplicated(my[c('phylum','phylum')]),c(1,1)]),
  as.matrix(my[!duplicated(my[c('phylum','class')]),c(1,2)]),
  as.matrix(my[!duplicated(my[c('phylum','order')]),c(1,3)]),
  as.matrix(my[!duplicated(my[c('phylum','family')]),c(1,4)]),
  as.matrix(my[!duplicated(my[c('phylum','genus')]),c(1,5)]),
  as.matrix(my[!duplicated(my[c('phylum','species')]),c(1,6)])
  #as.matrix(my[!duplicated(my[c('phylum','strain_name')]),c(1,7)])
))

colnames(my_color)<-c('phylum','name')
sgb_info = read.table("./sunplot.group.tsv", sep="\t", header=T, check.names=F)
my_color = merge(sgb_info, my_color, by='name', all=T)

vertices<-left_join(vertices,my_color,by='name')

dm = merge(edges,vertices, by.x='to',by.y='name', all.x=T)
colnames(dm)[1:2] = c("node","parent")
rownames(dm) = paste(dm$parent,dm$node,sep="")
mm = paste(edges$from, edges$to, sep="")
dm = dm[mm,]
# dm = dm %>% filter(is.na(plot_label))


breaks <- c(1, 2, 6, 21, 101, 1001, Inf)
labels <- c("1", "2-5", "6-20", "21-100", "101-1000", "1000+")

dm$MAGs_cut <- cut(dm$MAGs_strains, 
                   breaks = breaks,
                   labels = labels,
                   include.lowest = TRUE,  
                   right = FALSE)        
dm$Isolated = as.numeric(dm$cultured_strains)
dm$Isolated_cut <- cut(dm$Isolated, 
                       breaks = breaks,
                       labels = labels,
                       include.lowest = TRUE,  # 包含最小值端点
                       right = FALSE)  


# write.table(dm, "temp.csv", row.names = F, sep = ",")



library(ggsunburst)
use_python("/share/data1/software/miniconda3/bin/python")
py_config()
sb <- sunburst_data("temp.csv", sep=",", 
                    type="node_parent",
                    node_attributes = c("plot_label","phylum","plot_label","group","MAGs_cut","Isolated_cut")
)

#sb$rects[!sb$rects$leaf,]$phylum <- sb$rects[!sb$rects$leaf,]$phylum


# 1、get the num of every taxo levels
list_names = c()
hide_names = c()

for (i in 1:(ncol(my)-3)){
  temp = aggregate(value ~ ., my[,c(i,7)],  sum) 
  
  tmpf <- subset(temp, value < 10)[,1]
  hide_names = c(hide_names, tmpf)
  
  x = temp[which(temp$value < 100), 1]
  list_names = c(list_names, x)
}

sb$node_labels[which(sb$node_labels$label %in% list_names),'pangle'] = sb$node_labels[which(sb$node_labels$label %in% list_names),'rangle']
sb$node_labels[which(sb$node_labels$label %in% list_names),'pvjust'] = sb$node_labels[which(sb$node_labels$label %in% list_names),'rhjust']
sb$node_labels[which(sb$node_labels$label %in% hide_names),'label'] = NA
sb$node_labels$label = gsub("^.__", "", sb$node_labels$label, perl=T)

# sb$node_labels$pangle[3] = sb$node_labels$rangle[3] 
####################
#   phylum
col_dt = read.table("./color.phylum.tsv", sep="\t", comment.char = "")
mycol = col_dt$V2
names(mycol) = col_dt$V1
p1 <- sunburst(sb, 
               rects.fill.aes = "phylum"
               ,rects.size = 0.5
               ,rects.color = "white"
               ,node_labels = T
               ,leaf_labels = F   
               ,leaf_labels.size=1
               ,node_labels.min = 0
               ,node_labels.size= 1.5 
)+
  scale_fill_manual(values=mycol)
p1


#######################
#       culture
colors.source = structure(c("#97afd7","#d98487","#f8c06b"),
                          names=c("Cultured(other)","Uncultured","Cultured(human)"))

p2 <- sunburst(sb, 
               rects.fill.aes = "group"
               ,rects.size = 0
               ,rects.color = "white"
               ,node_labels = T
               ,leaf_labels = F  
               #,leaf_labels.size=1
               ,node_labels.min = 0
               ,node_labels.size=0
)+
  scale_fill_manual(values=colors.source)
p2

#######################
#       MAGs-count
colors.counts = c("#f0f9e8", "#ccebc5", "#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac")
names(colors.counts) = c("1","2-5","6-20","21-100","101-1000","1000+")


p3 <- sunburst(sb, 
               rects.fill.aes = "MAGs_cut"
               ,rects.size = 0
               ,rects.color = "white"
               ,node_labels = T
               ,leaf_labels = F   # 最外层的叶子节点
               #,leaf_labels.size=1
               ,node_labels.min = 0
               ,node_labels.size=0
)+
  scale_fill_manual(values=colors.counts)
p3


#######################
#       Isolated-count
p4 <- sunburst(sb, 
               rects.fill.aes = "Isolated_cut"
               ,rects.size = 0
               ,rects.color = "white"
               ,node_labels = T
               ,leaf_labels = F   
               #,leaf_labels.size=1
               ,node_labels.min = 0
               ,node_labels.size=0
)+
  scale_fill_manual(values=colors.counts)
p4


p <- ggpubr::ggarrange(plotlist=list(p1,p2,p3,p4), ncol=2, nrow=2)
p

#bar
lab.ord <- dm %>%
  group_by(phylum) %>%
  count() %>%
  arrange(desc(n)) %>%
  ungroup() %>%
  mutate(label = paste(phylum, "(n = ",n,")", sep=""),
         label = gsub("p__", "", label))

# save(lab.ord, file="phylum.ord.RData")
pdt <- dm %>%
  group_by(phylum, group) %>%
  count()

pdt <- merge(pdt, lab.ord, by='phylum') %>%
  mutate(label = factor(label, levels=rev(lab.ord$label)),
         rate = n.x/n.y * 100)

p.phy <- ggplot(pdt, aes(y=label, x=rate, fill=group))+
  geom_bar(stat="identity", color="black", width=0.8)+
  scale_fill_manual(values=colors.source)+
  theme_bw()+
  theme(panel.grid = element_blank())

pdt %>%
  filter(group=="Uncultured") %>%
  arrange(rate)



