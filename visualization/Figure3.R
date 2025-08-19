source("OAPGC_function.R")

#heatmap
dt = read.table("./family.kegg.model_comp", sep="\t", header=T, check.names=F, row.names=1)
model.map = read.table("/share/data1/Database/KEGG/20230401/module.A_B_C.describe.tsv", sep="\t", header=F, comment.char = "", quote = "")
row.names(model.map) = model.map$V1
colnames(model.map) = c("model","A","B","C","desc")

info = read.table("./family.show.list", sep="\t", header=T, check.names = F, row.names=1)
info$family = row.names(info)



info = info[colnames(dt),]


order.phy = info %>%
  group_by(phylum) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(phylum)


order.family = info %>%
  arrange(desc(sgb.count), desc(abun.mean)) %>%
  pull("family")

count.family = length(order.family)


dtf = t(t(dt) / info$sgb.count)
dtf = dtf[rowSums(dtf) !=0, ] 

## core
x1 <- melt(dtf[rowMeans(dtf) == 1,]) %>%
  merge(model.map[,c("model","B")], by.y='model', by.x='Var1') %>%
  group_by(Var2,B) %>%
  summarise(value=sum(value)) %>%
  dcast(Var2 ~ B, value.var='value')
rownames(x1) = x1[,1]; x1 = x1[,-1]

## common
x2 <- dtf[rowSums(dtf>0)/count.family > 0.5 & rowMeans(dtf) !=1 ,] %>%
  melt() %>%
  merge(model.map[,c("model","B")], by.y='model', by.x='Var1') %>%
  group_by(Var2,B) %>%
  summarise(value=sum(value)) %>%
  dcast(Var2 ~ B, value.var='value')
rownames(x2) = x2[,1]; x2 = x2[,-1]

## rare
x3 <- melt(dtf[rowSums(dtf>0)/count.family <0.5,]) %>%
  merge(model.map[,c("model","B")], by.y='model', by.x='Var1') %>%
  group_by(Var2,B) %>%
  summarise(value=sum(value)) %>%
  dcast(Var2 ~ B, value.var='value')
rownames(x3) = x3[,1]; x3 = x3[,-1]

mo_status = data.frame(
  core = colSums(dtf[rowSums(dtf>0)/count.family == 1,]),# core
  common = colSums(dtf[rowSums(dtf>0)/count.family > 0.5,]), # common
  rare = colSums(dtf[rowSums(dtf>0)/count.family <= 0.5,]) # rare
)

info = merge(info, mo_status, by='row.names')
rownames(info) = info[,1]; info = info[,-1]
dtf = dtf[rowSums(dtf>0)/count.family > 0.5 & rowMeans(dtf)<1,]

heat_matrix = dtf[, order.family]
colors.taxo = c("#8dd3c7", "#08306b", "#08519c", "#4292c6", "#9ecae1", "#deebf7", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffed6f", "#b15928")
names(colors.taxo) = c("p__Actinomycetota", "p__Bacillota", "p__Bacillota_A", "p__Bacillota_B", "p__Bacillota_C", "p__Bacillota_I", "p__Bacteroidota", "p__Campylobacterota", "p__Chloroflexota", "p__Desulfobacterota", "p__Fusobacteriota", "p__Methanobacteriota", "p__Patescibacteria", "p__Pseudomonadota", "p__Spirochaetota", "p__Synergistota")

#----------------------------
#       col_anno
info = info[order.family, ] 
info$phylum = factor(info$phylum, levels = order.phy)
info$family = factor(info$family, levels = order.family)


color.col.model =c("#6a3c9a","#737373","#e2191c","white", "#253494", "#f781bf", "#377eb8", "#b30000", "#bebada", "#33a02c", "#d9d9d9", "#ff7f00")
names(color.col.model) = c("Glycan metabolism","Biosynthesis of other secondary metabolites","Xenobiotics biodegradation","Module set","Amino acid metabolism", "Biosynthesis of terpenoids and polyketides", "Carbohydrate metabolism", "Energy metabolism", "Gene set", "Lipid metabolism", "Metabolism of cofactors and vitamins", "Nucleotide metabolism")

x1 = x1[order.family,]
x2 = x2[order.family,]
x3 = x3[order.family,]

col_anno_bottom <-
  HeatmapAnnotation(which = "column" 
                    , core = anno_barplot(x1, gp = gpar(fill = color.col.model[colnames(x1)]), 
                                          bar_width = 1, height = unit(6, "cm"))
                    , common = anno_barplot(x2, gp = gpar(fill = color.col.model[colnames(x2)]), 
                                            bar_width = 1, height = unit(6, "cm"))
                    , rare = anno_barplot(x3, gp = gpar(fill = color.col.model[colnames(x3)]), 
                                          bar_width = 1, height = unit(6, "cm"))
                    # , sgb.count = anno_barplot(info$sgb.count, axis_param = list(at=c(0,50,100, 200)),  bar_width=1, height=unit(3,'cm'))
                    , abund.mean = anno_barplot(info$abun.mean,  bar_width=1, height=unit(3,'cm'))
                    # , rate.kegg = anno_barplot(info$kegg.rate.mean, bar_width=1, height=unit(3,'cm')) # 后面要删除掉
                    # , gene.count = anno_barplot(info$gene.count.mean, bar_width=1, height=unit(3,'cm'))
                    # , obs.kegg = anno_barplot(info$obs.kegg.mean, bar_width=1, height=unit(3,'cm'))
                    # , sgb.size = anno_barplot(info$len.mean, bar_width=1, height=unit(3,'cm'))
                    # , core = anno_barplot(info$core, bar_width=1, height=unit(3,'cm'))
                    # , common = anno_barplot(info$common, bar_width=1, height=unit(3,'cm'))
                    # , rare = anno_barplot(info$rare, bar_width=1, height=unit(3,'cm'))
                    , phylum = info$phylum 
                    , col = list(
                      phylum = colors.taxo
                      ,core = color.col.model
                    )
  )


color.kegg =c("#253494", "#f781bf", "#377eb8", "#b30000", "#bebada", "#33a02c", "#d9d9d9", "#ff7f00", "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#fccde5", "#f781bf", "#80b1d3", "#377eb8", "#d7301f", "#fc8d59", "#bebada", "#b3de69", "#d9d9d9", "#fdb462", "#ff7f00")
names(color.kegg) = c("Amino acid metabolism", "Biosynthesis of terpenoids and polyketides", "Carbohydrate metabolism", "Energy metabolism", "Gene set", "Lipid metabolism", "Metabolism of cofactors and vitamins", "Nucleotide metabolism", "Arginine and proline metabolism", "Lysine metabolism", "Cysteine and methionine metabolism", "Serine and threonine metabolism", "Branched-chain amino acid metabolism", "Aromatic amino acid metabolism", "Histidine metabolism", "Polyamine biosynthesis", "Terpenoid backbone biosynthesis", "Polyketide sugar unit biosynthesis", "Central carbohydrate metabolism", "Other carbohydrate metabolism", "Carbon fixation", "Methane metabolism", "Drug resistance", "Fatty acid metabolism", "Cofactor and vitamin metabolism", "Purine metabolism", "Pyrimidine metabolism")

model.map = model.map[rownames(heat_matrix),]
rownames(heat_matrix) = model.map$desc
order_level_B = table(model.map$B) %>% sort(decreasing = T) %>% names()
model.map$B = factor(model.map$B, levels=order_level_B)

row_anno_right <-
  HeatmapAnnotation(which = "row", 
                    levelB = model.map$B,
                    levelC = model.map$C,
                    col = list(
                      levelC = color.kegg,
                      levelB = color.kegg
                    )
  )



#------------------------------------
mycolor = colorRamp2(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     colors = c("#ffffff", "#ffffd9", "#c7e9b4","#41b6c4","#1d91c0","#225ea8"))

nh = nrow(heat_matrix)
nw = ncol(heat_matrix)
myunit = 4
max_nchar = max(nchar(rownames(heat_matrix)))



row.names(heat_matrix) = paste(model.map$model, ", ", sapply(strsplit(row.names(heat_matrix), ","), function(x) x[1]), sep="")
# row.names(heat_matrix) = paste(rownames(heat_matrix), ",", model.map$model)


ha.model <- Heatmap(heat_matrix
                    
                    # row
                    , show_row_dend = F
                    , show_row_names = T
                    , row_split = model.map$B 
                    , row_title_rot = 0 
                    , cluster_row_slices = F 
                    , row_gap = unit(0,'mm') 
                    
                    ## anno
                    ,right_annotation = row_anno_right,
                    
                    # column
                    ## anno
                    , top_annotation = col_anno_bottom # 
                    
                    ## heatmap
                    , column_split = info$phylum 
                    , cluster_column_slices = F 
                    , column_gap = unit(0,'mm') 
                    , cluster_columns = F 
                    , show_column_dend = F 
                    , show_column_names = T 
                    # , column_title = NULL,
                    , column_title_rot = 90 
                    
                    ## global
                    , border = T
                    ,col = mycolor
                    ,rect_gp = gpar(col = "grey", lwd = 0.1) 
                    ,height = nh * unit(myunit,"mm"), width = nw * unit(myunit,"mm"), 
                    ,row_names_max_width = unit(max_nchar, "char") 
                    
                    ,heatmap_legend_param = list(
                      at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                      title='rate',
                      legend_height = unit(8,'cm'),
                      ncol=5,nrow=2,
                      border=T
                    )
)
pdf("family.heatmap.function-model.pdf", width=32, height=32)
draw(ha.model)
dev.off()
