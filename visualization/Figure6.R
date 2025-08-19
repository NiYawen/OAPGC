source("OAPGC_function.R")

#prevalence in oral and gut
limits = c(0,100)
### 
p1 <- ggplot(pdt1f, aes(x=x.preva*100, y=y.preva*100))+
  # geom_point(aes(fill=.data[[color_tax]]), shape=21, size=3, color="black")+
  geom_point(aes(fill=shared_with_oral, size=y.mean, color='black', group=as.factor(shared_with_oral)), shape=21)+
  scale_color_manual(values=structure(c("black","red"), names=c("non-sig","sig")))+
  scale_fill_manual(values=colors.sgb.shared)+
  scale_size_continuous(trans = "sqrt", name = "Avg. gut", breaks = c(4,3,2,1,0.1,0.01,0.001, 0.0001, 0.00001))+
  scale_x_continuous(limits = limits)+
  scale_y_continuous(limits = limits)+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed")+
  geom_vline(xintercept = 20, color = "red", linetype = "dashed")+
  theme_bw()+
  coord_fixed()+
  labs(x = "Prevalence in oral", y = "Prevalence in gut")
p1

p2 <- ggplot(pdt1f, aes(x=x.preva, fill=shared_with_oral))+
  geom_histogram(bins=100, position = "stack", alpha=0.5)+
  # geom_histogram(bins=100)+
  # geom_density(alpha=0.5)+
  scale_fill_manual(values=colors.sgb.shared)+
  # facet_wrap(.~shared_with_oral, scales="free")+
  #scale_x_continuous(limits=c(0,1))+
  theme_bw()
p2

p3 <- ggplot(pdt1f, aes(y=y.preva, fill=shared_with_oral))+
  geom_histogram(bins=100, position = "stack", alpha=0.5)+
  # geom_density(alpha=0.5)+
  scale_fill_manual(values=colors.sgb.shared)+
  #facet_wrap(.~shared_with_oral, scales="free")+
  #scale_y_continuous(limits=c(0,1))+
  theme_bw()
p3


p <- ggpubr::ggarrange(
  plotlist=
    list(p2, NULL,
         p1, p3), 
  ncol=2, nrow=2, widths=c(3,1), heights=c(1,3), legend = 'right')
p



#fecal prevalence and MRA
plot_dt <- subset(pdt1f, preference == "oral.high")

sgb.taxo = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
             "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
             "#cab2d6", "#6a3d9a", "#b15928", "#8dd3c7",
             "blue","darkblue")

top.sgb = plot_dt %>%
  count(genus) %>%
  arrange(desc(n)) %>%
  filter(row_number()<=12 | genus %in% c("g__Bacteroides_D", "g__Bacteroides"))
names(sgb.taxo) = top.sgb$genus

lab1 <- plot_dt %>%
  filter(y.mean > 0.01) %>%
  pull(scientific_name)

lab2 <- plot_dt %>%
  slice_max(y.preva, n=10) %>%
  pull(scientific_name)
mylabs = c(lab1, lab2)

plot_dt <- plot_dt %>%
  mutate(labels = ifelse(scientific_name %in% mylabs, scientific_name, NA),
         color = ifelse(genus %in% names(sgb.taxo), genus, NA))



p1 <- ggplot(plot_dt, aes(x=y.preva, y = y.mean, fill=genus))+
  geom_point(aes(size=y.preva*y.preva),shape=21, color='black', alpha=0.9)+
  geom_text_repel(aes(label=labels), size=3, min.segment.length = 0
  )+
  #geom_point(shape=21, color='black', alpha=1, size=2)+
  scale_y_continuous(trans='log10')+
  geom_vline(xintercept = c(0.2), lty="dashed", color="red")+
  scale_size_continuous(name="mean * preva", breaks=c(0.4,0.2,0.3,0.1,0.01,0.001,0.0001,0.00001,0.000001) )+
  scale_fill_manual(values=sgb.taxo, na.value = "#d9d9d9")+
  scale_x_continuous(trans='sqrt', breaks = c(0.01,0.05,0.1,0.25,0.5,0.75,1), limits=c(0,1))+
  theme_bw()+
  labs(x="gut.prevalence", y="Avg. abundance(%)")+
  guides(
    fill = guide_legend(override.aes = list(size = 5))  
  )
p1

px = ggplot(plot_dt, aes(x=y.preva))+
  geom_histogram(fill="#ff9800",  color="black", bins=100)+
  scale_x_continuous(trans='sqrt', breaks = c(0.01,0.05,0.1,0.25,0.5,0.75,1), limits=c(0,1))+
  scale_color_manual(values=sgb.taxo, na.value = "#d9d9d9")+
  theme_bw()+
  theme(legend.position = 'none')

px1 = ggplot(plot_dt, aes(x=y.preva))+
  scale_x_continuous(trans='sqrt', breaks = c(0.01,0.05,0.1,0.25,0.5,0.75,1), limits=c(0,1))+
  scale_color_manual(values=sgb.taxo, na.value = "#d9d9d9")+
  theme_bw()+
  theme(legend.position = 'none')

py = ggplot(plot_dt, aes(y=y.mean, color=color))+
  geom_histogram(fill="#ff9800", color="black", bins=100)+
  scale_y_log10()+
  scale_color_manual(values=sgb.taxo, na.value = "#d9d9d9")+
  theme_bw()+
  theme(legend.position = 'none')

p <- ggpubr::ggarrange(plotlist=list(px, NULL,p1, py), ncol=2, nrow=2, widths = c(3,1), heights = c(1,3), legend = 'right', common.legend = T)

p

