source("OAPGC_function.R")

#scatter plot
dt = read.table("./strain.group", sep="\t", header=T)

dt <- dt %>%
  mutate(adj.size = (GenomeSize*(1-Contamination/100))/(GenomeSize * Completeness/100 ) )
colors.strain.source = c("HQMAG"="#fdb462", "Isolated"="#8dd3c7", "shared"="#80b1d3")
dt$group = factor(dt$group, levels=rev(c("HQMAG","Isolated","shared")))


p.size.qs <- ggplot(dt, aes(x=GenomeSize, y=QS, color=group, fill=group))+
  scale_fill_manual(values=colors.strain.source)+
  scale_color_manual(values=colors.strain.source)+
  geom_point(alpha=0.6, size=.5, shape=20)+
  labs(x="GenomeSize", y="QS")+
  theme_bw()+
  theme(panel.grid =  element_blank())

ggsave("genome-qs.pdf",p.size.qs, width=8, height=6.5)

x <- ggplot(dt, aes(x=GenomeSize, fill=group))+
  scale_fill_manual(values=colors.strain.source)+
  geom_histogram(color="black", bins=100)+
  theme_bw()

y <- ggplot(dt, aes(y=QS, fill=group))+
  scale_fill_manual(values=colors.strain.source)+
  geom_histogram(color="black", bins=100)+
  theme_bw()+
  scale_y_continuous(limits=c(60,100))

p2 <- ggpubr::ggarrange(plotlist=list(x,y))
ggsave("genome-qs.margin.pdf", p2, width=12, height=12)

#pie plot
pdt <- dt %>%
  mutate(quality = ifelse(rRNA==3 & tRNA>=18, "near-complete","High"),
         fill=paste(group, " (",quality,")", sep="")) %>%
  group_by(fill) %>%
  count()

p <- zy_pie(pdt, value="n", fill="fill")

