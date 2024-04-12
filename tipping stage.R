library(stringr)
library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(lemon)
library(igraph)
library(ggrepel)
library(ggprism)
library(ggsignif)
library(ggpubr)

# max romve nodes 的转折点
library(NetSwan)
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S/edge_inf")
edge_files = list.files()

daya_ann = read.csv("../OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)
colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")

phy_ann_data = subset(daya_ann, Class %in% algae_class)


dino_genus = subset(phy_ann_data, Phylum == "Dinoflagellata")
dino_genus = unique(paste(dino_genus$Class, dino_genus$Genus))
diat_genus = subset(phy_ann_data, Class == "Bacillariophyta")
diat_genus = unique(paste(diat_genus$Class, diat_genus$Genus))
chlo_genus = subset(phy_ann_data, Phylum == "Chlorophyta")
chlo_genus = unique(paste(chlo_genus$Class, chlo_genus$Genus))

dino_genus = as.data.frame(dino_genus)
diat_genus = as.data.frame(diat_genus)
chlo_genus = as.data.frame(chlo_genus)

len = length(edge_files)
df = as.data.frame(matrix(nrow = len, ncol = 12))

colnames(df) = c("Temperature", "Edge number", "Positive edge number", 
                 "Negative edge number", "Max remove nodes", 
                 "Node number", "Weight mean",
                 "Core node number", "Mean of edge number of nodes",
                 "Core node number in Dinoflagellata",
                 "Core node number in Bacillariophyta",
                 "Core node number in Chlorophyta")


k = 0

edge_files = edge_files[grep(".csv", edge_files)]

for (i in edge_files) {
  k = k + 1
  Temperature = gsub(".csv", "", i)
  temp_file = read.csv(i)
  edge_num = length(temp_file$color)
  pos_num = length(subset(temp_file, color == "positive")$color)
  neg_num = length(subset(temp_file, color == "negative")$color)
  nodes = unique(temp_file$Source)
  nodes = append(nodes, unique(temp_file$Target))
  nodes = unique(nodes)
  edge_num_set = c()
  core_node_num = 0
  dino_cores = 0
  diat_cores = 0
  chlo_cores = 0
  gra = graph.edgelist(as.matrix(temp_file[,c("Source", "Target")]), 
                       directed = FALSE)
  f4<-swan_combinatory(gra,10)
  f4 = as.data.frame(f4)
  colnames(f4) = c("remove", "betweeness", "degree", "cascading", "random")
  f5 = subset(f4, random<0.3)
  max_remove = max(f5$remove)*length(nodes)
  for (j in nodes) {
    nodes_edge_num = length(subset(temp_file, Source == j)$color)
    nodes_edge_num = nodes_edge_num + length(subset(temp_file, Target == j)$color)
    edge_num_set = append(edge_num_set, nodes_edge_num)
    if (nodes_edge_num > 10) {
      core_node_num = core_node_num + 1
      genu = gsub("_LF", "", j)
      genu = gsub("_SF", "", genu)
      if (genu %in% dino_genus$dino_genus) dino_cores = dino_cores + 1
      if (genu %in% diat_genus$diat_genus) diat_cores = diat_cores + 1
      if (genu %in% chlo_genus$chlo_genus) chlo_cores = chlo_cores + 1
    }
  }
  nodes_edge_num_mean = mean(edge_num_set)
  df$Temperature[k] = Temperature
  df$`Edge number`[k] = edge_num
  df$`Positive edge number`[k] = pos_num
  df$`Negative edge number`[k] = neg_num
  df$`Max remove nodes`[k] = max_remove
  df$`Node number`[k] = length(nodes)
  df$`Weight mean`[k] = mean(temp_file$weight)
  df$`Core node number`[k] = core_node_num
  df$`Mean of edge number of nodes`[k] = nodes_edge_num_mean
  df$`Core node number in Dinoflagellata`[k] = dino_cores
  df$`Core node number in Bacillariophyta`[k] = diat_cores
  df$`Core node number in Chlorophyta`[k] = chlo_cores
}

head(df)
library(segmented)



plot_data = df[, c("Temperature", "Node number")]



colnames(plot_data)[2] = "value"
plot_data$value[plot_data$value<0] = 0
plot_data$Temperature = as.numeric(plot_data$Temperature)

o = lm(value~Temperature, data=plot_data)
os = segmented(o, seg.Z=~Temperature, npsi=2,var.type="bootstrap", 
               ci.bootstrap.size=100)
summary(os)
davies.test(o,seg.Z = ~Temperature)
AIC(os)
ob = vcov(os, B=50, ret.b=TRUE)
confint(os, obj.boot=ob)


library(chngpt)
fit=chngptm(formula.1= value~1, formula.2=~Temperature, 
            plot_data, 
            type="stegmented", family="gaussian",
            var.type="bootstrap", ci.bootstrap.size=100)
print(summary(fit))
summary(os)
davies.test(o, seg.Z = ~Temperature) # 25.622, p-value = 5.785e-06



plot_data_1 = plot_data

plot_data_1$absT = plot_data_1$Temperature - 25.622
plot_data_1$Group = ifelse(plot_data_1$absT > 0,
                           "Region 3",
                           "Region 1")
plot_data_1$absT = abs(plot_data_1$absT)
plot_data_1$Group = ifelse(plot_data_1$absT < 0.5, "Region 2", plot_data_1$Group)

plot_data_1 = subset(plot_data_1, Group %in% c("Region 1", "Region 2",
                                               "Region 3"))

p1 = ggplot(plot_data_1,
            aes(x = Temperature,
                y = value)) +
  scale_color_brewer(palette = "Set2", name = "Region") +
  labs(x = "Temperature [°C]", y = "Nodes number") +
  geom_rect(aes(xmin = 23.5, xmax = 25.622-0.5, ymin = 35, ymax = 75), 
            fill = "#C2E7DB", alpha = 0.2) +
  geom_rect(aes(xmin = 25.622-0.5, xmax = 25.622+0.5, ymin = 35, ymax = 75), 
            fill = "#FED1C0", alpha = 0.2) +
  geom_rect(aes(xmin = 25.622+0.5, xmax = 29, ymin = 35, ymax = 75), 
            fill = "#D1D9EA", alpha = 0.2) +
  xlim(23.5, 29) +
  annotate("text",
           label = "Tipping point = 25.622°C, p-value = 5.785e-06",
           x = 26, y = 70,
           size = 6,
           col = "red") +
  geom_smooth(data = plot_data_1, aes(Temperature, value),
              method = "gam", formula = y ~ x + I((x - 25.622) * (x > 25.622)),
              se = T, color='blue', fill = "blue", alpha = 0.1)+
  theme_bw() +
  geom_point(size = 5) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
p1


daya_abu = read.csv("../OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

colnames(daya_abu)[1] = "OTU_id"
phy_ann_data$OTU_id = rownames(phy_ann_data)

head(daya_abu[1:7])

phy_ann_data = phy_ann_data[c("OTU_id", "Class", "Genus")]

head(phy_ann_data,12)

phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL


phy_ann_abu$Class = NULL
phy_ann_abu$Genus = NULL

head(phy_ann_abu)[1:5]
phy_ann_abu = as.data.frame(t(phy_ann_abu))

phy_ann_abu$Fraction = ifelse(grepl("0.2um",rownames(phy_ann_abu)),
                              "0.2-3 μm",
                              "3-200 μm")
head(phy_ann_abu[1:6])

phy_ann_abu$Stations = gsub("_", "-", rownames(phy_ann_abu))
phy_ann_abu$Stations = gsub("-3um-[1,2,3]", "", phy_ann_abu$Stations)
phy_ann_abu$Stations = gsub("-0.2um-[1,2,3]", "", phy_ann_abu$Stations)




tmp = phy_ann_abu[,1:(length(phy_ann_abu)-2)]

tmp = as.data.frame(t(apply(tmp, 1, function(x) x/sum(x))))
rowSums(tmp)
phy_ann_abu[,1:(length(phy_ann_abu)-2)] = tmp

phy_ann_abu = aggregate(phy_ann_abu[,1:(length(phy_ann_abu)-2)],
                        by = list(Fraction = phy_ann_abu$Fraction,
                                  Stations = phy_ann_abu$Stations),
                        mean)



# algae PCA 轴
phy_ann_abu_L = subset(phy_ann_abu, Fraction == "3-200 μm")
phy_ann_abu_S = subset(phy_ann_abu, Fraction == "0.2-3 μm")

data = phy_ann_abu_L[,3:length(phy_ann_abu)]


otu.distance = vegdist(data, method = 'bray')
pcoa = cmdscale(otu.distance,eig=TRUE)
pc12 = as.data.frame(pcoa$points[,1:2])
head(pc12)

pc = round(pcoa$eig/sum(pcoa$eig)*100,digits=2) # 37.08

df1 = pc12 


df1 = cbind(phy_ann_abu_L[c("Stations", "Fraction")], df1)


plot_data = merge(df1,
                  daya_env[, c("Stations", "Temperature")],
                  by = "Stations")

plot_data = plot_data[order(plot_data$Temperature), ]

plot_data_2 = plot_data
plot_data_2$absT = plot_data_2$Temperature - 25.622
plot_data_2$Group = ifelse(plot_data_2$absT > 0,
                           "Region 3",
                           "Region 1")
plot_data_2$absT = abs(plot_data_2$absT)
plot_data_2$Group = ifelse(plot_data_2$absT < 0.5, "Region 2", plot_data_2$Group)



cv = function(data) {
     cv = sd(data) / mean(data) * 100
      return(cv)
}


plot_data_2_1 = subset(plot_data_2, Group == "Region 1")
sd(plot_data_2_1$V1) # 0.12
length(plot_data_2_1$Stations) # 18

plot_data_2_2 = subset(plot_data_2, Group == "Region 2")
sd(plot_data_2_2$V1) # 0.20
length(plot_data_2_2$Stations) # 9


plot_data_2_3 = subset(plot_data_2, Group == "Region 3")
sd(plot_data_2_3$V1) # 0.22
length(plot_data_2_3$Stations) # 23



p2 = ggplot(plot_data_2,  aes(x = Temperature, y = V1, fill = Group)) +
  scale_fill_brewer(palette = "Set2", name = "Region")+
  scale_color_gradient(high = "red", low = "blue")+
  geom_point(size = 1) +
#  geom_text_repel(aes(label = round(Temperature, digits = 2))) +
  stat_ellipse(geom = "polygon",level=0.9,
               linetype = 2, size=0.5, col = "black",
               alpha=0.4)+
  geom_smooth(data = subset(plot_data_2, Group == "Region 1"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  geom_smooth(data = subset(plot_data_2, Group == "Region 2"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  geom_smooth(data = subset(plot_data_2, Group == "Region 3"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  annotate("text",
           label = "SD = 0.11, n = 18",
           x = 24, y =0,
           size = 4,
           col = "red") +
  annotate("text",
           label = "SD = 0.20, n = 9",
           x = 25.622, y =0.55,
           size = 4,
           col = "red") +
  annotate("text",
           label = "SD = 0.22, n = 23",
           x = 28, y = 0.55,
           size = 4,
           col = "red") +
  theme_bw() +
  labs(x= "Temperature [°C]", y = "Algae Community in 3-200 μm PCoA1 (38.46%)") +
  theme(legend.position = c(0.12, 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))
p2

library(ggpubr)
p3 = ggplot(plot_data_2,aes(x=Group, y=V1))+
  theme_bw() +
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=Group), alpha = 0.4,
               outlier.colour="white",size=0.5)+
  theme(axis.text= element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  ylab("Algae Community in 3-200 μm PCoA1 (38.46%)")+
  scale_fill_brewer(palette = "Set2", name = "Region")+
  geom_signif(comparisons = list(c("Region 1","Region 2"),
                                 c("Region 2","Region 3"),
                                 c("Region 1","Region 3")),
              map_signif_level = F,
              test = t.test,
              y_position = c(1, 0.8, 0.6),
              tip_length = c(c(0.02, 0.02, 0.02)),
              size=2, textsize = 5,
              color="black")
p3


p5 = p2 + p3 + plot_layout(widths = c(4, 1))
p5


data = phy_ann_abu_S[,3:length(phy_ann_abu)]


otu.distance = vegdist(data, method = 'bray')
pcoa = cmdscale(otu.distance,eig=TRUE)
pc12 = as.data.frame(pcoa$points[,1:2])
head(pc12)

pc = round(pcoa$eig/sum(pcoa$eig)*100,digits=2) # 24.40

df3 = pc12 

df3 = cbind(phy_ann_abu_S[c("Stations", "Fraction")], df3)


plot_data = merge(df3,
                  daya_env[, c("Stations", "Temperature")],
                  by = "Stations")




plot_data_3 = plot_data
plot_data_3$absT = plot_data_3$Temperature - 25.622
plot_data_3$Group = ifelse(plot_data_3$absT > 0,
                           "Region 3",
                           "Region 1")
plot_data_3$absT = abs(plot_data_3$absT)
plot_data_3$Group = ifelse(plot_data_3$absT < 0.5, "Region 2", plot_data_3$Group)

plot_data_3_1 = subset(plot_data_3, Group == "Region 1")
sd(plot_data_3_1$V1) # 0.14
length(plot_data_3_1$Stations) # 18

plot_data_3_2 = subset(plot_data_3, Group == "Region 2")
sd(plot_data_3_2$V1) # 0.22
length(plot_data_3_2$Stations) # 9


plot_data_3_3 = subset(plot_data_3, Group == "Region 3")
sd(plot_data_3_3$V1) # 0.24
length(plot_data_3_3$Stations) # 23



p6 = ggplot(plot_data_3,  aes(x = Temperature, y = V1, fill = Group)) +
  scale_fill_brewer(palette = "Set2", name = "Region")+
  scale_color_gradient(high = "red", low = "blue")+
  geom_point(size = 1) +
#  geom_text_repel(aes(label = round(Temperature, digits = 2))) +
  stat_ellipse(geom = "polygon",level=0.9,
               linetype = 2, size=0.5, col = "black",
               alpha=0.4)+
  geom_smooth(data = subset(plot_data_3, Group == "Region 1"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  geom_smooth(data = subset(plot_data_3, Group == "Region 2"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  geom_smooth(data = subset(plot_data_3, Group == "Region 3"),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  annotate("text",
           label = "SD = 0.14, n = 18",
           x = 24, y = 0.05,
           size = 4,
           col = "red") +
  annotate("text",
           label = "SD = 0.22, n = 9",
           x = 25.622, y =0.55,
           size = 4,
           col = "red") +
  annotate("text",
           label = "SD = 0.24, n = 23",
           x = 28, y = 0.55,
           size = 4,
           col = "red") +
  theme_bw() +
  labs(x= "Temperature [°C]", y = "Algae Community in 0.2-3 μm PCoA1 (24.40%)") +
  theme(legend.position = c(0.12, 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))
p6

library(ggpubr)
p7 = ggplot(plot_data_3,aes(x=Group, y=V1))+
  theme_bw() +
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=Group), alpha = 0.4,
               outlier.colour="white",size=0.5)+
  theme(axis.text= element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  ylab("Algae Community in 0.2-3 μm PCoA1 (24.40%)")+
  scale_fill_brewer(palette = "Set2", name = "Region")+
  geom_signif(comparisons = list(c("Region 1","Region 2"),
                                 c("Region 2","Region 3"),
                                 c("Region 1","Region 3")),
              map_signif_level = F,
              test = t.test,
              y_position = c(1, 0.8, 0.6),
              tip_length = c(c(0.02, 0.02, 0.02)),
              size=2, textsize = 5,
              color="black")
p7


p8 = p6 + p7 + plot_layout(widths = c(4, 1))
p8

p9 = p2/p6
p9

p10 = p3|p7
p10

p11 = (p1 / p10) | p9 
p11

ggsave("tipping point range.tiff", dpi = 300)



























