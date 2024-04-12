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

phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))


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

plot_data = df[, c("Temperature", "Max remove nodes")]
colnames(plot_data)[2] = "value"

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
davies.test(o, seg.Z = ~Temperature)



davies.test(o, seg.Z = ~Temperature) # 25.622, p-value = 1.407e-06

plot_data_1 = plot_data

p1 = ggplot(plot_data_1,
           aes(x = Temperature,
               y = value)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Temperature [°C]", y = "Max remove nodes") +
  annotate("text",
           label = "T = 25.622°C, p-value = 1.407e-06",
           x = 26.5, y =18,
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

head(phy_ann_data$OTU_id)

phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")




rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

# phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
#                               list(Class = phy_ann_abu$Class,
#                                    Genus = phy_ann_abu$Genus),
#                               sum)



phy_ann_genus_abu = phy_ann_abu

head(phy_ann_genus_abu)[1:4]
length(phy_ann_genus_abu$Class)

# phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
# phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


# rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
#                                     phy_ann_genus_abu$Genus,
#                                     sep = " ")

phy_ann_genus_abu$Class = NULL
phy_ann_genus_abu$Genus = NULL

LF_phy_ann_genus_abu = phy_ann_genus_abu[grepl(pattern = "3um", 
                                               colnames(phy_ann_genus_abu))]

SF_phy_ann_genus_abu = phy_ann_genus_abu[grepl(pattern = "0.2um", 
                                               colnames(phy_ann_genus_abu))]


LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))
SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))


LF_phy_ann_genus_abu = LF_phy_ann_genus_abu/rowSums(LF_phy_ann_genus_abu)
SF_phy_ann_genus_abu = SF_phy_ann_genus_abu/rowSums(SF_phy_ann_genus_abu)

# 大粒径的
LF_phy_ann_genus_abu$Stations = gsub("_", "-", rownames(LF_phy_ann_genus_abu))
LF_phy_ann_genus_abu$Stations = gsub("-3um-[1,2,3]", "", LF_phy_ann_genus_abu$Stations)
head(LF_phy_ann_genus_abu$Stations[1:6])

data = aggregate(LF_phy_ann_genus_abu[1:(length(LF_phy_ann_genus_abu)-1)],
                 by = list(Stations = LF_phy_ann_genus_abu$Stations),
                 mean)
rownames(data) = data$Stations
data$Stations = NULL
data = data[, colSums(data != 0) > 0]
head(data[1:4])
length(data)
# algae PCA 轴
# com1 = prcomp(data, center = TRUE,scale. = TRUE)
# df1 = com1$x
# summ = summary(com1)
# summ$importance[2,1]*100 # 17.471

# # algae NMDS1 轴
# 
# genus.distance = vegdist(data, method = 'bray')
# 
# #NMDS排序分析——vegan包中的metaMDS函数
# df_nmds = metaMDS(genus.distance, k = 2)
# #结果查看——关注stress、points及species三个指标
# summary(df_nmds)
# #应力函数值（<=0.2合理）
# df_nmds_stress = df_nmds$stress
# df_nmds_stress
# #检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
# stressplot(df_nmds)
# #提取作图数据
# df_points <- as.data.frame(df_nmds$points)
# #添加samp1es变量
# df_points$Stations <- row.names(df_points)
# #修改列名
# names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
# head(df_points)
# df1 = df_points

otu.distance = vegdist(data, method = 'bray')
pcoa = cmdscale(otu.distance,eig=TRUE)
pc12 = as.data.frame(pcoa$points[,1:2])
head(pc12)

pc = round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

df1 = pc12 # 37.08



# 环境因子 PCA 轴
head(daya_env)
env = daya_env[, c("SiO32-", "PO43-", "Temperature",
                   "NO3-", "NO2-", "NH4+",
                   "Stations")]



rownames(env) = env$Stations
env$Stations = NULL
env = scale(env)

com2 = prcomp(env, center = TRUE,scale. = TRUE)
df2 = com2$x
summ2 = summary(com2)
summ2$importance[2,1]*100 # 55.43


# 画图
library(ggrepel)
df1 = as.data.frame(df1)
df1$Stations = rownames(df1)
head(df1)


df2 = as.data.frame(df2)
df2$Stations = rownames(df2)

plot_data = merge(df2[, c("Stations", "PC1")],
                  df1[, c("Stations", "V1")],
                  by = "Stations")
plot_data = merge(plot_data,
                  daya_env[, c("Stations", 
                               "Temperature",
                               "PO43-",
                               "TChla")],
                  by = "Stations")

plot_data = plot_data[order(plot_data$Temperature), ]

plot_data_2 = plot_data


# gridding
library(gstat)
library(sp)

plot_data_2$Group = ifelse(plot_data_2$absT > 0,
                           "T < 25.622°C",
                           "T > 25.622°C")


p2 = ggplot(plot_data_2,  aes(x = Temperature, y = V1, 
                              col = Temperature, fill = Group)) +
  scale_fill_manual(values=c("#1597A5","#FEB3AE"))+
  scale_color_gradient(high = "red", low = "blue")+
  geom_point(size = 5) +
  geom_text_repel(aes(label = round(Temperature, digits = 2))) +
  stat_ellipse(geom = "polygon",level=0.9,
               linetype = 2, size=0.5, col = "black",
               alpha=0.1)+
  geom_path(data = subset(plot_data_2, Temperature > 25 & Temperature < 26)) +
#  geom_path() +
  geom_smooth(data = subset(plot_data_2, Temperature <= 25.622),
              aes(x = Temperature, y = V1), se = F, col = "black",
              linewidth = 2,
              method = "lm") +
  geom_smooth(data = subset(plot_data_2, Temperature >= 25.622),
              aes(x = Temperature, y = V1), se = F, col = "black",
                  linewidth = 2,
                  method = "lm") +
  theme_bw() +
  labs(x= "Temperature [°C]", y = "Algae Community PCoA1 (37.08%)") +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))
p2

library(ggpubr)
p3 = ggplot(plot_data_2,aes(x=Group, y=NMDS1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=Group),
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(),
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#1597A5","#FEB3AE"))+
  geom_signif(comparisons = list(c("T > 25.622°C","T < 25.622°C")),
              map_signif_level = T,
              test = t.test,
              y_position = c(1),
              tip_length = c(c(0,0)),
              size=2, textsize = 8,
              color="black")
p3
# 绘制y轴为PC1值的分组箱线图
p4 = ggplot(plot_data_2,aes(y=Group, x= Temperature))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=Group),
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(),
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank())+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#1597A5","#FEB3AE"))+
  geom_signif(comparisons = list(c("T > 25.622°C","T < 25.622°C")),
              map_signif_level = T,
              test = t.test,
              y_position = c(30),
              tip_length = c(c(0,0)),
              size=2, textsize = 8,
              color="black")
p4

# ggpubr::ggarrange()函数对图进行拼接
library(lemon)

p5 = p2 + p3 + plot_layout(widths = c(4, 1))
p6 = p4 + p5 + plot_layout(heights = c(1, 4))
p6













































































