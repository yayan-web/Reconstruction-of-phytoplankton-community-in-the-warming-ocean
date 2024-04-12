setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

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

daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann = read.csv("OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)


colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")


algae_class = subset(algae_class, algae_class != "")
algae_class = subset(algae_class, algae_class != "Embryophyceae")


phy_ann_data = subset(daya_ann, Class %in% algae_class)


unique(phy_ann_data[, c('Class')])



phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
colnames(daya_abu)[1] = "OTU_id"


phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")

head(phy_ann_abu)[1:5]

colSums(phy_ann_abu[4:295])


count_data = as.data.frame(t(phy_ann_abu[4:295]))
observed_species = estimateR(count_data)[1,]
chao1 = estimateR(count_data)[2,]
shannon <- vegan::diversity(count_data, index = 'shannon', base = exp(1))
simpson <- vegan::diversity(count_data, index = 'simpson', base = exp(1))
divs = cbind.data.frame(observed_species, chao1, shannon, simpson)

divs["Stations"] = gsub("_", "-", rownames(divs))
divs["Stations"] = gsub("-3um-[1,2,3]", "", divs$Stations)
divs["Stations"] = gsub("-0.2um-[1,2,3]", "", divs$Stations)

divs$Sample_id = rownames(divs)
divs["Fraction"] = ""

for (i in c(1:length(divs$Fraction))) {
  if (str_count(rownames(divs)[i], "3um") > 0) {
    divs$Fraction[i] = "3-200 μm"
  }
  if (str_count(rownames(divs)[i], "0.2um") > 0) {
    divs$Fraction[i] = "0.2-3 μm"
  }
}

head(divs)
divs$Sample_id = NULL

divs = aggregate(divs[1:4],
                 by = list(Fraction = divs$Fraction,
                           Stations = divs$Stations),
                 mean)


divs = merge(divs, daya_env[c("Stations", "Temperature")], by = "Stations")

divs_long = melt(divs, measure.vars = c("observed_species", "chao1", "shannon", "simpson"))

x11()

colnames(divs_long)
divs_long = subset(divs_long, variable %in% c("observed_species", "shannon"))

windowsFonts(Arial = windowsFont("Arial"))


divs_long$variable = gsub("observed_species", "Observed Species", divs_long$variable)
divs_long$variable = gsub("shannon", "Shannon index", divs_long$variable)
p1 = ggplot(divs_long, aes(x = Temperature,
                           y = value,
                           col = Fraction)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(Fraction~variable, scales = "free_y", 
             nrow = 1) +
  theme(legend.position = "None",
        strip.text = element_text(size = 12, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))


p1

divs_long_1 = subset(divs_long, Fraction == "0.2-3 μm"&variable=="Observed Species")
sd(divs_long_1$value) # 29.41
 
divs_long_2 = subset(divs_long, Fraction == "3-200 μm"&variable=="Observed Species")
sd(divs_long_2$value) # 58.73

divs_long_3 = subset(divs_long, Fraction == "0.2-3 μm"&variable=="Shannon index")
sd(divs_long_3$value) # 0.53

divs_long_4 = subset(divs_long, Fraction == "3-200 μm"&variable=="Shannon index")
sd(divs_long_4$value) # 0.69



p2 = ggplot(divs_long, aes(x = Fraction, y = value, fill = Fraction)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_text(data = data.frame(Fraction = c("0.2-3 μm",  "3-200 μm", "0.2-3 μm", "3-200 μm"),
                              value = c(450, 450, 4, 4),
                              variable = c("Observed Species", "Observed Species",
                                           "Shannon index", "Shannon index"),
                              label = c("SD = 29.31", "SD = 57.90",
                                        "SD = 0.53", "SD = 0.70")),
            aes(label = label),
            size = 5, col = "red") +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  theme(legend.position = "None",
        strip.text = element_text(size = 12, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

p2

p = p1 + p2 + plot_layout(widths = c(2, 1))
p

ggsave("daya_divs.tiff", dpi = 300)


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
                              list(Class = phy_ann_abu$Class,
                                   Genus = phy_ann_abu$Genus),
                              sum)

phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
                                    phy_ann_genus_abu$Genus,
                                    sep = " ")

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

LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))
SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))

LF_phy_ann_genus_abu$Sum = rowSums(LF_phy_ann_genus_abu)
SF_phy_ann_genus_abu$Sum = rowSums(SF_phy_ann_genus_abu)

attach(LF_phy_ann_genus_abu)
attach(SF_phy_ann_genus_abu)
LF_phy_ann_genus_abu = LF_phy_ann_genus_abu[order(LF_phy_ann_genus_abu$Sum,
                                                  decreasing = T),]
SF_phy_ann_genus_abu = SF_phy_ann_genus_abu[order(SF_phy_ann_genus_abu$Sum,
                                                  decreasing = T),]
detach(LF_phy_ann_genus_abu)
detach(SF_phy_ann_genus_abu)

LF_phy_ann_genus_abu$Sum = NULL
SF_phy_ann_genus_abu$Sum = NULL


LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))

colnames(LF_phy_ann_genus_abu)[1:13]


LF_phy_ann_genus_abu_top = LF_phy_ann_genus_abu[c("Dinophyceae Heterocapsa",
                                                  "Dinophyceae Dinophyceae_XXX",
                                                  "Dinophyceae Pyrophacus",
                                                  "Dinophyceae Tripos",
                                                  "Bacillariophyta Thalassiosira",
                                                  "Bacillariophyta Polar-centric-Mediophyceae_X")]


SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))

colnames(SF_phy_ann_genus_abu)[1:20]
SF_phy_ann_genus_abu_top = SF_phy_ann_genus_abu[c("Mamiellophyceae Ostreococcus",
                                                  "Mamiellophyceae Micromonas",
                                                  "Mamiellophyceae Bathycoccus",
                                                  "Trebouxiophyceae Picochlorum",
                                                  "Prymnesiophyceae Chrysochromulina",
                                                  "Eustigmatophyceae Eustigmatophyceae_XXX")]

head(SF_phy_ann_genus_abu_top)

LF_phy_ann_genus_abu_top$Stations = gsub("_", "-", 
                                         rownames(LF_phy_ann_genus_abu_top))
LF_phy_ann_genus_abu_top$Stations = gsub("-3um-[1,2,3]", "",
                                         LF_phy_ann_genus_abu_top$Stations)

SF_phy_ann_genus_abu_top$Stations = gsub("_", "-", 
                                         rownames(SF_phy_ann_genus_abu_top))
SF_phy_ann_genus_abu_top$Stations = gsub("-0.2um-[1,2,3]", "",
                                         SF_phy_ann_genus_abu_top$Stations)



LF_phy_ann_genus_abu_top = merge(daya_env[c("Stations", "Temperature")],
                                 LF_phy_ann_genus_abu_top,
                                 by = "Stations")

LF_phy_ann_genus_abu_top_long = melt(LF_phy_ann_genus_abu_top,
                                     id.vars = c("Stations", "Temperature"))

LF_phy_ann_genus_abu_top_long[c("Class", "Genus")] = str_split_fixed(LF_phy_ann_genus_abu_top_long$variable,
                                                                      pattern = " ",
                                                                      2)

SF_phy_ann_genus_abu_top = merge(daya_env[c("Stations", "Temperature")],
                                 SF_phy_ann_genus_abu_top,
                                 by = "Stations")

SF_phy_ann_genus_abu_top_long = melt(SF_phy_ann_genus_abu_top,
                                     id.vars = c("Stations", "Temperature"))

SF_phy_ann_genus_abu_top_long[c("Class", "Genus")] = str_split_fixed(SF_phy_ann_genus_abu_top_long$variable,
                                                                      pattern = " ",
                                                                      2)

x11()

LF_phy_ann_genus_abu_top_long = aggregate(LF_phy_ann_genus_abu_top_long$value,
                                          list(Stations = LF_phy_ann_genus_abu_top_long$Stations,
                                               Temperature = LF_phy_ann_genus_abu_top_long$Temperature,
                                               variable = LF_phy_ann_genus_abu_top_long$variable,
                                               Class = LF_phy_ann_genus_abu_top_long$Class,
                                               Genus = LF_phy_ann_genus_abu_top_long$Genus),
                                          mean)

colnames(LF_phy_ann_genus_abu_top_long)[6] = "value"

p = ggplot(LF_phy_ann_genus_abu_top_long,
           aes(x = Temperature,
               y = log2(value),
               col = Class)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Log2(Relative abundance in 3-200 μm)") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(~Genus, scales = "free_y") +
  theme(legend.position = "right",
        strip.text = element_text(size = 8, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 8))
p

SF_phy_ann_genus_abu_top_long = aggregate(SF_phy_ann_genus_abu_top_long$value,
                                          list(Stations = SF_phy_ann_genus_abu_top_long$Stations,
                                               Temperature = SF_phy_ann_genus_abu_top_long$Temperature,
                                               variable = SF_phy_ann_genus_abu_top_long$variable,
                                               Class = SF_phy_ann_genus_abu_top_long$Class,
                                               Genus = SF_phy_ann_genus_abu_top_long$Genus),
                                          mean)

colnames(SF_phy_ann_genus_abu_top_long)[6] = "value"


p1 = ggplot(SF_phy_ann_genus_abu_top_long,
            aes(x = Temperature,
                y = log2(value),
                col = Class)) +
  scale_color_brewer(palette = "Set1") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Log2(Relative abundance in 0.2-3 μm)") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(~Genus, scales = "free_y") +
  theme(legend.position = "right",
        strip.text = element_text(size = 8, family = "Arial", face = "bold"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8))
p1

p2 = p1 / p

p2

ggsave("daya_genus.tiff", dpi = 300)


# tara divs

tara_env = read.csv("tara_env.txt", sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

tara_phy = subset(data, chloroplast == "yes")
tara_phy = as.data.frame(tara_phy)

tara_phy = tara_phy[, c(4:(length(tara_phy)-9), length(tara_phy))]

tara_phy = subset(tara_phy, id %in% ann_data$id)


tara_phy_abu = tara_phy[1:(length(tara_phy)-1)]

tara_phy_abu = as.data.frame(t(tara_phy_abu))
observed_species = estimateR(tara_phy_abu)[1,]
chao1 = estimateR(tara_phy_abu)[2,]

tara_phy_abu[1:5, 1:5]

shannon = vegan::diversity(tara_phy_abu, index = 'shannon', base = exp(1))
simpson = vegan::diversity(tara_phy_abu, index = 'simpson', base = exp(1))
divs = cbind.data.frame(observed_species, chao1, shannon, simpson)
divs$Barcode = rownames(divs)

tara_env_divs = merge(tara_env[c("Barcode", "Temperature",
                                 "Depth_Nominal", "Fraction_Lower",
                                 "Fraction_Upper")],
                      divs,
                      by = "Barcode")


tara_env_divs$Fraction = paste(tara_env_divs$Fraction_Lower,
                               tara_env_divs$Fraction_Upper,
                               sep = "-")
tara_env_divs$Fraction = paste(tara_env_divs$Fraction,
                               "μm")
unique(tara_env_divs$Fraction)
tara_env_divs_2 = subset(tara_env_divs, Fraction %in% c("0.8-5.00 μm",
                                                        "5-20.00 μm",
                                                        "20-180.00 μm"))

tara_env_divs_2 = subset(tara_env_divs_2, Depth_Nominal == "SRF")

tara_divs_long = melt(tara_env_divs_2, measure.vars = c("observed_species", "chao1", 
                                                        "shannon", "simpson"))

tara_divs_long = subset(tara_divs_long, Temperature > 23&Temperature < 30)

# x11()
tara_divs_long$Fraction = factor(tara_divs_long$Fraction,
                                 levels = c("0.8-5.00 μm",
                                            "5-20.00 μm",
                                            "20-180.00 μm"))

p = ggplot(tara_divs_long, aes(x = Temperature,
                               y = value,
                               col = variable)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, linewidth = 2) +
  facet_wrap(Fraction~variable, scales = "free_y", nrow = 3) +
  theme(legend.position = "None",
        strip.text = element_text(size = 15, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))


p
ggsave("tara_divs.tiff", dpi = 300)

# 网络
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann = read.csv("OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)


colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")

phy_ann_data = subset(daya_ann, Class %in% algae_class)

phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")


colnames(daya_abu)[1] = "OTU_id"
phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
                              list(Class = phy_ann_abu$Class,
                                   Genus = phy_ann_abu$Genus),
                              sum)

phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
                                    phy_ann_genus_abu$Genus,
                                    sep = " ")

phy_ann_genus_abu$Class = NULL
phy_ann_genus_abu$Genus = NULL

phy_ann_genus_abu = as.data.frame(t(phy_ann_genus_abu))

phy_ann_genus_abu = phy_ann_genus_abu/rowSums(phy_ann_genus_abu)

phy_ann_genus_abu = phy_ann_genus_abu[!grepl(pattern = "Unknown", 
                                             colnames(phy_ann_genus_abu))]



phy_ann_genus_abu["Stations"] = gsub("_", "-", rownames(phy_ann_genus_abu))
phy_ann_genus_abu["Stations"] = gsub("-3um-[1,2,3]", "", phy_ann_genus_abu$Stations)
phy_ann_genus_abu["Stations"] = gsub("-0.2um-[1,2,3]", "", phy_ann_genus_abu$Stations)


phy_ann_genus_abu$Sample_id = rownames(phy_ann_genus_abu)
phy_ann_genus_abu["Fraction"] = ""

for (i in c(1:length(phy_ann_genus_abu$Fraction))) {
  if (str_count(rownames(phy_ann_genus_abu)[i], "3um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "3-200 μm"
  }
  if (str_count(rownames(phy_ann_genus_abu)[i], "0.2um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "0.2-3 μm"
  }
}

phy_ann_genus_env = merge(daya_env[c("Stations", "Temperature")], 
                          phy_ann_genus_abu,
                          by = "Stations")


phy_ann_genus_env = aggregate(phy_ann_genus_env[2:(length(phy_ann_genus_env)-2)],
                              by = list(Stations = phy_ann_genus_env$Stations,
                                        Fraction = phy_ann_genus_env$Fraction),
                              mean)

phy_ann_genus_env_sum = aggregate(phy_ann_genus_env[4:(length(phy_ann_genus_env)-2)],
                                  by = list(Fraction = phy_ann_genus_env$Fraction),
                                  sum)

rownames(phy_ann_genus_env_sum) = phy_ann_genus_env_sum$Fraction
phy_ann_genus_env_sum$Fraction = NULL
phy_ann_genus_env_sum = as.data.frame(t(phy_ann_genus_env_sum))
phy_ann_genus_env_sum$ratio = phy_ann_genus_env_sum$`3-200 μm`/phy_ann_genus_env_sum$`0.2-3 μm`

LF_genus = rownames(subset(phy_ann_genus_env_sum, ratio >= 1))
SF_genus = rownames(subset(phy_ann_genus_env_sum, ratio < 1))

LF_genus = append(x = "Stations", LF_genus)
SF_genus = append(x = "Stations", SF_genus)

LF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "3-200 μm")
SF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "0.2-3 μm")

LF_phy_ann_genus_env = LF_phy_ann_genus_env[LF_genus]
SF_phy_ann_genus_env = SF_phy_ann_genus_env[SF_genus]

phy_ann_genus_env = merge(LF_phy_ann_genus_env,
                          SF_phy_ann_genus_env,
                          by = "Stations")


phy_ann_genus_env = merge(daya_env[c("Stations", "Temperature")],
                          phy_ann_genus_env,
                          by = "Stations")



setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")
library(igraph)
library(NetSwan)
library(psych)
library(dplyr)



genus_env = phy_ann_genus_env
rownames(genus_env) = genus_env$Stations
genus_env$Stations = NULL

attach(genus_env)
genus_env = genus_env[order(genus_env$Temperature),]
detach(genus_env)



i = 0
len = length(genus_env$Temperature) - 3
df = as.data.frame(matrix(nrow = len-3, ncol = 4))
colnames(df) = c("Temperature", "Edge number")
df$Temperature = genus_env$Temperature[4:len]


for (i in c(4:len)) {
  j = i - 3
  genus_env_part = genus_env[c(i-3,i-2,i-1, i, i+1, i+2,i+3),]
  for (k in colnames(genus_env_part)) {
    n = nrow(genus_env_part[genus_env_part[k]>0,])
    if (n < 3) {
      genus_env_part[k] = NULL
    }
  }
  genus_env_part$Temperature = NULL
  genus_env_part = scale(genus_env_part)
  corr = corr.test(genus_env_part, method = "spearman", 
                   ci = FALSE, adjust = "fdr")
  relation = corr$r
  pvalue = corr$p
  relation[abs(relation)<0.6|pvalue>0.05] = 0
  dr = as.data.frame(relation)
  dr$species1 = rownames(dr)
  edge = melt(dr, id.var = "species1")
  colnames(edge)[2:3] = c("species2", "correction")
  edge$color = ""
  edge$color = ifelse(edge$correction>0, "positive", "negative")
  edge$correction = abs(edge$correction)
  edge = edge[edge$species1 != edge$species2,]
  edge = subset(edge, abs(correction) > 0)
  colnames(edge) = c("Source", "Target", "weight", "color")
  edge$unique = ""
  for (i in c(1:length(edge$Source))) {
    if (as.character(edge$Source[i]) >= as.character(edge$Target[i])) edge$unique[i] = paste(edge$Source[i], edge$Target[i])
    if (as.character(edge$Source[i]) < as.character(edge$Target[i])) edge$unique[i] = paste(edge$Target[i], edge$Source[i])
  }
  edge = edge %>% distinct(unique, .keep_all =TRUE)
  edge$unique = NULL
  file_name = paste("edge_inf/", df$Temperature[j], ".csv", sep = "")
  write.csv(edge, file_name, row.names=FALSE, quote=FALSE)
  edge_number = length(edge$Source)
  df$`Edge number`[j] = edge_number
}


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


dino_genus = subset(phy_ann_data, Class == "Dinophyceae")
dino_genus = unique(paste(dino_genus$Class, dino_genus$Genus))

pry_genus = subset(phy_ann_data, Class == "Prymnesiophyceae")
pry_genus = unique(paste(pry_genus$Class, pry_genus$Genus))


diat_genus = subset(phy_ann_data, Class == "Bacillariophyta")
diat_genus = unique(paste(diat_genus$Class, diat_genus$Genus))
chlo_genus = subset(phy_ann_data, Phylum == "Chlorophyta")
chlo_genus = unique(paste(chlo_genus$Class, chlo_genus$Genus))

dino_genus = as.data.frame(dino_genus)
pry_genus = as.data.frame(pry_genus)
diat_genus = as.data.frame(diat_genus)
chlo_genus = as.data.frame(chlo_genus)


len = length(edge_files)


df = as.data.frame(matrix(nrow = len, ncol = 13))

colnames(df) = c("Temperature", "Edge number", "Positive edge number", 
                 "Negative edge number", "Max remove nodes", 
                 "Node number", "Weight mean",
                 "Core node number", "Mean of edge number of nodes",
                 "Core node number in Dinophyceae",
                 "Core node number in Prymnesiophyceae",
                 "Core node number in Bacillariophyta",
                 "Core node number in Chlorophyta")


k = 0

edge_files = edge_files[grep(".csv", edge_files)]

library(igraph)
library(NetSwan)
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
  pry_cores = 0
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
      if (genu %in% pry_genus$pry_genus) pry_cores = pry_cores + 1
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
  df$`Core node number in Dinophyceae`[k] = dino_cores
  df$`Core node number in Prymnesiophyceae`[k] = pry_cores
  df$`Core node number in Bacillariophyta`[k] = diat_cores
  df$`Core node number in Chlorophyta`[k] = chlo_cores
}


head(df)

write.csv(df, "../net_stat.csv")


library(ggplot2)
library(reshape2)

env_data = read.csv("../env_data_final.csv", header = T)

df2 = merge(df, env_data[, c("Temperature", "TChla")], by = "Temperature")

colnames(df)


df_long = melt(df2, id.vars = c("Temperature", "TChla"))

df_long$Temperature = as.numeric(df_long$Temperature)


df_long$Region = ''

for (i in c(1:length(df_long$Temperature))) {
  if (df_long$Temperature[i] <= 26.6) df_long$Region[i] = "Low Temperature"
  if (df_long$Temperature[i] > 26.6) df_long$Region[i] = "High Temperature"
}


color = c("#03090C", "#03090C", "#03090C", "#03090C",
          "#03090C", "#03090C", "#03090C", "#03090C",
          "#8DA0CB", "#A6D854", "#66C2A5", "#FC8D62")



p = ggplot(df_long, aes(x = Temperature, y = value, col = variable)) +
  labs(x = "Temperature [°C]",
       y = "Value") +
  scale_color_manual(values = color, name = "Region") +
  theme_bw() +
  #  geom_vline(xintercept = 26.6, linewidth = 2) +
  geom_point(alpha = 0.5) +
  labs(y = "") +
  facet_wrap(~variable, scales = "free_y") +
  #  geom_smooth(method = "glm", linewidth = 2, se = F, formula = y ~ x+I(x^2)) +
  geom_smooth(linewidth = 2,) +
  theme(legend.position = "none",
        text = element_text(family = "serif", face = "bold"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 15)) 
# x11()
p

ggsave("../network_temperatrue_spearman.pdf", dpi = 300)

library(network)

edge_data = read.csv("28.083.csv")
colnames(edge_data)
edge_data$color = gsub("positive", "green", edge_data$color)
edge_data$color = gsub("negative", "red", edge_data$color)

node_data = rbind(edge_data$Source, edge_data$Target)
node_data = str_split_fixed(node_data, " ", 2)
node_data = as.data.frame(node_data)
node_data = unique(node_data)

colnames(node_data) = c("class", "genus")
node_data = node_data[c("genus", "class")]
node_data$genus = paste(node_data$class, node_data$genus)
head(node_data)

unique(node_data$class)

key_class = c("Dinophyceae", "Prymnesiophyceae", "Bacillariophyta", "Chlorophyta")

node_data$class = ifelse(!(node_data$class %in% key_class), "others", key_class)


unique(node_data$class)



# X11()
library(ggplot2)
library(ggnet)

net_1 = network(edge_data, node_data)
p1 = ggnet2(net_1, size = 5,
            color = "class", 
            alpha = 0.8,
            edge.alpha = 0.5,
            palette = "Set2",
            mode = "circrand",
            edge.color = "color") +
  
  theme(legend.position = "none")
p1

net_2 = network(edge_data, node_data)
p2 = ggnet2(net_2, size = 5,
            color = "class", 
            alpha = 0.8,
            edge.alpha = 0.5,
            palette = "Set2",
            mode = "circrand",
            edge.color = "color") +
  theme(legend.position = "none")
p2

net_3 = network(edge_data, node_data)
p3 = ggnet2(net_3, size = 5,
            color = "class", 
            alpha = 0.8,
            edge.alpha = 0.5,
            palette = "Set2",
            mode = "circrand",
            edge.color = "color")
p3

library(patchwork)
p4 = (p1 + p2 + p3)/p +
  plot_layout(widths = c(5, 4, 6), heights = c(3,7))

p4


ggsave("../network_temperatrue_spearman_2.tiff", dpi = 300)


## CCA

# LASSO
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

library(glmnet) # LASSO

daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_pig = read.csv("pigment.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")


env_data = daya_env[c("Stations", "Temperature", "NO3-", "NO2-", 
                      "PO43-", "SiO32-", "NH4+")]

pig_data = daya_pig
pig_data$Longitude = NULL
pig_data$Latitude = NULL


rownames(env_data) = env_data$Stations
env_data$Stations = NULL

length(pig_data)
pig_data = aggregate(pig_data[2:20],
                     by = list(Stations = pig_data$Stations),
                     mean)

rownames(pig_data) = pig_data$Stations
pig_data$Stations = NULL

setdiff(rownames(env_data), rownames(pig_data))
setdiff(rownames(pig_data), rownames(env_data))

env_data = subset(env_data, rownames(env_data) != "M27-1")
env_data$Tabs = abs(env_data$Temperature - 24.8805)


pig_env = merge(env_data,
                pig_data,
                by = "row.names")

pig_env$Row.names = NULL
pig_env_long = melt(pig_env, id.vars = "Temperature")
p = ggplot(data = pig_env_long,
           aes(x = Temperature,
               y = value)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~variable, scale = "free_y") +
  geom_smooth(method = "lm", formula = y~x+I(x^2))

p

colnames(pig_env)

data = pig_env[1:8]
model = lm(`Chlorophyll c3`~`Tabs`+ `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
           data = data)
a = summary(model)

a$fstatistic[1]

model = glm(pigment ~`Temperature`+ `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
            data = data)
a = summary(model)


model = lm(`Chlorophyll c3`~`Temperature`+ `NO3-` + `NO2-` + `SiO32-` + `NH4+`,
           data = data)
a = summary(model)
a$aic

model = glm(`Chlorophyll c3`~`Tabs`+ `NO3-` + `NO2-` + `SiO32-` + `NH4+`,
            data = data)
a = summary(model)
a$fstatistic[1]


df = as.data.frame(matrix(nrow = 19, ncol = 3))
row = 0
for (i in colnames(pig_data)) {
  row = row + 1
  data = merge(env_data, pig_data[c(i)], by = "row.names")
  data$Row.names = NULL
  colnames(data)[8] = "pigment"
  model = lm(pigment ~  Temperature + `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
             data = data)
  a = summary(model)
  df[row, 1] = as.numeric(a$fstatistic[1])
  model = lm(pigment ~ `Temperature`,
             data = data)
  a = summary(model)
  df[row, 2] = as.numeric(a$fstatistic[1])
  model = lm(pigment ~ `Temperature` + I(Temperature^2),
             data = data)
  a = summary(model)
  df[row, 3] = as.numeric(a$fstatistic[1])
}

rownames(df) = colnames(pig_data)
colnames(df) = c("T + Inorganic nutrient",
                 "T",
                 "T + T^2")

write.csv(df, "lm_model.csv")


df$pigment = rownames(df)
df_long = melt(df, id.vars = "pigment")
p = ggplot(df_long, aes(x = variable, 
                        y = value,
                        fill = variable)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_col() +
  theme_bw() +
  labs(x = "Combination of variables",
       y = "F-statistic") +
  facet_wrap(~pigment, scale = "free_y") +
  theme(legend.position = c(0.9, 0.1),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))
p

ggsave("lm_model.pdf")

# CLASS
phy_ann_abu[1:10, 1:10]
length(phy_ann_abu)
fun_abu = aggregate(phy_ann_abu[4:length(phy_ann_abu)],
                    by = list(phy_ann_abu$Class),
                    sum)
rownames(fun_abu) = fun_abu$Group.1

fun_abu = as.data.frame(t(fun_abu))

fun_abu$Stations = gsub("_", "-", rownames(fun_abu))
fun_abu$Stations = gsub("-3um-[1,2,3]", "", fun_abu$Stations)
fun_abu$Stations = gsub("-0.2um-[1,2,3]", "",fun_abu$Stations)



fun_abu$Fraction = rownames(fun_abu)
fun_abu$Fraction = ifelse(grepl("0.2um", fun_abu$Fraction), 
                          "0.2-3 μm", "3-200 μm")

colnames(daya_env)
fun_abu_env = merge(daya_env[c("Stations", "Temperature")], 
                    fun_abu,
                    by = "Stations")


lf_fun_abu_env = subset(fun_abu_env, Fraction == "3-200 μm")
Sf_fun_abu_env = subset(fun_abu_env, Fraction == "0.2-3 μm")

length(lf_fun_abu_env)
head(lf_fun_abu_env)


lf_fun_abu_env[3:27] <- apply(lf_fun_abu_env[3:27], 2,
                              as.numeric)
lf_fun_abu_env[3:27] = lf_fun_abu_env[3:27]/rowSums(lf_fun_abu_env[3:27])

colSums(lf_fun_abu_env[3:27])

data = lf_fun_abu_env[3:27]

lf_fun_abu_env[3:27] = data[, order(-colSums(data))]

colnames(lf_fun_abu_env)[3:27] = colnames(data[, order(-colSums(data))])

colSums(lf_fun_abu_env[3:27])

lf_fun_abu_env = lf_fun_abu_env[1:8]
lf_fun_abu_env_long = melt(lf_fun_abu_env, id.vars = c("Stations","Temperature"))



Sf_fun_abu_env[3:27] = apply(Sf_fun_abu_env[3:27], 2, as.numeric)
Sf_fun_abu_env[3:27] = Sf_fun_abu_env[3:27]/rowSums(Sf_fun_abu_env[3:27])


data = Sf_fun_abu_env[3:27]

Sf_fun_abu_env[3:27] = data[, order(-colSums(data))]

colnames(Sf_fun_abu_env)[3:27] = colnames(data[, order(-colSums(data))])

colSums(Sf_fun_abu_env[3:27])

Sf_fun_abu_env = Sf_fun_abu_env[1:8]
Sf_fun_abu_env_long = melt(Sf_fun_abu_env, id.vars = c("Stations","Temperature"))

colnames(lf_fun_abu_env_long)
lf_fun_abu_env_long = aggregate(lf_fun_abu_env_long$value,
                                by = list(Stations = lf_fun_abu_env_long$Stations,
                                          Class = lf_fun_abu_env_long$variable,
                                          Temperature = lf_fun_abu_env_long$Temperature),
                                mean)
colnames(lf_fun_abu_env_long)
Sf_fun_abu_env_long = aggregate(Sf_fun_abu_env_long$value,
                                by = list(Stations = Sf_fun_abu_env_long$Stations,
                                          Class = Sf_fun_abu_env_long$variable,
                                          Temperature = Sf_fun_abu_env_long$Temperature),
                                mean)

colnames(Sf_fun_abu_env_long)


lf_fun_abu_env_long$Region = ''
for (i in c(1:length(lf_fun_abu_env_long$Temperature))) {
  if (lf_fun_abu_env_long$Temperature[i] <= 25) lf_fun_abu_env_long$Region[i] = "Low Temperature"
  if (lf_fun_abu_env_long$Temperature[i] > 28) lf_fun_abu_env_long$Region[i] = "High Temperature"
  if (lf_fun_abu_env_long$Temperature[i] > 25&lf_fun_abu_env_long$Temperature[i] <= 28) {
    lf_fun_abu_env_long$Region[i] = "Middle Temperature"
  }
}

lf_fun_abu_env_long$Region = factor(lf_fun_abu_env_long$Region,
                                    levels = c("Low Temperature",
                                               "Middle Temperature",
                                               "High Temperature"))


p1 = ggplot(lf_fun_abu_env_long,
            aes(x = Temperature,
                y = log2(x),
                col = Region)) +
  geom_point(alpha = 0.5) +
  #  geom_boxplot(aes(group = Region, y = x*100, col = Region)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Temperature [°C]",
       y = "log2(Relative abundance)",
       title = "3-200 μm Fraction") +
  facet_wrap(~Class, scale = "free_y") +
  #  geom_smooth(method = "glm", se = F, linewidth = 2, col = "red", formula = y ~ x+I(x^2)) +
  geom_smooth(se = T, linewidth = 2, col = "red") +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

#p1

Sf_fun_abu_env_long$Region = ''
for (i in c(1:length(Sf_fun_abu_env_long$Temperature))) {
  if (Sf_fun_abu_env_long$Temperature[i] <= 25) Sf_fun_abu_env_long$Region[i] = "Low Temperature"
  if (Sf_fun_abu_env_long$Temperature[i] > 28) Sf_fun_abu_env_long$Region[i] = "High Temperature"
  if (Sf_fun_abu_env_long$Temperature[i] > 25&Sf_fun_abu_env_long$Temperature[i] <= 28) {
    Sf_fun_abu_env_long$Region[i] = "Middle Temperature"
  }
}

p2 = ggplot(Sf_fun_abu_env_long,
            aes(x = Temperature,
                y = log2(x),
                col = Region)) +
  geom_point(alpha = 0.5) +
  #  geom_boxplot(aes(group = Region, y = x*100, col = Region)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Temperature [°C]",
       y = "log2(Relative abundance)",
       title = "0.2-3 μm Fraction") +
  facet_wrap(~Class, scale = "free_y") +
  #  geom_smooth(method = "glm", se = F, linewidth = 2, col = "Blue", formula = y ~ x+I(x^2)) +
  geom_smooth(se = T, linewidth = 2, col = "blue") +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p = p2 + p1 
p

ggsave("Daya_class.tiff", dpi = 300)

