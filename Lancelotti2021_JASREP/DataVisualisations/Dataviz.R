### Data visualisation for publication
#### Lancelotti et al. Phytolith analyses from Khil and Kaf Taht el-Ghar (Western Maghreb) - JASREP

library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(tidypaleo)

# Load and organise datasets
EKH <- read.csv("EKH_phytoliths-new.csv")
KTG <- read.csv("KTG_phytoliths.csv")

colnames(EKH)[1] <- "sample"
colnames(KTG)[1] <- "sample"

EKH$sample <- factor(EKH$sample, levels = rev(EKH$sample))
KTG$sample <- factor(KTG$sample, levels = rev(KTG$sample))


#### Regression concentration versus number of morphotypes with linear trend and confidence interval
EKHconc <-  ggplot(EKH, aes(x=concentration, y=morphotype_n)) +
  geom_point(color="steelblue1") +
  geom_smooth(method=lm , color="red", fill="steelblue1", se=TRUE) +
  ylab("number of morphotypes") +
  scale_x_continuous(labels = comma, expand = c(0,100000))
     

KTGconc <- ggplot(KTG, aes(x=concentration, y=morphotype_n)) +
  geom_point(color="springgreen3") +
  geom_smooth(method=lm , color="red", fill="springgreen3", se=TRUE) +
  ylab("number of morphotypes") +
  scale_x_continuous(labels = comma)


#### Field values are represented by their own row and identified through the variable column
EKH.m <- melt(EKH, id.var="sample")
KTG.m <- melt(KTG, id.var="sample")

#### Subset EKH.m to isolate variables under consideration
EKH.cells <- EKH.m[33:80,]  # Proportion of elongates cells, short cells and silica skeletons
KTG.cells <- KTG.m[16:60,]

EKH.poa <- EKH.m[81:128,] # Proportion of Pooideae, Chloridoideae, Panicoideae
KTG.poa <- KTG.m[61:105,]

EKH.main <- EKH.m[129:160,] # Main groups of phytoliths (Inflorescence, Leaf/culm, Woody taxa, Palms)
KTG.main <- KTG.m[106:135,]

#### Barcharts "elongate cells", "short cells", "silica skeletons"
EKHcells <- ggplot(EKH.cells, aes(x = sample, y = as.numeric(value), fill = as.factor(variable))) + 
  geom_bar(stat = "identity", position = "fill") +
  coord_flip()+
  ylab("") +
  scale_y_continuous(label=percent) +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")

KTGcells <- ggplot(KTG.cells, aes(x = sample, y = as.numeric(value), fill = as.factor(variable))) + 
  geom_bar(stat = "identity", position = "fill") +
  coord_flip()+
  ylab("") +
  scale_y_continuous(label=percent) +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")

#### Barcharts "Pooideae", "Chloridoideae", "Panicoideae"
EKHpoa <- ggplot(EKH.poa, aes(x = sample, y = as.numeric(value), fill = as.factor(variable))) + 
  geom_bar(stat = "identity", position = "fill") +
  coord_flip()+
  ylab("") +
  scale_y_continuous(label=percent) +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")

KTGpoa <- ggplot(KTG.poa, aes(x = sample, y = as.numeric(value), fill = as.factor(variable))) + 
  geom_bar(stat = "identity", position = "fill") +
  coord_flip()+
  ylab("") +
  scale_y_continuous(label=percent) +
  theme(legend.title=element_blank(), legend.position = "bottom") +
    scale_fill_brewer(palette = "Paired")


#### Stratigraphic changes in main groups
EKH.m1 <- mutate(EKH.m, variable = recode(variable, "inflorescence" = "Inflorescence", 
                                          "leaf.culm" = "Leaf/Culm",
                                          "woody_taxa" = "Woody Taxa", 
                                          "palms" = "Palms"))
EKH.facet <- rbind(
  cbind(subset(EKH.m1, EKH.m1$variable == "Inflorescence"), period = EKH$period),
  cbind(subset(EKH.m1, EKH.m1$variable == "Leaf/Culm"), period = EKH$period),
  cbind(subset(EKH.m1, EKH.m1$variable == "Woody Taxa"), period = EKH$period),
  cbind(subset(EKH.m1, EKH.m1$variable == "Palms"), period = EKH$period))
EKH.facet$value <- as.numeric(EKH.facet$value)
EKH.facet$period <- factor(EKH.facet$period, 
                           levels = c("Post-Neolithic", "Middle-Neolithic", "Early-Neolithic-B", "Early-Neolithic-A"))
colnames(EKH.facet)[4] <- "Period"

EKH.facet.plot <- ggplot(EKH.facet, aes(x = value, y = sample, group = "all")) + 
  geom_lineh(position = "identity") +
  geom_point(aes(color = Period), size = 4) +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL, y = NULL) 

KTG.m1 <- mutate(KTG.m, variable = recode(variable, "inflorescence" = "Inflorescence", 
                                          "leaf.culm" = "Leaf/Culm",
                                          "woody_taxa" = "Woody Taxa", 
                                          "palms" = "Palms"))
KTG.facet <- rbind(
  cbind(subset(KTG.m1, KTG.m1$variable == "Inflorescence"), period = KTG$period),
  cbind(subset(KTG.m1, KTG.m1$variable == "Leaf/Culm"), period = KTG$period),
  cbind(subset(KTG.m1, KTG.m1$variable == "Woody Taxa"), period = KTG$period),
  cbind(subset(KTG.m1, KTG.m1$variable == "Palms"), period = KTG$period))
KTG.facet$value <- as.numeric(KTG.facet$value)
KTG.facet$period <- factor(KTG.facet$period, 
                           levels = c("Historic", "Middle-Neolithic", "Early-Neolithic-B", "Early-Neolithic-A", 
                                      "Transition", "Epipaleolithic", "Paleolithic"))
colnames(KTG.facet)[4] <- "Period"

KTG.facet.plot <- ggplot(KTG.facet, aes(x = value, y = sample, group = "all")) + 
  geom_lineh(position = "identity") +
  geom_point(aes(color = Period), size = 4) +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL, y = NULL) 
