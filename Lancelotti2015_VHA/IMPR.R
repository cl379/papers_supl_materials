#####################################################
### Script used for analyses related to the paper ###
### Lancelotti et al. Investigating the use of fuel through a combination of phytoliths and multi-element analysis. An ethnographic experiment
### Authors:    Carla Lancelotti - CaSES UPF - carla.lancelotti@upf.edu
###             Javier Ruiz-Pérez - CaSES and UAB - jruizperez1991@gmail.com
### Created:    01 June 2015
### URL:        https://github.com/cl379/papers_supl_materials/tree/master/Lancelotti2015impr
#####################################################

require(FactoMineR)
require(MASS)
require(ggplot2)
require(plyr)
require(reshape2)


setwd("~/Documents/R/JAN/phy")
phy <- read.table("morphotype_percent.txt", header = TRUE, dec = ",", row.names = 1)
str(phy)

#### Descriptive statistics ####
summary(phy)

### Subset of the three fireplaces with their control ###
fp.in.r <- phy[ which(phy$Sample_1 == 'Inside_right'), ]
fp.in.l <- phy[ which(phy$Sample_1 == 'Inside_left'), ]
fp.out <- phy[ which(phy$Sample_1 == 'Veranda'), ]
str(fp.in.r)
str(fp.in.l)
str(fp.out)

### Boxplots of concentration ###
ggplot(phy, aes(x=Sample_1, y=Concentration, fill=Sample_1)) + geom_boxplot() + theme(axis.title.x = element_blank()) + theme(legend.title=element_blank())
ggplot(phy, aes(x=Type, y=Concentration, fill=Sample_1)) + geom_boxplot() + theme(axis.title.x = element_blank()) + theme(legend.title=element_blank())

#### Multivariate statistics ####
### MANOVA all samples, all groups ###
X <- as.matrix(phy[, 4:28])
man <- manova(X ~ Type*Sample_1, data = phy) ### Type III SS for Type ###
summary(man, "Pillai")
X <- as.matrix(phy[, 4:28]) ### Type III SS for Sample_1 ###
> man <- manova(X ~ Sample_1*Type, data = phy)
> summary(man)

### MANOVA fireplaces vs. control samples ###
X <- as.matrix(fp.in.r[2:7, 4:28]) # Fireplace Inside_right
aov <- aov(X ~ Type, data = fp.in.r[2:7, ])
summary(aov)

X <- as.matrix(fp.in.l[2:7, 4:28]) # Fireplace Inside_left
aov <- aov(X ~ Type, data = fp.in.l[2:7, ])
summary(aov)

X <- as.matrix(fp.out[2:7, 4:28]) # Fireplace Veranda
aov <- aov(X ~ Type, data = fp.out[2:7, ])
summary(aov)

####Perform the PCA (developed for datasets with qualitative/nominal supplementary variable/s)
result_dataset<-PCA(dataset, quali.sup=0, graph=FALSE) #Put the column/s number/s of the qualitative/categorical variable/s in quali.sup

pdf("URL_FOR_THE_OUTPUT.pdf") #Export the graphics in a single PDF file
plot.PCA(result_dataset, axes=c(1,2), choix="var") #Specify the components to plot in axes
concat<-cbind.data.frame(dataset[0], result_dataset$ind$coord) #Substitute the "0" with the same value entered in quali.sup
ellipse.coord<-coord.ellipse(concat, bary=TRUE) #bary=TRUE calculates the coordinates of the ellipses around the barycentre of individuals
plot.PCA(result_dataset, habillage=0, ellipse=ellipse.coord) #Substitute the "0" with the same value entered in quali.sup
dev.off() #End the pdf function

print(x=result_dataset, file="URL_FOR_THE_OUTPUT.csv", sep=";") #Export the PCA results in a CSV file (columns separated by ";")

