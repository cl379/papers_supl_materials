# PCAs - grinding stones from SKP, LTS, DTR and VHV - phytoliths (phy) and starch (st)

##### SWD & PACKAGES #####
# Set working directory
setwd("F:/R")

# Install and load required packages and libraries
install.packages("FactoMineR")
library(FactoMineR)

##### LEGEND #####
# res.pca<-PCA(Sites_phy_count[,1:13] -> [,1:13] = columns read by R from the original table
# quali.sup = qualitative variable
# habillage = idem
# concat<-cbind.data.frame(Sites_phy_count[,13], res.pca$ind$coord) -> [,13] = idem
# bary = TRUE makes small confidence ellipses

##### PHYTOLITHS ALL SITES #####
Sites_phy<-read.table("Phytoliths.txt", header=TRUE)
res.pca<-PCA(Sites_phy[,1:14], scale.unit=TRUE, ncp=5, quali.sup=14, graph=FALSE)
pdf("F:/R/Graphs/Phytoliths.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=14, title="Phytoliths")
concat<-cbind.data.frame(Sites_phy[,14], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE) 
plot.PCA(res.pca, habillage=14, ellipse=ellipse.coord, cex=0.8, title="Phytoliths")
dev.off()
dimdesc(res.pca, axes=c(1,2))
print(x=res.pca, file="Phytoliths_results.csv", sep=",")

##### STARCH ALL SITES #####
Sites_st<-read.table("Starch.txt", header=TRUE)
res.pca<-PCA(Sites_st[,1:13], scale.unit=TRUE, ncp=5, quali.sup=13, graph=FALSE)
pdf("F:/R/Graphs/Starch.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=13, title="Starch")
concat<-cbind.data.frame(Sites_st[,13], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE) 
plot.PCA(res.pca, habillage=13, ellipse=ellipse.coord, cex=0.8, title="Starch")
dev.off()
dimdesc(res.pca, axes=c(1,2))
print(x=res.pca, file="Starch_results.csv", sep=",")