# Install and load required packages and libraries
install.packages("FactoMineR")
library(FactoMineR)

# set wd
setwd("~/Documents/R/phd/percentiles")

# Import dataset (all data were transformed into percentiles and ranked)
charc.site <- read.table("all_charc.csv", header = TRUE)
phy.site <- read.table("all_phyt.csv", header = TRUE)
chem.site <- read.table("all_chem.csv", header = TRUE)

dun.phy <- read.table("dung_phyt.csv", header = TRUE)
dun.chem <- read.table("dung_chem.csv", header = TRUE)

har.char <- read.table("harp_charc.csv", header = TRUE)
har.phy <- read.table("harp_phyt.csv", header = TRUE)
har.chem <- read.table("harp_chem.csv", header = TRUE)

kmr.char <- read.table("kmr_charc.csv", header = TRUE)
kmr.phy <- read.table("kmr_phyt.csv", header = TRUE)
kmr.chem <- read.table("kmr_chem.csv", header = TRUE)

skp.char <- read.table("skp_charc.csv", header = TRUE)
skp.phy <- read.table("skp_phyt.csv", header = TRUE)
skp.chem <- read.table("skp_chem.csv", header = TRUE)

alm.char <- read.table("alm_charc.csv", header = TRUE)
alm.phy <- read.table("alm_phyt.csv", header = TRUE)
alm.chem <- read.table("alm_chem.csv", header = TRUE)



##### PCA #####


### Charcoal ###
 
## PCA on charcoal data per site (all sites) ##
res.pca<-PCA(charc.site[,1:8], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/charcoal_sites.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Anthraco-remains")
concat<-cbind.data.frame(charc.site[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Anthraco-remains")
mtext("", outer=TRUE, line=-1.5, font=2, adj=0.5, cex=1.25)
dev.off()

# Harappa charcoal
res.pca<-PCA(har.char[,1:6], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/charcoal_harappa.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Harappa Anthraco-remains")
concat<-cbind.data.frame(har.char[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Harappa Anthraco-remains")
dev.off()

# Kanmer charcoal
res.pca<-PCA(kmr.char[,1:7], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/charcoal_kanmer.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Kanmer Anthraco-remains")
concat<-cbind.data.frame(kmr.char[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Kanmer Anthraco-remains")
dev.off()

# Shikarpur charcoal
res.pca<-PCA(skp.char[,1:6], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/charcoal_shikarpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Shikarpur Anthraco-remains")
concat<-cbind.data.frame(skp.char[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Shikarpur Anthraco-remains")
dev.off()

# Alamgirpur charcoal
res.pca<-PCA(alm.char[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/charcoal_alamgirpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Alamgirpur Anthraco-remains")
concat<-cbind.data.frame(alm.char[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Alamgirpur Anthraco-remains")
dev.off()

### Phytoliths ###

## All sites (per site) ##
res.pca<-PCA(phy.site[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_sites.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths")
concat<-cbind.data.frame(phy.site[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths")
dev.off()

# Dung phytoliths
res.pca<-PCA(dun.phy[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_dung.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths dung")
concat<-cbind.data.frame(dun.phy[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths dung")
dev.off()

# Harappa phytoliths
res.pca<-PCA(har.phy[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_harappa.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths Harappa")
concat<-cbind.data.frame(har.phy[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths Harappa")
dev.off()

# Kanmer phytoliths
res.pca<-PCA(kmr.phy[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_Kanmer.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths")
concat<-cbind.data.frame(kmr.phy[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths Kanmer")
dev.off()

# Shikarpur phytoliths
res.pca<-PCA(skp.phy[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_shikarpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths Shikarpur")
concat<-cbind.data.frame(skp.phy[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths Shikarpur")
dev.off()

# Alamgirpur phytoliths
res.pca<-PCA(alm.phy[,1:5], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/phytoliths_alamgirpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Phytoliths Alamgirpur")
concat<-cbind.data.frame(alm.phy[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Phytoliths Alamgirpur")
dev.off()


### Chemicals ###

## All sites (per site) ##
res.pca<-PCA(chem.site[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_all.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement")
concat<-cbind.data.frame(chem.site[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement")
dev.off()


# Dung multielement
res.pca<-PCA(dun.chem[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_dung.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement dung")
concat<-cbind.data.frame(dun.chem[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement dung")
dev.off()

# Harappa multielement
res.pca<-PCA(har.chem[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_harappa.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement Harappa")
concat<-cbind.data.frame(har.chem[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement Harappa")
dev.off()

# Kanmer multielement
res.pca<-PCA(kmr.chem[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_kanmer.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement Kanmer")
concat<-cbind.data.frame(kmr.chem[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement Kanmer")
dev.off()

# Shikarpur multielement
res.pca<-PCA(skp.chem[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_shikarpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement Shikarpur")
concat<-cbind.data.frame(skp.chem[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement Shikarpur")
dev.off()

# Alamgirpur multielement
res.pca<-PCA(alm.chem[,1:36], scale.unit=TRUE, ncp=5, quali.sup=1, graph=FALSE)
pdf("~/Documents/R/phd/percentiles/plots2/chemical_alamgirpur.pdf")
plot.PCA(res.pca, axes=c(1, 2), choix="var", habillage=1, title="Multielement Alamgirpur")
concat<-cbind.data.frame(alm.chem[,1], res.pca$ind$coord)
ellipse.coord<-coord.ellipse(concat, bary=TRUE)
plot.PCA(res.pca, habillage=1, ellipse=ellipse.coord, cex=0.8, title="Multielement Alamgirpur")
dev.off()

