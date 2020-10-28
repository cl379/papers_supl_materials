#### Code used for the analysis of chemical elements measured by pXRF
#### Associated to the paper:
#### Biagetti et al. 2020. "Using pXRF on anthropogenic sediments to  map geo-chemical signatures at Seoke (Botswana)"
#### Submitted to PLOS ONE.
#### Author of the code: Jonas Alcaina-Mateos (jonas.alcaina@upf.edu)

#### All datasets used in this analysis are provided with the paper


### Load required packages 
require(sp)
require(gstat)
require(maptools)
require(rgdal)
require(raster)
require(ggplot2)
require(grid)    
require(gridExtra)
require(RColorBrewer)

### Load Goulard & Voltz (1992) algorithm    
source("./data/GV_functions.R")

### Load dataset 
Data <- read.csv("./data/dataset.csv", row.names=1, na.strings = "LOD")
Data[is.na(Data)] <- 0 # set "LOD" to 0 
# create a transformed dataset matrix 
dataset <- scale(log1p(Data[,10:ncol(Data)])) # (z-scores of log-transformed values)

### Load control data
ControlData <- read.csv("./data/control.csv", row.names=1, na.strings = "LOD")
ControlData[is.na(ControlData)] <- 0

### Load additional data 
load("./data/Grid.RData")   # real distanced grid 
load("./data/GridM.RData")  # modified distanced grid (for plot) 
load("./data/Topo.RData")   # topography 

### Build spatial objects 
# build Spatial Points 
spdf <- Data
coordinates(spdf) <- ~X+Y
spdf@proj4string <- CRS("+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs")
spdf@data <- as.data.frame(dataset) # add transformed values
# build gstat object 
g <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf)
for (i in 2:length(names(spdf))) g <- gstat(g, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf)
# compute empirical variogram 
v <- variogram(g, cutoff=20, width=2)

### Fit linear model of corregionalization (LMC)
# fill initial model
g <- gstat(g, model = rbind(vgm(0.5, "Sph", range=5, 0), vgm(0.5, "Exp", range=15)), fill.all = T)
# fit by Goulard & Voltz algorithm (see GV_functions.R file)
Model <- fitGV(v, g$model)
# plot variogram and LMC
plot(v, model=Model, ylim=c(-0.75, 1.5), cex=0.4)

### Build separate gstat objects for each of the areas
gA <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Data$area == "A"),])
for (i in 2:length(names(spdf))) gA <- gstat(gA, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Data$area == "A"),])
gA$model <- Model
gB <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Data$area == "B"),])
for (i in 2:length(names(spdf))) gB <- gstat(gB, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Data$area == "B"),])
gB$model <- Model
gC <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Data$area == "C"),])
for (i in 2:length(names(spdf))) gC <- gstat(gC, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Data$area == "C"),])
gC$model <- Model

### Predictions by Co-Kriging 
ckA <- predict(gA, newdata = Grid[Grid@data$id == "A",])
ckB <- predict(gB, newdata = Grid[Grid@data$id == "B",])
ckC <- predict(gC, newdata = Grid[Grid@data$id == "C",])

### Principal Component Analysis
# obtain marix of variance-covariance from the LMC parameters
S <- buildMat(rowSums(t(sapply(names(Model), function(x) Model[[x]]$psill))))
# apply pca
ei <- eigen(S)
U <- ei$vectors
L <- diag(ei$values)
expVar <- diag(L) / sum(diag(L)) # explained variance

### PCA plot of archaeological samples 
# set objects for plot 
dPoints <- as.data.frame(dataset %*% U %*% solve(L^0.5))
dArrows <- as.data.frame(U %*% L^0.5)
Const <- max(abs(range(dArrows[,1:2]))) / max(abs(range(dPoints[,1:2])))
SIZ.aux <- ifelse(Data$element == "undetermined", 1, 2)
# build plot with ggplot2
pcaplot1 <- ggplot(data=dArrows, aes(x=V1, y=V2)) + 
    geom_hline(yintercept = 0, colour = "gray50", linetype = "dashed") + 
    geom_vline(xintercept = 0, colour = "gray50", linetype = "dashed") + 
    geom_point(data= dPoints, aes(x=Const*V1, y=Const*V2, colour=Data$area, fill=Data$area, size=SIZ.aux, shape=Data$element)) + 
    scale_y_continuous(sec.axis = sec_axis(~ . * Const^-1)) + 
    scale_x_continuous(sec.axis = sec_axis(~ . * Const^-1)) + 
    geom_segment(aes(x = 0, y = 0, xend = V1, yend = V2), 
            arrow = arrow(length = unit(1/2, 'picas')), color = "royalblue4", size=0.6) + 
    geom_text(aes(x=V1*1.075, y=V2*1.075, label = names(spdf@data), fontface = "bold"), color = "royalblue4", size=4) + 
    labs(x = paste0("PC1 (", round(expVar[1]*100, 2), "% explained var.)"), 
            y = paste0("PC2 (", round(expVar[2]*100, 2), "% explained var.)")) + 
    theme_classic(base_size=12) + 
    coord_cartesian(xlim = c(-max(abs(dArrows[,1:2]))*1.02, max(abs(dArrows[,1:2]))*1.02), 
            ylim = c(-max(abs(dArrows[,1:2]))*1.02, max(abs(dArrows[,1:2]))*1.02)) + 
    scale_colour_brewer(palette = "Set1") + 
    scale_fill_brewer(palette = "Set1") +    
    scale_size(range = c(1, 4)) + 
    scale_shape_manual(values = c(24, 25, 7, 12, 10, 13, 15, 16)) + 
    guides(fill = FALSE, size = FALSE) + 
    labs(color = "Area", shape="Element") + 
    theme(legend.title = element_text(face="bold")) + 
    theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="gray24"))
    
### Project the control samples on PCA space 
# transform control data
Control <- log1p(ControlData[,6:ncol(ControlData)])
Control <- sapply(1:ncol(Control), function(x) (Control[,x] - attr(dataset, "scaled:center")[x]) / attr(dataset, "scaled:scale")[x])
# set objects for plot 
dPoints1 <- as.data.frame(dataset %*% U)
dPoints2 <- as.data.frame(Control %*% U)
# build plot with ggplot2
pcaplot2 <- ggplot(dPoints1) +
    geom_point(aes(x = V1, y = V2), size=2, color="blue", alpha=0.2) + 
    geom_point(data=dPoints2, aes(x = V1, y = V2), size=3, color="red1") + 
    geom_text(data=dPoints2, aes(label = row.names(ControlData), x = V1, y = V2), size = 4, angle = 0, vjust = 2) + 
    labs(x = paste0("PC1 (", round(expVar[1]*100, 2), "% explained var.)"), 
            y = paste0("PC2 (", round(expVar[2]*100, 2), "% explained var.)")) + 
    theme_linedraw(base_size=12)

### Plot the Co-Kriging interpolation 
# bind the results from the three areas
df <- rbind(ckA@data[,seq(1, 34, 2)], ckB@data[,seq(1, 34, 2)], ckC@data[,seq(1, 34, 2)])
colnames(df) <- names(spdf)
# build a Spatial Pixels object with a modified-distance grid
CK <- SpatialPixelsDataFrame(as(GridM, 'SpatialPixels'), data=df)
# compute variable breaks for plot 
brks <- lapply(1:17, function(x) seq(min(CK@data[,x])-0.1, max(CK@data[,x])+0.1, length.out=9))
# obtain a list of individual plots 
LIST <- lapply(order(colnames(df)), function(x) {
    spplot(CK[,x], cex=10, 
        par.settings=list(axis.line=list(lwd=4), strip.border=list(lwd=4)), 
        par.strip.text=list(cex=5),
        main=textGrob(colnames(df)[x], gp=gpar(fontsize=40,font=4)),
        sp.layout=list(Topo[[1]], Topo[[2]], Topo[[3]], first=F, lwd = 3), 
        col.regions = unlist(lapply(rev(brewer.pal(8, 'RdYlBu')), function(y) rep(y, 2))),  
        at=brks[[x]],
        colorkey = list(at=brks[[x]], labels=list(at=brks[[x]], labels=round(brks[[x]],1), cex=2), width=3))
})
# build a panel with all plots 
grid.arrange(LIST[[1]], LIST[[2]], LIST[[3]], LIST[[4]], LIST[[5]], LIST[[6]], LIST[[7]], 
             LIST[[8]], LIST[[9]], LIST[[10]], LIST[[11]], LIST[[12]], LIST[[13]], 
             LIST[[14]], LIST[[15]], LIST[[16]], LIST[[17]], ncol=5)

