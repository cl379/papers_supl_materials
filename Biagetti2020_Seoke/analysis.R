### Load required packages 
require(sp)
require(gstat)
require(maptools)
require(rgdal)
require(raster)
require(rasterVis)
require(ggplot2)
require(grid)    
require(gridExtra)
require(latticeExtra)
require(easyCODA)
require(ggpubr)
require(stringr) 
require(inlmisc)
source("./Robjects/GV_functions.R")

### ====================================================================================
### Load and prepare data 
### ====================================================================================

### Load data --------------------------------------------------------------------------
# load dataset 
Data <- read.csv("./BOTSWANA_dataset.csv", na.strings = "< LOD")
# load spatial points data 
Points <- readOGR("./shp/Sample_grid.shp")
# load spatial grid for interpolatoin 
load("./Robjects/Grid.RData")
# load topographic maps (to add in the plots) 
load("./Robjects/SHP.RData")

### repeated samples -------------------------------------------------------------------
# Average the repeated measures of samples
# identify repeated measures 
Index <- str_sub(Data$SAMPLE, -1, -1) %in% c("a", "b", "c", "d")
# split by sample 
RepeatedList <- split(Data[Index,], f=str_sub(Data[Index,]$SAMPLE, 1, -2))
# average 
Averaged <- do.call("rbind", lapply(1:length(RepeatedList), function(x) {
                aux <- cbind(RepeatedList[[x]][1,1:12], t(colMeans(RepeatedList[[x]][,13:90])))
                aux$SAMPLE <- names(RepeatedList)[x]
                return(aux)}))
# add to data 
Data <- rbind(Averaged, Data[!Index,])
# sort by Reading.No
Data <- Data[order(Data$Reading.No),] 

### Select varables --------------------------------------------------------------------
# selection based on the error (RDS) provided by the device 
# "Bal" refers to the "residual" part (aso reported by the device) 
VARS <- c("Al", "Si", "P", "S", "Cl", "K", "Ca", "Ti", "Cr", "Mn", "Fe", "Zr", "Bal")

### <LOD replacement -------------------------------------------------------------------
# repplace LOD to half the detection limit (only for selected variables) 
for (i in VARS) Data[is.na(Data[,i]),i] <- min(Data[,i], na.rm=T)/2

### divide grid and control samples ----------------------------------------------------
# grid samples (match to points data)
Dataset <- Data[match(Points@data$OBJNAME, Data$SAMPLE),VARS]
rownames(Dataset) <- Points@data$OBJNAME
# control samples 
Control <- Data[478:nrow(Data),VARS]
rownames(Control) <- Data[478:nrow(Data),]$SAMPLE

### add misisng values to Bal ----------------------------------------------------------
# Sum up to 100 and add to "Bal" for each row 
Dataset$Bal <- Dataset$Bal + (100 - rowSums(Dataset))
Control$Bal <- Control$Bal + (100 - rowSums(Control))

### raw summary statistics -------------------------------------------------------------
# summary statistics
write.csv(summary(Dataset[,1:(ncol(Dataset)-1)]), file="summary.csv")
write.csv(summary(Control[,1:(ncol(Dataset)-1)]), file="summary_control.csv")
# correlation 
aux <- round(cor(Dataset[,1:(ncol(Dataset)-1)]), 3)
write.csv(aux, file="correlation.csv")

### ====================================================================================
### LOG-RATIO analysis
### ====================================================================================

### Archaeological samples--------------------------------------------------------------

### exclude Bal and re-clausure 
Dataset.aux <- t(apply(Dataset[,1:(ncol(Dataset) -1)], 1, function(x) x / sum(x)))

### Unweighted LDA
lra <- LRA(Dataset.aux, weight=F)
PLOT.LRA(lra, map="contribution")

### loadings table 
aux <- lra$colcoord
rownames(aux) <- colnames(Dataset.aux)
aux <- round(aux[,1:5], 4)
write.csv(aux, file="lra_loadings.csv")

### prepare objects for contribution biplot
Col <- rep(1 / ncol(Dataset.aux), ncol(Dataset.aux))
expVar <- lra$sv^2 / sum(lra$sv^2)
dPoints <- as.data.frame(lra$rowcoord[,1:2] %*% diag(lra$sv[1:2]))
dArrows <- as.data.frame(diag(Col)^0.5 %*% lra$colcoord[,1:2])
sel <- which(Points@data$element == "undetermined" | Points@data$element == "enclosure_2")
Const <- 2

### build plot with ggplot2 (elipses) 
pcaplot1 <- ggplot(data=dArrows, aes(x=V1, y=V2)) + 
    geom_hline(yintercept = 0, colour = "gray50", linetype = "dashed") + 
    geom_vline(xintercept = 0, colour = "gray50", linetype = "dashed") + 
    geom_segment(aes(x = 0, y = 0, xend = V1, yend = V2), 
            arrow = arrow(length = unit(1/2, "picas")), color = "gray60", size=0.6) + 
    geom_point(data = dPoints, aes(x=Const*V1, y=Const*V2), alpha=0.1) + 
    stat_conf_ellipse(data = dPoints[-sel,], 
            aes(x=Const*V1, y=Const*V2, colour=Points@data$element[-sel], fill=Points@data$element[-sel]), 
            level=0.95, geom = "polygon", alpha = 0.2) + 
    geom_text(aes(x=V1*1.075, y=V2*1.075, label = colnames(Dataset.aux), 
                fontface = "bold"), color = "gray20", size=4) + 
    labs(x = paste0("LRA dimension 1 (", round(expVar[1]*100, 2), "%)"), 
            y = paste0("LRA dimension 2 (", round(expVar[2]*100, 2), "%)")) +  
    scale_y_continuous(sec.axis = sec_axis(~ . * Const^-1)) + 
    scale_x_continuous(sec.axis = sec_axis(~ . * Const^-1)) + 
    coord_fixed() + 
    scale_colour_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette = "Dark2") +  
    theme_bw() +     
    labs(colour = "Feature", fill="Feature") + 
    theme(legend.title = element_text(face="bold")) + 
    theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="gray24")) + 
    theme(legend.position = c(0.1, 0.75)) 

### Control samples --------------------------------------------------------------------

### exclude Bal and re-clausure 
Control.aux <- t(apply(Control[,1:(ncol(Dataset) -1)], 1, function(x) x / sum(x)))
DATA <- rbind(Dataset.aux, Control.aux)

### Unweighted LDA
lra <- LRA(DATA, weight=F)
PLOT.LRA(lra, map="contribution")

### prepare objects for contribution biplot
Col <- rep(1 / ncol(DATA), ncol(DATA))
expVar <- lra$sv^2 / sum(lra$sv^2)
dPoints <- as.data.frame(lra$rowcoord[,1:2] %*% diag(lra$sv[1:2])) 
dArrows <- as.data.frame(diag(Col)^0.5 %*% lra$colcoord[,1:2])
dPoints1 <- dPoints[1:nrow(Dataset),]
dPoints2 <- dPoints[478:522,]

### build plot with ggplot2
pcaplot2 <- ggplot(dPoints1) +
    geom_point(aes(x = V1, y = V2), size=2, color="blue", alpha=0.1) + 
    geom_segment(data=dArrows, aes(x = 0, y = 0, xend = V1, yend = V2), 
            arrow = arrow(length = unit(1/2, "picas")), color = "gray50", size=0.6) + 
    geom_text(data=dPoints2, aes(label = row.names(Control), x = V1, y = V2), size = 3, fontface = "italic") + 
    geom_text(data=dArrows, aes(x=V1*1.075, y=V2*1.075, label = colnames(Dataset.aux), fontface = "bold"), 
            color = "gray20", size=4) + 
    labs(x = paste0("LRA dimension 1 (", round(expVar[1]*100, 2), "%)"), 
            y = paste0("LRA dimension 2 (", round(expVar[2]*100, 2), "%)")) +  
    theme_bw()     

### ====================================================================================
### Geostats
### ====================================================================================

### clr transformation
# "Bal" nedded to backtransform the interpolated composition 
CLRm <- t(apply(Dataset, 1, function(x) log(x / prod(x)^(1/length(x)))))

### build spatial object (using CLR as data) 
spdf <- cbind(Points@data[,2:3], as.data.frame(CLRm))
coordinates(spdf) <- ~X+Y
spdf@proj4string <- CRS("+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs")

### build gstat object 
g <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf)
for (i in 2:length(names(spdf))) g <- gstat(g, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf)
# compute empirical variogram 
v <- variogram(g, cutoff=20, width=2)

### Fit linear model of corregionalization (LMC)
# fill initial model
g <- gstat(g, model = rbind(vgm(0.5, "Sph", range=7, 0), vgm(0.5, "Exp", range=10)), fill.all = T)
# fit by Goulard & Voltz algorithm (see GV_functions.R file)
Model <- fitGV(v, g$model)
# plot variogram and LMC
plot(v, model=Model, cex=0.4)

### Build separate gstat objects for each area (using the same geostatistical model) 
gA <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Points@data$area == "A"),])
for (i in 2:length(names(spdf))) gA <- gstat(gA, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Points@data$area == "A"),])
gA$model <- Model
gB <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Points@data$area == "B"),])
for (i in 2:length(names(spdf))) gB <- gstat(gB, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Points@data$area == "B"),])
gB$model <- Model
gC <- gstat(id=names(spdf)[1], formula=as.formula(paste(names(spdf)[1], " ~ 1")), data=spdf[which(Points@data$area == "C"),])
for (i in 2:length(names(spdf))) gC <- gstat(gC, names(spdf)[i], as.formula(paste(names(spdf)[i], " ~ 1")), spdf[which(Points@data$area == "C"),])
gC$model <- Model

### predict by Co-Kriging
ckA <- predict(gA, newdata = Grid[Grid@data$id == "A",])
ckB <- predict(gB, newdata = Grid[Grid@data$id == "B",])
ckC <- predict(gC, newdata = Grid[Grid@data$id == "C",])
CKlist <- list(A=ckA, B=ckB, C=ckC)

### ====================================================================================
### Mapping 
### ====================================================================================

### prepare data 
# bind the results of the three areas and back-transform 
dfb <- rbind(invCLR(ckA@data[,seq(1, ncol(Dataset)*2, 2)]),
             invCLR(ckB@data[,seq(1, ncol(Dataset)*2, 2)]),
             invCLR(ckC@data[,seq(1, ncol(Dataset)*2, 2)]))
# express in percentage
dfb <- dfb * 100

# set breaks for colorkey
brks <- lapply(1:ncol(dfb), function(x) seq(min(dfb[,x]), max(dfb[,x]), length.out=40))
# set breaks for labels  
brksL <- lapply(1:ncol(dfb), function(x) {
                        aux <- round(seq(min(dfb[,x]), max(dfb[,x]), length.out=5), 4)
                        aux[1] <- aux[1] + 0.0001; aux[length(aux)] <- aux[length(aux)] - 0.0001; 
                        return(aux)})

### spatial plot -----------------------------------------------------------------------
# set zone ("A", "B" or "C")
SET <- "A"
# values for plot
df <- invCLR(CKlist[[SET]]@data[,seq(1, ncol(Dataset)*2, 2)])
df <- df * 100
colnames(df) <- names(spdf)
# build spatial object 
CKbrk <- brick(SpatialPixelsDataFrame(as(Grid[Grid$id == SET,], "SpatialPixels"), data=df))
# list of plots, one for each variable  
LIST <- lapply(order(colnames(df)), function(x) {   
    levelplot(CKbrk[[x]], 
        at=brks[[x]], 
        col.regions = GetColors(255, scheme = "sunset"), 
        colorkey = list(at=brks[[x]], labels=list(cex=1.5, at=brksL[[x]], labels=brksL[[x]]), width=3, height=0.7),        
        par.settings=list(axis.line=list(lwd=3), strip.border=list(lwd=3)), 
        par.strip.text=list(cex=5),
        margin=F, main = textGrob(colnames(df)[x], gp=gpar(fontsize=30,font=4))) +  
    latticeExtra::layer(sp.points(SHP, col = "black", lwd = 3))
})
# build a panel with all plots 
grid.arrange(LIST[[1]], LIST[[3]], LIST[[4]], LIST[[5]], LIST[[6]], LIST[[7]], 
             LIST[[8]], LIST[[9]], LIST[[10]], LIST[[11]], LIST[[12]], LIST[[13]], ncol=4)


