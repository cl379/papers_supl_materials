library(sf)
library(sp)
library(rgdal)
library(dplyr)
library(vegan)
library(raster)
library(readxl)
library(maptools)
library(textshape)
library(spatialEco)
library(spatial.tools)
library(PresenceAbsence)

#Extract spatial data for analysis

setwd() #Set directory with raster data.
list <- list.files(path = "D:/", pattern = ".tif$", all.files = TRUE, full.names = TRUE) # Import rasters
r0 <- raster("D:/") # Load reference raster
ext <- as(extent(-20, 150, -35, 60), "SpatialPolygons") # Select reference extent for study area
crs(ext) <- crs(r0) # Establish CRS for reference extent
r0_ext <- crop(r0, ext) # Crop reference raster to match extent

# Raster homogenization:
rlist_h <- list() 
for (i in 1:length(rlist)) { 
  print(rlist[[i]])
  ri <- raster(rlist[[i]])
  crs(ri) <- crs(r0)
  if (any(res(ri)) != any(res(r0))) {
    ri_r <- resample(ri, r0, "bilinear") # Resample to the same grid
    ri_c <- crop(ri_r, ext) # Set raster extents to match by cropping them w/ new spatial object (ext)
    ri_m <- mask(ri_c, r0_ext) # Remove data which falls outside one of the rasters (e.g. oceans)
    rlist_h[[i]] <- ri_m
  }
  else {
    ri_c <- crop(ri, ext) # Set raster extents to match by cropping them w/ new spatial object (ext)
    ri_m <- mask(ri_c, r0_ext) # Remove data which falls outside one of the rasters (e.g. oceans)
    rlist_h[[i]] <- ri_m
  }
} 

# Extract spatial data by XY coordinates
data <- read_excel("D:/")
a <- textshape::column_to_rownames (data, loc=1)
xy <- data.frame(a[1:2]) # Load XY coordinates
coordinates(xy) <- cbind(xy$LONG, xy$LATI) # Convert into spatial points
crs(xy) <- crs(r0) # Set CRS (units of "buffer" in extract depends on chosen CRS)

b <- data.frame()
b <- textshape::column_to_rownames (data.frame(rownames(a)), loc=1)
for (i in 1:length(rlist_h)) {
  print(i)
  print(rlist_h[[i]])
  ri <- rlist_h[[i]]
  e1 <- raster::extract(ri, xy, buffer = 50000, fun = mean, na.rm = TRUE)
  e2 <- raster::extract(ri, xy, buffer = 50000, fun = var, na.rm = TRUE)
  b <- cbind(b, (as.data.frame(cbind(e1, e2), row.names = rownames(xy))))
  colnames(b)[names(b) == c("e1", "e2")] <- c(paste((names(rlist_h[[i]])), "_m", sep = ""),
                                              paste((names(rlist_h[[i]])), "_v", sep = ""))
}
print(b)

# Extract spatial data successive rasters by all polygons:
data <- read_excel("D:/") # Load data frame w/ cultures as cases
a <- textshape::column_to_rownames (data, loc=1) # Set cultures names as rows
shp <- st_read("D:/") # Load .shp file
crs(shp) <- crs(r0) # Establish CRS for polygons

c2 <- data.frame
c2 <- textshape::column_to_rownames (data.frame(rownames(a)), loc=1)
for (i in 1:length(rlist_h)) {
  print(rlist_h[[i]])
  ri <- rlist_h[[i]]
  b2 <- data.frame()
  for (j in 1:nrow(shp)) {
    print(j)
    pi <- remove.holes(shp[j,])
    e3 <- raster::extract(ri, pi, fun = mean, na.rm = TRUE)
    e4 <- raster::extract(ri, pi, fun = var, na.rm = TRUE)
    b2 <- rbind(b2, as.data.frame(cbind(e3, e4)))
  }
  colnames(b2)[names(b2) == c("V1", "V2")] <- c(paste((names(rlist_h[[i]])), "_m", sep = ""),
                                                paste((names(rlist_h[[i]])), "_v", sep = ""))
  c2 <- cbind(c2, b2)
}
print(c2)


#Data modelling

setwd("C:/Users/User/Downloads/") #Set working directory with SI files

#Load and transform data
##Training dataset
train.dataset <- read.csv("1-eHRAF-dataset.csv")
data <- textshape::column_to_rownames(train.dataset[,-1], loc = 1)

df.cro <- data[, colnames(data)[c(4:6, 21:22, 29:ncol(data))]] #Selecting crop data.
colnames(df.cro)[1:3] <- c("FM", "PM", "SB")
cro_trainMean <- apply((df.cro[,6:ncol(df.cro)]), 2, mean) #Find training data mean by column
cro_trainSD <- apply((df.cro[,6:ncol(df.cro)]), 2, sd) #Find training data SD by column
df.cro[,6:ncol(df.cro)] <- scale(df.cro[,6:ncol(df.cro)]) #Scale training explanatory data

df.fm <- data[data$FM == 1, colnames(data)[c(7:10, 21:24, 29:ncol(data))]] #Selecting FM cultivation data.
colnames(df.fm)[1:4] <- c("EXT", "INT", "RF", "IRR")
fm_trainMean <- apply((df.fm[,7:ncol(df.fm)]), 2, mean) #Find training data mean by column
fm_trainSD <- apply((df.fm[,7:ncol(df.fm)]), 2, sd) #Find training data SD by column
df.fm[,7:ncol(df.fm)] <- scale(df.fm[,7:ncol(df.fm)]) #Scale training explanatory data

df.pm <- data[data$PM == 1, colnames(data)[c(11:14, 21:22, 25:26, 29:ncol(data))]] #Selecting FM cultivation data.
colnames(df.pm)[1:4] <- c("EXT", "INT", "RF", "IRR")
pm_trainMean <- apply((df.pm[,7:ncol(df.pm)]), 2, mean) #Find training data mean by column
pm_trainSD <- apply((df.pm[,7:ncol(df.pm)]), 2, sd) #Find training data SD by column
df.pm[,7:ncol(df.pm)] <- scale(df.pm[,7:ncol(df.pm)]) #Scale training explanatory data

df.sb <- data[data$SB == 1, colnames(data)[c(15:22, 27:ncol(data))]] #Selecting SB cultivation data.
colnames(df.sb)[1:6] <- c("CAS", "EXT", "INT", "RF", "DEC", "IRR")
sb_trainMean <- apply((df.sb[,9:ncol(df.sb)]), 2, mean) #Find training data mean by column
sb_trainSD <- apply((df.sb[,9:ncol(df.sb)]), 2, sd) #Find training data SD by column
df.sb[,9:ncol(df.sb)] <- scale(df.sb[,9:ncol(df.sb)]) #Scale training explanatory data

df.list <- list(df.cro, df.fm, df.pm, df.sb) %>%
  sapply(function (x) cbind(decostand(x[1:(which(colnames(x) == "LONG") - 1)], 
                                      method = "hellinger"),
                            sepcol = rep(NA, nrow(x)),
                            x[(which(colnames(x) == "LATI") + 1):ncol(x)]))
names(df.list) <- c("cro", "fm", "pm", "sb")
list2env(df.list, .GlobalEnv)

##Testing dataset (interviews)
test.i.dataset <- read.csv("2-interviews.csv")
data.i <- textshape::column_to_rownames(test.i.dataset, loc = 1)

df.cro.i <- data.i[, colnames(data.i)[c(4:6, 2:3, 27:ncol(data.i))]] #Selecting crop data.
colnames(df.cro.i)[1:3] <- c("FM", "PM", "SB")
for(i in 1:length(cro_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- cro_trainMean[i]
  SD <- cro_trainSD[i]
  df.cro.i[,(i+5)] <- (df.cro.i[,(i+5)] - M) / SD
}
df.cro.i.res <- data.frame(df.cro.i[,1:3])
df.cro.i.ind <- data.frame(df.cro.i[,6:ncol(df.cro.i)])

df.fm.i <- data.i[data.i$FM_cultiv == 1, colnames(data.i)[c(7:10, 2:3, 21:22, 27:ncol(data.i))]] #Selecting FM cultivation data.
colnames(df.fm.i)[1:4] <- c("EXT", "INT", "RF", "IRR")
for(i in 1:length(fm_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- fm_trainMean[i]
  SD <- fm_trainSD[i]
  df.fm.i[,(i+6)] <- (df.fm.i[,(i+6)] - M) / SD
}
df.fm.i.res <- data.frame(df.fm.i[,1:4])
df.fm.i.ind <- data.frame(df.fm.i[,7:ncol(df.fm.i)])

df.pm.i <- data.i[data.i$PM_cultiv == 1, colnames(data.i)[c(11:14, 2:3, 23:24, 27:ncol(data.i))]] #Selecting FM cultivation data.
colnames(df.pm.i)[1:4] <- c("EXT", "INT", "RF", "IRR")
for(i in 1:length(pm_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- pm_trainMean[i]
  SD <- pm_trainSD[i]
  df.pm.i[,(i+6)] <- (df.pm.i[,(i+6)] - M) / SD
}
df.pm.i.res <- data.frame(df.pm.i[,1:4])
df.pm.i.ind <- data.frame(df.pm.i[,7:ncol(df.pm.i)])

df.sb.i <- data.i[data.i$SB_cultiv == 1, colnames(data.i)[c(15:20, 2:3, 25:ncol(data.i))]] #Selecting SB cultivation data.
colnames(df.sb.i)[1:6] <- c("CAS", "EXT", "INT", "RF", "DEC", "IRR")
for(i in 1:length(sb_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- sb_trainMean[i]
  SD <- sb_trainSD[i]
  df.sb.i[,(i+8)] <- (df.sb.i[,(i+8)] - M) / SD
}
df.sb.i.res <- data.frame(df.sb.i[,1:6])
df.sb.i.ind <- data.frame(df.sb.i[,9:ncol(df.sb.i)])

##Testing dataset (aggregated data)
test.c.dataset <- read.csv("3-inteviews_aggregate.csv")
data.c <- textshape::column_to_rownames(test.c.dataset[,-1], loc = 1)

df.cro.c <- data.c[, colnames(data.c)[c(2:4, 25:ncol(data.c))]] #Selecting crop data.
colnames(df.cro.c)[1:3] <- c("FM", "PM", "SB")
for(i in 1:length(cro_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- cro_trainMean[i]
  SD <- cro_trainSD[i]
  df.cro.c[,(i+3)] <- (df.cro.c[,(i+3)] - M) / SD
}
df.cro.c.res <- data.frame(df.cro.c[,1:3])
df.cro.c.ind <- data.frame(df.cro.c[,4:ncol(df.cro.c)])

df.fm.c <- data.c[data.c$FM_cultiv == 1, colnames(data.c)[c(5:8, 19:20, 25:ncol(data.c))]] #Selecting FM cultivation data.
colnames(df.fm.c)[1:4] <- c("EXT", "INT", "RF", "IRR")
for(i in 1:length(fm_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- fm_trainMean[i]
  SD <- fm_trainSD[i]
  df.fm.c[,(i+4)] <- (df.fm.c[,(i+4)] - M) / SD
}
df.fm.c.res <- data.frame(df.fm.c[,1:4])
df.fm.c.ind <- data.frame(df.fm.c[,5:ncol(df.fm.c)])

df.pm.c <- data.c[data.c$PM_cultiv == 1, colnames(data.c)[c(9:12, 21:22, 25:ncol(data.c))]] #Selecting FM cultivation data.
colnames(df.pm.c)[1:4] <- c("EXT", "INT", "RF", "IRR")
for(i in 1:length(pm_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- pm_trainMean[i]
  SD <- pm_trainSD[i]
  df.pm.c[,(i+4)] <- (df.pm.c[,(i+4)] - M) / SD
}
df.pm.c.res <- data.frame(df.pm.c[,1:4])
df.pm.c.ind <- data.frame(df.pm.c[,5:ncol(df.pm.c)])

df.sb.c <- data.c[data.c$SB_cultiv == 1, colnames(data.c)[c(13:18, 23:ncol(data.c))]] #Selecting SB cultivation data.
colnames(df.sb.c)[1:6] <- c("CAS", "EXT", "INT", "RF", "DEC", "IRR")
for(i in 1:length(sb_trainMean)){ #Scaling using mean and SD from training dataset.
  M <- sb_trainMean[i]
  SD <- sb_trainSD[i]
  df.sb.c[,(i+6)] <- (df.sb.c[,(i+6)] - M) / SD
}
df.sb.c.res <- data.frame(df.sb.c[,1:6])
df.sb.c.ind <- data.frame(df.sb.c[,7:ncol(df.sb.c)])

#Model building
##Redundancy Analysis and Forward Selection

x = cro #Change "cro" by "fm", "pm" or "sb".

x.dep <- x[1:(which(colnames(x) == "sepcol") - 1)] #Subset response variables
x.ind <- x[(which(colnames(x) == "sepcol") + 1):ncol(x)] #Subset predictor variables
x.rda <- rda(x.dep ~ ., x.ind) ; print(x.rda) #Perform RDA
x.rda.r2 <- RsquareAdj(x.rda)$adj.r.squared ; print (x.rda.r2) #Adjusted R2
x.rda0 <- rda(x.dep ~ 1, x.ind)
x.model <- ordiR2step (x.rda0, x.rda)
(x.model.r2 <- RsquareAdj(x.model)$adj.r.squared) #Adjusted R2
anova.cca(x.model, step = 1000) #Test selected variables. If significant, draw triplots.
anova.cca(x.model, step = 1000, by = "axis") #Test RDA axes. Draw triplots using the significant ones.
plot(x.model, choices = 1:2, scaling = 1, main = "Triplot RDA") #Distance triplot using dep. column weighted sums.

##dbMEM analysis

x.xy = df.cro[,c("LONG", "LATI")] #Change "cro" by "fm", "pm" or "sb".

source('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
quickMEM(x.dep, x.xy)

##Variation partitioning with spatial components

x.model
x.xy.rda0 <- rda(x.dep ~ 1, x.xy)
x.xy.rda <- rda(x.dep ~ ., x.xy)
x.xy.model <- ordiR2step(x.xy.rda0, x.xy.rda)

#For crops use:
x.vp <- varpart(x.dep, ~ x.ind$BIO10_v + x.ind$GHI_m + x.ind$BIO15_m, ~ x.ind$WC11_m + x.ind$PH1_m + x.ind$CLAY2_v, x.xy)
#For FM use: 
#x.vp <- varpart(x.dep, ~ x.ind$PCI_m, ~ x.ind$TS2_m + x.ind$TP1_m)
#For PM use:
#x.vp <- varpart(x.dep,~ x.ind$PMcycle_m,~ x.ind$BIO4_v +x.ind$BIO8_m + x.ind$BIO9_m, ~ x.ind$WC21_v + x.ind$GRAV2_m + x.ind$CLAY1_m + x.ind$OC1_v + x.ind$SILT1_v,x.xy)
#For SB use:
#x.vp <- varpart(x.dep, ~ x.ind$SBcycle_m, ~ x.ind$CEC2_v + x.ind$CEC1_v + x.ind$SOC_m)

print(x.vp)
plot(x.vp)

#Model validation
## Model fitting assessment by Performance Measures.

x = df.cro #Change "cro" by "fm", "pm" or "sb".
x.res = df.cro[1:(which(colnames(df.cro) == "LONG") - 1)] #Change "df.cro" by "df.fm", "df.pm" or "df.sb".

x.model
x.pred <- fitted(x.model, type = "response")
x.thres <- data.frame((matrix(nrow = 1, ncol=0)))
for(i in 1:ncol(x.pred)) {
  thres_data <- data.frame(cbind(seq(1:(nrow(x.res))), x.res[,i], x.pred[,i])) 
  t <- optimal.thresholds(thres_data, opt.methods = "MaxSens+Spec")
  colnames(t)[names(t) == "X3"] <- colnames(x.res[i])
  x.thres <- cbind(x.thres, t[2])
  print(x.thres)
}

x.predicted <- as.data.frame(x.pred)
for (i in 1:ncol(x.predicted)) {
  x.predicted[i] <- ifelse(x.predicted[i] < x.thres[,i], 0, 1) 
}

cf <- table(factor(as.matrix(x.res)), factor(as.matrix(x.predicted))) # Confusion matrix
TP <- cf[2,2] # True positives.
TN <- cf[1,1] # True negatives.
FP <- cf[1,2] # False positives.
FN <- cf[2,1] # False negatives.
acc <- (TP+TN)/(TP+FP+FN+TN) # Accuracy = (TP+TN)/(TP+FP+FN+TN)
prec <- TP/(TP+FP) # Precision = TP/(TP+FP)
rec <- TP/(TP+FN) # Recall = TP/(TP+FN)
F1_score <- 2*((rec*prec)/(rec+prec)) # F1 Score = 2*((Recall*Precision)/(Recall+Precision))
r <- cbind(acc, prec, rec, F1_score)
colnames(r) <- c("Acc.(All)", "Prec.(TP Rate)", "Rec.(% all TP)", "F1 (Class. Strength)")                  
x.fitrep <- r

print(x.fitrep)


##Cross validation against ethnographic data.

x.model
x.pred <- fitted(x.model, type = "response")
x.thres <- data.frame((matrix(nrow = 1, ncol=0)))
for(i in 1:ncol(x.pred)) {
  thres_data <- data.frame(cbind(seq(1:(nrow(x.res))), x.res[,i], x.pred[,i])) 
  t <- optimal.thresholds(thres_data, opt.methods = "MaxPCC")
  colnames(t)[names(t) == "X3"] <- colnames(x.res[i])
  x.thres <- cbind(x.thres, t[2])
  print(x.thres)
}

i.ind <- df.cro.i.ind #Change "df.cro" by "df.fm", "df.pm" or "df.sb".
i.res <- df.cro.i.res #Change "df.cro" by "df.fm", "df.pm" or "df.sb".
c.ind <- df.cro.c.ind #Change "df.cro" by "df.fm", "df.pm" or "df.sb".
c.res <- df.cro.c.res #Change "df.cro" by "df.fm", "df.pm" or "df.sb".

i.predicted <- as.data.frame(predict(x.model, i.ind, type = "response"))
c.predicted <- as.data.frame(predict(x.model, c.ind, type = "response"))
for (i in 1:ncol(i.predicted)) {
  i.predicted[i] <- ifelse(i.predicted[i] < x.thres[,i], 0, 1)
  c.predicted[i] <- ifelse(c.predicted[i] < x.thres[,i], 0, 1)
}

list.pred <- list(i.predicted, c.predicted)
list.real <- list(i.res, c.res)
x.report <- data.frame()
for (i in 1:length(list.pred)) {
  cf <- table(factor(as.matrix(list.real[[i]])), factor(as.matrix(list.pred[[i]]))) # Confusion matrix
  TP <- cf[2,2] # True positives.
  TN <- cf[1,1] # True negatives.
  FP <- cf[1,2] # False positives.
  FN <- cf[2,1] # False negatives.
  acc <- (TP+TN)/(TP+FP+FN+TN) # Accuracy = (TP+TN)/(TP+FP+FN+TN)
  prec <- TP/(TP+FP) # Precision = TP/(TP+FP)
  rec <- TP/(TP+FN) # Recall = TP/(TP+FN)
  F1_score <- 2*((rec*prec)/(rec+prec)) # F1 Score = 2*((Recall*Precision)/(Recall+Precision))
  r <- cbind(acc, prec, rec, F1_score)
  colnames(r) <- c("Acc.(All)", "Prec.(TP Rate)", "Rec.(% all TP)", "F1 (Class. Strength)")                  
  x.report <- rbind(x.report, r)
}

rownames(x.report) <- c("x.i", "x.c")
print(x.report)













