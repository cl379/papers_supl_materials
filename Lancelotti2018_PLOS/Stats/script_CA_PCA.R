### install packages -------------------------------------------------------------------

install.packages('ca')
install.packages('ggplot2')
install.packages('FactoMineR')
install.packages("devtools")
install.packages('gridExtra')
devtools::install_github("slowkow/ggrepel")
devtools::install_github("ggbiplot", "vqv")
    
### load packages ----------------------------------------------------------------------------

require(ca)
require(ggplot2)
require(FactoMineR)
library(ggrepel)
library(ggbiplot)
library(gridExtra)

### CA of combined charcoal and phytoliths data for all sites---------------------------------

### load data for CA
tb <- read.csv('phyto_charc_round.csv', row.names = 1)

### simple ca
my.ca <- ca(tb)

### CA plot with ggplot
### first generate a new table
TB <- as.matrix(tb)
df <- data.frame()
for (i in 1:nrow(TB)) {
    for (j in 1:ncol(TB)){
        if (TB[i,j] > 0) {
            for (n in 1:TB[i,j]) {
                df <- rbind(df, data.frame(var1= rownames(TB)[i], var2=colnames(TB)[j]))
            }
        }
    }
}

# manual ca build (prepare data)
cats <- apply(df, 2, function(x) nlevels(as.factor(x)))                     # number of factors
mca.ob <- MCA(df, graph = F)                                                # correspondence
gg.aux <- data.frame(mca.ob$var$coord, Variable = rep(names(cats), cats))   # variable coordinates
mcaObs = data.frame(mca.ob$ind$coord)

# ggplot CA 
ca1 <- ggplot(data = gg.aux, aes(x = Dim.1, y = Dim.2, label = rownames(gg.aux))) +
    geom_hline(yintercept = 0, colour = 'gray70') +
    geom_vline(xintercept = 0, colour = 'gray70') +
    geom_text_repel(aes(colour = Variable, fontface='bold'), size=5, nudge_x = 0.2) +  # added nudge_x = 0.2
    scale_colour_discrete(l=40) +
    geom_point(aes(colour = Variable, shape = Variable), size = 3) +
    theme(legend.position="none")
ca1

### PCA----------------------------------------------------------------------------------------

### load data for PCA
tb.h <- read.csv('chem_harappa.csv', row.names=1) # Data Harappa
tb.k <- read.csv('chem_kanmer.csv', row.names=1) # Data Kanmer
tb.s <- read.csv('chem_shikarpur.csv', row.names=1) # Data Shikarpur
tb.a <- read.csv('chem_alamgirpur.csv', row.names=1) # Data Alamgirpur

### generate a new table and a object for grouping (one for each site)
## Harappa
tb.pca.h<-prcomp(tb.h[, -1], center=TRUE, scale.=TRUE)
tb.context.h<-tb.h[, 1]

## Kanmer
tb.pca.k<-prcomp(tb.k[, -1], center=TRUE, scale.=TRUE)
tb.context.k<-tb.k[, 1]

#Shikarpur
tb.pca.s<-prcomp(tb.s[, -1], center=TRUE, scale.=TRUE)
tb.context.s<-tb.s[, 1]

## Alamgirpur
tb.pca.a<-prcomp(tb.a[, -1], center=TRUE, scale.=TRUE)
tb.context.a<-tb.a[, 1]

### Analyse PCA results (subsitute 'tb.pca' with the corersponding site, e.g. 'tb.pca.h' for Harappa)
print(tb.pca)                   # Returns standard deviation for each PC and their rotation
plot(tb.pca, type = "lines")    # Returns a plot of the variances (y-axis) associated with the PCs (x-axis). Types are "lines" or "barplot"
summary(tb.pca)                 # Describe the importance of the PCs

### Create common color palette for all graphs
group.colors<- c('Ash accumulation' = 'maroon2', 'Dung ash' = 'grey34', 'Dung fresh' = 'orange2', Floor = 'blueviolet', Fireplace = 'Blue4', Pit = 'Red1', Plaster = 'lightseagreen', Street = 'darkgreen')

### CReate PCA biplots for each site
## Harappa
g.h <- ggbiplot(tb.pca.h, obs.scale = 1, var.scale = 1, groups = tb.context.h, ellipse = TRUE, circle = TRUE)
g.h <- g.h + scale_color_manual(values = group.colors) + theme(legend.position="none")  + xlim(-5, 6.5) + ylim(-4.5, 4)

## Kanmer
g.k <- ggbiplot(tb.pca.k, obs.scale = 1, var.scale = 1, groups = tb.context.k, ellipse = TRUE, circle = TRUE)
g.k <- g.k + scale_color_manual(values = group.colors) + theme(legend.position="none")  + xlim(-5, 6.5) + ylim(-4.5, 4)

## Shikrpur
g.s <- ggbiplot(tb.pca.s, obs.scale = 1, var.scale = 1, groups = tb.context.s, ellipse = TRUE, circle = TRUE)
g.s <- g.s + scale_color_manual(values = group.colors) + theme(legend.position="none")  + xlim(-5, 6.5) + ylim(-4.5, 4)

## Alamgirpur
g.a <- ggbiplot(tb.pca.a, obs.scale = 1, var.scale = 1, groups = tb.context.a, ellipse = TRUE, circle = TRUE)
g.a <- g.a + scale_color_manual(values = group.colors) + theme(legend.position="none") + xlim(-5, 6.5) + ylim(-4.5, 4)

### Define groups' colours
Colors<- c('Ash accumulation', 'Dung ash', 'Dung fresh', 'Floor', 'Fireplace', 'Pit', 'Plaster', 'Street')

### CReate common legend for the 4 grobs
aux = data.frame(x=rnorm(8), y=rnorm(8), col=Colors)
g.aux <- ggplot(aux, aes(x=x, y=y, colour=Colors)) + geom_point() + scale_color_manual(values = group.colors) + 
  theme(legend.title=element_blank()) + theme(legend.direction = 'horizontal', legend.position = 'bottom')

tmp <- ggplot_gtable(ggplot_build(g.aux))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

### Plot the grobs in the same space with a common legend
pca.chem <- grid.arrange(arrangeGrob(g.h, g.k, g.s, g.a, ncol=2),
                   legend, nrow=2, heights=c(10, 1))


