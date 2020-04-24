# Statistics (ANOVA and PCA) on phytolith analysis from Els Trocs

### Install necessary packages and load libries
install.packages("FactoMineR")
install.packages("car")
library(car)
library(FactoMineR)

# Read the txt file containing the data and examine the data
etr=read.table("ETR_fum_stat.txt", header=TRUE, row.names=2)
summary(etr)

### Create a model for ANOVA based on the data
aov.model<-(lm(cbind(Concentration, Morphotypes, Unidentified, Inflorescence, Leaf.culm, Woody, Spherulites)~SU, data=etr))

# lm: linear model
# cbind: the variables you want to include in your analysis as they are called in the data file
# ~ indicates the factor
# data=: your data file

# Display results of the analysis
summary(aov.model)

# Perform Post-hoc tests Bonferroni
outlierTest(aov.model, cutoff=0.05)

# Boxplots of catergories
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 2, 4, byrow=TRUE))
boxplot(etr$Concentration~etr$SU, main= "Phytolith concentration", lab="10, 20, 53, 29")
boxplot(etr$Morphotypes~etr$SU, main= "Number of morphotypes", lab="10, 20, 53, 29")
boxplot(etr$Inflorescence~etr$SU, main= "Inflorescence", lab="10, 20, 53, 29")
boxplot(etr$Leaf.culm~etr$SU, main= "Leaf/culm", lab="10, 20, 53, 29")
boxplot(etr$Woody~etr$SU, main= "Woody dicotyledons", lab="10, 20, 53, 29")
boxplot(etr$Unidentified~etr$SU, main= "Unidentified", lab="10, 20, 53, 29")
boxplot(etr$Spherulites~etr$SU, main= "Spherulites", lab="10, 20, 53, 29")

### Perform PCA on the data
res.pca = PCA(etr[,2:8], scale.unit=TRUE, ncp=5, graph=FALSE)
#etr: the dataset used
#[,2:8]: the column analysed (scaled data)
# scale.unit: chose wether to scale the data or not (if not scaled a covariance matrix instead that a corrlation matrix is used to perform PCA)
#ncp: number of dimensions retained in the analysis
#graph: chose whether to display the plots (if you chose not to display the plot at this stage you canalways display them later using plotPCA() - see below)

# Check results
summary(res.pca)

#Display scores
res.pca$ind$coord

# Plot PCa graphs when graphs = FALSE is chosen during analysis
layout(matrix(c(1,2),1,2, byrow=TRUE))
# This option prints the samples
plot.PCA(res.pca, axes(1,2), choix="ind", col.ind="black", label=c("ind"), title="Samples")
# This option prints the variables
plot.PCA(res.pca, axes(1,2), choix="var", col.var="black", label=c("var"), title="Variables")


