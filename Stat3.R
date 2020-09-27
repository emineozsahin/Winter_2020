##=========================================================*
#**********************************************************#
## By: Emine Ozsahin 
#**********************************************************#
##=========================================================*
#loading the libraries ----

library(ggplot2)
library(GGally)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(tidyverse)

#setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/BINF6970_Stat")

##=========================================================*
#**********************************************************#

#Problem 1 ----

#Consider the data on bulls in data_1.csv. 
#Utilizing the seven variables (Yr Hgt, FtFrBody, PrctFFB, Frame, BkFat, SaleHt, and Sale Wt) 
#to perform a principal component analysis using the covariance matrix S and the correlation matrix R. 

# Variables are defined as follows:
# Col. 1: Breed (1 = Angus, 5 = Hereford, 8 = Simental)
# Col. 2: SalePr = selling price ($)
# Col. 3: YrHgt = yearling height at shoulder (in)
# Col. 4: FtFrBody = fat free body wt (lbs)
# Col. 5: PrctFFB = percent fat-free body weight (%)
# Col. 6: Frame = size scale (1 = small to 8 = large)
# Col. 7: BkFat = back fat (in)
# Col. 8: SaleHt = sale height at shoulder (in)
# Col. 9: SaleWt = sale weight (lbs)

#**********************************************************#
##=========================================================*

#Loading in the data as a matrix.
bulls <- as.matrix(read.csv("Data_1.csv"))
head(bulls) 

#Remove the variables that are not indicated in the question.
bulls2 <- bulls[,c(-1, -2)]

#**********************************************************#
# Problem 1 (a) ----
# Determine the appropriate number of components to effectively summarize the sample variability. 
# Construct a scree plot to aid your determination.

#First, the covariance (S) and correlation (R) matrices are generated. 
S <- var(bulls2) 
R <- cor(bulls2)

#Eigen vectors and values for correlation and covariance matrices are found.
S_eigen <- eigen(S)
R_eigen <- eigen(R)

#The results are summarized in a matrix. The first 7 rows of the matrix will contain the eigen vectors and 
#the last 2 rows will contain the eigen values and cumulative percentage of each PC respectively.
resultmat_var <- matrix(NA, ncol=7, nrow = 9)
resultmat_cor <- matrix(NA, ncol=7, nrow = 9)

#Row names are added.
rownames(resultmat_var) <- c(colnames(bulls2)  , "lambda.hat"   , "cumulative %" )
rownames(resultmat_cor) <- c(colnames(bulls2)  , "lambda.hat"   , "cumulative %" )

#Column names are added.
colnames(resultmat_var) <- c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7")
colnames(resultmat_cor) <- c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7")

#Eigen values are added to the matrix under the column lambda.hat
resultmat_var[8,] <- S_eigen$values
resultmat_cor[8,] <- R_eigen$values

#Cumulative sums of PC variances are added as well.
resultmat_var[9,] <- cumsum(S_eigen$values/sum(S_eigen$values)*100)
resultmat_cor[9,] <- cumsum(R_eigen$values/sum(R_eigen$values)*100)

#Finally, the eigen vectors are added.
resultmat_var[1:7,1:7] <- S_eigen$vectors
resultmat_cor[1:7,1:7] <- R_eigen$vectors

#The result summary matrices for both correlation and variance are complete.There are 7 Principal Components.
round( resultmat_var, 3)
round( resultmat_cor, 3)

# we used the following line of code to obtain the results as a table to represent in the report.
#write.csv(round(resultmat_cor, 3),'Resultmat.csv') 
#write.csv(round(resultmat_var, 3),'Resultmat_cov.csv')

### The Scree Plot using the covariance matrix (non-scaled data).
#png("scree_bulls_S_matrix.png", units="in", width=5, height=5, res=300)

plot(as.ts(S_eigen$values) , ylab="Eigen values" , xlab="PC" , main = "Scree Plot using Covariance Matrix", col='brown',lwd=2)

#dev.off()

### The Scree PLOT using the correlation matrix (scaled data).
#png("scree_bulls_R_matrix.png", units="in", width=5, height=5, res=300)

plot(as.ts (R_eigen$values) , ylab="Eigen values" , xlab="Number of Principle Components" , main = "Scree Plot using Correlation Matrix", col='brown',lwd=2)

#dev.off()

#Calculation of the Principal components using princomp() function.

#The cor=FALSE option uses the covariance matrix (S) and the results represent the non-scaled data.
bulls.pca <- princomp(bulls2, cor = FALSE, scores = TRUE) 

#The cor=TRUE option uses the correlation matrix (R) and the results represent the scaled the data.
bulls.pca_sc <- princomp(bulls2, cor = TRUE, scores = TRUE) 

# when  using the covariance matrix (non-scaled), check the eigen values which is same with 
# the principle components computed with non-scaled data. 
eigen(S)[1]

bulls.pca$sdev^2 #Similar with eigen values of the covariance matrix.

# when using the correlation matrix (scaled) check the values which is same with the 
# principle components computed with scaled data.
eigen(R)[1]

bulls.pca_sc$sdev^2 # Identical to eigen values of the correlation matrix.

# Scree plots generated by eigen values calculated from the covariance and correlation matrices 
# are representative of non-scaled and scaled data respectively. 

summary(bulls.pca) # NON-SCALED: FIRST TWO PCs EXPLAIN MOST OF THE VARIATION.
resultmat_var[9,]

summary(bulls.pca_sc) # SCALED: FIRST FOUR PCs EXPLAIN MOST OF THE VARIATION.
resultmat_cor[9,]

#**********************************************************#
# Problem 1 (b) Interpret the sample principal components ----

#**********************************************************#
# Problem 1 (c) ----
# Do you think it is possible to develop a "body size" or 
# "body configuration" index from the data on the seven variables above? Explain.

#Correlation between all the variables and frame. 
cor(bulls, bulls[,6])

# Size Groups 
unique(bulls[,6])
length(which(bulls[,6]==5)) #16
length(which(bulls[,6]==6)) #28
length(which(bulls[,6]==7)) #24
length(which(bulls[,6]==8)) #8
# Col. 6: Frame = size scale (1 = small to 8 = large)
bulls.size <- c(rep("5", 16), rep("6", 28), rep("7", 24), rep("8", 8))

# A biplot was generated by using first two PCs with scaled data and grouped by frame.
#png("ggbiplot_bulls_frame.png", units="in", width=5, height=5, res=300)



ggbiplot(bulls.pca_sc, ellipse=TRUE, obs.scale = 1, var.scale = 1, labels=bulls[,6], choices=c(1,2), groups=bulls.size) +
  scale_colour_manual(name="Body Size", values= c("forest green", "red4", "dark blue", "orange")) +
  theme_minimal()+
  theme(legend.direction = 'horizontal', legend.position = 'top') 

#dev.off()



#**********************************************************#
# Problem 1 (d) Plot the principle components ----
# Using the values for the first two principal components, 
# plot the data in a two-dimensional space with PC1 along the 
# vertical axis and PC2 along the horizontal axis. 
# Can you distinguish groups representing the three breeds of cattle? 
# Are there any outliers?

# Breed Groups
length(which(bulls[,1]==1)) #32
length(which(bulls[,1]==5)) #17
length(which(bulls[,1]==8)) #27
# Breed (1 = Angus, 5 = Hereford, 8 = Simental)
bulls.breed <- c(rep("Angus", 32), rep("Hereford",17), rep("Simental", 27))

# biplot was generated by using first two PCs with scaled data and grouped by breed
#png("ggbiplot_bulls_breed.png", units="in", width=5, height=5, res=300)

ggbiplot(bulls.pca_sc, ellipse=TRUE, obs.scale = 1, var.scale = 1, labels=bulls[,1], groups=bulls.breed) +
  scale_colour_manual(name="Breed", values= c("forest green", "red4","dark blue")) +
  theme_minimal() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  coord_flip()

#dev.off()

# biplot was generated by using first two PCs with non-scaled data and grouped by breed
#png("ggbiplot_notscaledbulls.png", units="in", width=5, height=5, res=300)

ggbiplot(bulls.pca, ellipse=TRUE, obs.scale = 1, var.scale = 1, labels=bulls[,1], groups=bulls.breed) +
  scale_colour_manual(name="Breed", values= c("forest green", "red4","dark blue")) +
  theme_minimal()+
  theme(legend.direction = 'horizontal', legend.position = 'top') + coord_flip()

#dev.off()

# plot using covariation and correlation matrix which is not neserrary for report

ggbiplot(princomp(S, cor = FALSE)) # SAME AS NON-SCALED DATA

ggbiplot(princomp(R, cor = FALSE)) # SAME AS SCALED DATA

#**********************************************************#
# Problem 1 (e) QQ-Plot ----

# Construct a Q-Q plot using the first principal component (scaled, correlation matrix).

ysc <- bulls.pca_sc$scores[,1]

quantile(ysc)

par(mfrow=c(1,2))
#png("quantile_bulls_sc.png", units="in", width=5, height=5, res=300)
qqnorm(ysc); qqline(ysc, col = "steelblue", lwd = 2)
#dev.off()

#png("hist_bulls_sc.png", units="in", width=5, height=5, res=300)
hist(ysc, main = "Histogram of quantiles", xlab= "Quantiles")
#dev.off()
par(mfrow=c(1,1))


# Construct a Q-Q plot using the first principal component (not scaled, variance matrix).

y <- bulls.pca$scores[,1]

quantile(y)

par(mfrow=c(1,2))
#png("quantile_bulls.png", units="in", width=5, height=5, res=300)
qqnorm(y); qqline(y, col = "steelblue", lwd = 2)
#dev.off()

#png("hist_bulls.png", units="in", width=5, height=5, res=300)
hist(y, main = "Histogram of quantiles", xlab= "Quantiles")
#dev.off()
par(mfrow=c(1,1))

##=========================================================*
#**********************************************************#

#Problem 2 ----

#The data on national track records for women are provided (women.csv).

#**********************************************************#
##=========================================================*

#Load in the data.
women <- read.csv("women.csv")
str(women)

#**********************************************************#
# Problem 2 (a) correlation matrix ----
# Obtain the sample correlation matrix R for these data 
# and determine its eigenvalue and eigenvectors.

#Correlation matrix is created excluding the countries column and eigen values and vectors are found.

R_women <- cor(women[-1])

#write.csv(round(R_women, 3),'R_women.csv')

women_value <- eigen(R_women)$values

women_vector <- eigen(R_women)$vectors

#**********************************************************#
# Problem 2 (b) first two principal components ----
# Determine the first two principal components for the standardized variables. 
# Prepare a table showing the principal components, and their individual and 
# cumulative percentage contribution to the total (standardized) sample variance. 
# Comment on the results.

# Principle components by princomp() function.
# Correlation Matrix is equavalent to the scaled variables in the dataset. 
# Choosing the option cor = TRUE in the princomp function allows us to use the standardized variables. 
# Therefore, we did not use the scale() function to scale the data.

women.pca <- princomp(women[,-1], cor = TRUE, scores = TRUE) 

women.pca$loadings [,1:2] # First two principle components

summary(women.pca) # cumulative variance could be obtained

ggbiplot(women.pca) # visualization 

## Results matrix.

resultmat_women_cor <- matrix(NA, ncol=7, nrow = 9)

rownames(resultmat_women_cor) <- c( colnames(women[,-1]), "lambda.hat", "cumulative %" )

colnames(resultmat_women_cor) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7")

## Eigen values added.
resultmat_women_cor[8,] <- women_value

## cumulative sum of PC variances added.
resultmat_women_cor[9,] <- cumsum(women_value/sum(women_value)*100)

# Eigen vectors added.
resultmat_women_cor[1:7,1:7] <- women_vector

#Result summary matrix complete.
round(resultmat_women_cor, 3)

#write.csv(round(resultmat_women_cor, 3),'Resultmat_women.csv')

#**********************************************************#
# Problem 2 (c) ----
# Interpret the two principal components obtained in Part b. 
# (Note that the first component is essentially a normalized unit vector and might measure the athletic excellence of a given nation. 
# The second component might measure the relative strength of a nation at the various running distances)

women.pca$loadings[,1:2]
#write.csv(round(t(women.pca$loadings[,1:2]), 3),'PC1_PC2_women.csv')

str(women)
unique(women[,1]) # There are 54 countries that are unique therefore there are no groups (we may group the countries by continent)

# Biplot was generated using the first two PCs found with the scaled data.
#png("ggbiplot_women.png", units="in", width=5, height=5, res=300)

ggbiplot(women.pca, labels=women[,1], choices=c(1, 2)) + 
  theme_minimal()

#dev.off()

#**********************************************************#
# Problem 2 (d) ----
# Rank the nations based on their score on the first principal component. 
# Does this ranking correspond with your intuitive notion of athletic excellence for the various countries? 
# List top ten ranked nations and the bottom five ranked nations 
# (try to visualize the information instead of reporting as tabular summary).

# A dataframe containing the countries and their respective scores is created. 
# These scores are multiplied by -1.

scores_women <- data.frame(country = women[,1], score = women.pca$scores[,1]*-1)
scores_sorted_women <- scores_women[order(-scores_women$score),] 
scores_sorted_women <- cbind (scores_sorted_women, rank = 1:54)
scores_sorted_women

# The top ten and last five countries are extracted.
topten_lastfive <- rbind(scores_sorted_women[1:10,], scores_sorted_women[50:54,])
topten_lastfive

# Visualization of the rankings.

#png("barplot_women.png", units="in", width=5, height=5, res=300)

barplot(topten_lastfive[order(topten_lastfive[,3],decreasing = TRUE),][,2],
        space=c(0,00),
        legend.text=TRUE,
        beside=TRUE,
        main = "Athletic Excellence of Women",
        xlab = "PC1 Score",
        ylab = "Country",
        names.arg = topten_lastfive[order(topten_lastfive[,3],decreasing = TRUE),][,1],
        col = "grey",
        horiz = TRUE,
        density=NA,
        axes=TRUE, 
        cex.names=0.6, 
        las=1)

#dev.off()

##=========================================================*
#**********************************************************#

# Problem 3 ----

#The data on national track records for men are provided (men.csv).
#Repeat the principal component analysis outlined Problem 2 for men’s dataset. 
#Are the results consistent with those obtained from the women's data? 
#Page Limit for Problems 2-4: Four pages

#**********************************************************#
##=========================================================*

#Load in the data.
men <- read.csv("men.csv")
str(men)

#**********************************************************#
# Problem 3 (a) correlation matrix ----
# Obtain the sample correlation matrix R for these data 
# and determine its eigenvalue and eigenvectors.

#Correlation matrix is created excluding the countries column. 
R_men <- cor(men[-9])

# Eigen values and vectors of the correlation matrix were found.
men_value <- eigen(R_men)$values

men_vector <- eigen(R_men)$vectors

#**********************************************************#
# Problem 3 (b) Determine the first two principal components for the standardized variables. 
# Prepare a table showing the principal components, and their individual and 
# cumulative percentage contribution to the total (standardized) sample variance. 
# Comment on the results.

# The results matrix

resultmat_men_cor <- matrix(NA, ncol=8, nrow = 10)

rownames(resultmat_men_cor) <- c(colnames(men[,-1]), "lambda.hat", "cumulative %")

colnames(resultmat_men_cor) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")

## Eigen values added.
resultmat_men_cor[9,] <- men_value

## Cumulative sum of PC variances added.
resultmat_men_cor[10,] <- cumsum(men_value/sum(men_value)*100)

#Eigen vectors added.
resultmat_men_cor[1:8,1:8] <- men_vector

#Results matrix complete.
round(resultmat_men_cor, 3)

#write.csv(round(resultmat_men_cor, 3),'R_men.csv') 

#Principle components found using princomp function.

men.pca <- princomp(men[,-9], cor = TRUE, scores = TRUE)

#**********************************************************#
#(c) Interpret the two principal components obtained in Part b. 
# (Note that the first component is essentially a normalized unit vector and might measure the athletic excellence of a given nation. 
# The second component might measure the relative strength of a nation at the various running distances)

# First two principle components
men.pca$loadings[,1:2]

#write.csv(round(t(men.pca$loadings[,1:2]), 3),'PC1_PC2_men.csv')
summary(men.pca)

#Observe the dataset.
str(men)  

unique(men[,9]) #There are 55 countries (all unique), therefore there are no groups.

# Biplot was generated using the first two PCs found with the scaled data.
#png("ggbiplot_men.png", units="in", width=5, height=5, res=300)

ggbiplot(men.pca, labels=men[,9], choices=c(1, 2)) + 
  theme_minimal()

#dev.off()

#**********************************************************#
#(d) Rank the nations based on their score on the first principal component. Does this ranking correspond with your intuitive notion of athletic excellence for the various countries? List top ten ranked nations and the bottom five ranked nations (try to visualize the information instead of reporting as tabular summary).

# A dataframe containing the countries and their respective scores is created. 
# These scores are multipled by -1 and ranked.

scores_men <- data.frame(country = men[,9], score = men.pca$scores[,1]*-1)
scores_sorted_men <- scores_men[order(-scores_men$score),] 
scores_sorted_men <- cbind (scores_sorted_men, rank = 1:55)

#The top ten and last five countries are extracted.

topten_lastfive_men <- rbind(scores_sorted_men[1:10,], scores_sorted_men[51:55,])

# Country names are changed to the abreviated form
topten_lastfive_men$country <- as.character(topten_lastfive_men$country)

read.csv("Country_names.csv")

topten_lastfive_men$country <- topten_lastfive_men$country %>%
  str_replace_all(c("usa" = "USA", "gbni" = "GBR", "italy" = "ITA",  "ussr" = "RUS",  "gdr" = "East GER",  "frg" = "West GER",  
                    "australi" = "AUS", "france" = "FRA",  "kenya" = "KEN", "belgium" = "BEL", "singapor" = "SIN",
                    "png" = "PNG", "mauritiu" = "MRI", "wsamoa" = "SAM", "cookis" = "COK"))

topten_lastfive_men

#Visualization of the rankings.

#png("barplot_men.png", units="in", width=5, height=5, res=300)

barplot(topten_lastfive_men[order(topten_lastfive_men[,3],decreasing = TRUE),][,2],
        space=c(0,0),
        legend.text=TRUE,
        beside=TRUE,
        main = "Athletic Excellence of Men",
        xlab = "PC1 Score",
        ylab = "Country",
        names.arg = topten_lastfive_men[order(topten_lastfive_men[,3],decreasing = TRUE),][,1],
        col = "grey",
        horiz = TRUE,
        density=NA,
        axes=TRUE, 
        cex.names=0.6, 
        las=1)


barplot(GO$V3,
        col=c(rep("red",6), rep("blue",7), rep("green", 3)),
        space=c(0,0),
        legend.text=TRUE,
        beside=TRUE,
        main = "ene Ontology",
        xlab = "Number of Genes",
        ylab = "",
        names.arg = GO$V2,
        col = "grey",
        horiz = TRUE,
        density=NA,
        axes=TRUE, 
        cex.names=0.6, 
        las=1)

library(ggpubr)
ggbarplot(GO, x = "GO", y = "Number of Genes",
          fill = "GO class",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          
)
#dev.off()


## Extra ----

men_value #<- eigen(R_men)$values
diag(var(men.pca$scores))	

### The Scree PLOT

#png("scree_women_cor_matrix.png", units="in", width=5, height=5, res=300)

plot(as.ts(women_value) , ylab="variance" , xlab="PC" , main = " Scree Plot using Correlation Matrix", col='brown',lwd=2)

#dev.off()


### The Scree PLOT using correlation matrix (scaled data) for men.

#png("scree_men_cor_matrix.png", units="in", width=5, height=5, res=300)

plot(as.ts(men_value) , ylab="variance" , xlab="PC" , main = " Scree Plot using Correlation Matrix", col='brown',lwd=2)

#dev.off()

### Country name 
country <- read.csv("Country_names.csv")
length(country[,1]) #56 # last row is NA
length(men[,9]) #55 
dim(women) #54 8
cbind(country[-56,], country_men= men[,9], country_women=rbind(women, c(1:8))[,1])
# country name order does not match starting from 7th row in the files !!!!

