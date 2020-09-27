##=========================================================*
#**********************************************************#

## BINF 6970 Assignment 4
## By: Emine Ozsahin (1104289), 
##     Pasha Talebi Charmchi (1100889), and 
##     Rami Baghdan (1116137)
## Due Friday, April 10, 11:59 pm

#**********************************************************#
##=========================================================*
# Genes and their regions
# ZBTB38: 3q23 --> 3:141324186-141449792
# HMGA2 : 12q14 --> 12:65824460-65966291
# CYP19 : 15q21.2 --> 15:51208057-51338596

# Populations, Regions and Super Populations
# CHB --> CHINA --> EAS 
# FIN --> FINLAND --> EUR
# ESN --> NIGERIA --> AFR
# MXL --> MEXICAN --> AMR
# PJL --> PAKISTAN --> SAS
#**********************************************************#
##=========================================================*
# load package

setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/BINF6970_Stat")

# install_github("vqv/ggbiplot")
# install.packages("clues")
# install.packages("clustertend")
# BiocManager::install("snpStats")
library(snpStats)
library(factoextra)
library(ggbiplot)
library(clustertend)
library(clues)
library(dendextend)
library(colorspace)
library(cluster)
library(tidyverse)

##*********************************************
## Data Preprocessing ----
##*********************************************

# Load the vcf files
vcf_ZBTB38 <-read.table("3.141324186-141449792.ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", stringsAsFactors = FALSE)
vcf_HMGA2 <-read.table("12.65824460-65966291.ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", stringsAsFactors = FALSE)
vcf_CYP19 <-read.table("15.51208057-51338596.ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", stringsAsFactors = FALSE)

#Check the dimensions
dim(vcf_ZBTB38) #3899  470
dim(vcf_HMGA2) #3997  470
dim(vcf_CYP19) #3376  470

# Load in the names and get the ids for individuals
vcf <-readLines("3.141324186-141449792.ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf")

# colnames
vcf <- vcf[-(grep("#CHROM", vcf)+ 1):-(length(vcf))]
ids <- unlist(strsplit(vcf[length(vcf)],"\t"))
names(vcf_ZBTB38) <- ids
names(vcf_HMGA2) <- ids
names(vcf_CYP19) <- ids

plot(vcf_ZBTB38$QUAL)
plot(vcf_HMGA2$QUAL)
plot(vcf_CYP19$QUAL) # All qualities are 100

vcf_ZBTB38$INFO

# Extract the Allelle frequency
AF_ZBTB38 <- unlist(strsplit(vcf_ZBTB38$INFO,";"))
AF_HMGA2 <- unlist(strsplit(vcf_HMGA2$INFO,";"))
AF_CYP19 <- unlist(strsplit(vcf_CYP19$INFO,";"))

head(AF_ZBTB38)

AF_ZBTB38 <- data.frame(freq=gsub("(^AF=)", "", AF_ZBTB38[grep("^AF=", AF_ZBTB38)]))
AF_HMGA2 <- data.frame(freq=gsub("(^AF=)", "", AF_HMGA2[grep("^AF=", AF_HMGA2)]))
AF_CYP19 <- data.frame(freq=gsub("(^AF=)", "", AF_CYP19[grep("^AF=", AF_CYP19)]))


head(AF_ZBTB38)

hist(as.numeric(AF_ZBTB38$freq), main = "Allele frequency spectrum of the gene ZBTB38", 
     xlab="Frequency in Population", 
     ylab = "Proportion of Variants", 
     col = "grey")

hist(as.numeric(AF_HMGA2$freq), main = "Allele frequency spectrum of the gene HMGA2", 
     xlab="Frequency in Population", 
     ylab = "Proportion of Variants", 
     col = "grey")

hist(as.numeric(AF_CYP19$freq), main = "Allele frequency spectrum of the gene CYP19", 
     xlab="Frequency in Population", 
     ylab = "Proportion of Variants", 
     col = "grey")

# Filter the Alelle frequency less than 0.001 and and more than 0.999
vcf_ZBTB38 <- vcf_ZBTB38[which(as.character(AF_ZBTB38[,1]) > 0.001  &  as.character(AF_ZBTB38[,1]) < 0.999),]
vcf_HMGA2 <- vcf_HMGA2[which(as.character(AF_HMGA2[,1]) > 0.001  &  as.character(AF_HMGA2[,1]) < 0.999),] #0.01 0.9
vcf_CYP19 <- vcf_CYP19[which(as.character(AF_CYP19[,1]) > 0.001  &  as.character(AF_CYP19[,1]) < 0.999),]

# transpose the matrix only for variant data 
Variant_ZBTB38 <- t(vcf_ZBTB38[,-c(1:9)])
Variant_HMGA2 <- t(vcf_HMGA2[,-c(1:9)])
Variant_CYP19 <- t(vcf_CYP19[,-c(1:9)])

# Combine the variant data
Variant <- cbind(Variant_ZBTB38, Variant_HMGA2, Variant_CYP19) 
dim(Variant) #461 3828

# Clean the environment
rm(ids, vcf,
   Variant_ZBTB38, Variant_HMGA2, Variant_CYP19, 
   vcf_ZBTB38, vcf_HMGA2, vcf_CYP19,
   AF_ZBTB38, AF_HMGA2, AF_CYP19)

# Check if there are any NAs 
length(Variant[is.na(Variant)]) #0

# Remove the columns which has only the same values for each individuals meaning that the uniqe variables number is 1
Variant <- Variant[,apply(Variant, 2, function(x) length(unique(x))) != 1]
dim(Variant) #461 3408

#Check the values in the data 
Variant[1,]

# Convert the values to a numeric value and assign the alternative forms for more than 2 to NAs 
Variant <- ifelse(Variant =='0|0', 0, 
                  ifelse(Variant =='0|1', 1, 
                         ifelse(Variant =='1|0', 1, 
                                ifelse(Variant =='1|1', 2, NA))))

# look for NAs 
length(which(colSums(is.na(as.matrix(Variant))) > 0)) # 22 columns have NAs

# Remove the columns contain NA
Variant <- Variant[, -which(colSums(is.na(as.matrix(Variant))) > 0)]
dim(Variant) #461 3386

# Data check
# Calculate the number of 0s which means homozygous for reference sequence found each column
f <- data.frame()
for (i in 1:dim(Variant)[2]) { num <- length(grep(0, Variant[,i])) 
  temp <- data.frame(num)
  f <- rbind(f, temp)
  }
head(f)

#Calculate the frequency of the O's
frequency <- as.data.frame(apply(f, 1, function(x) x=x/dim(Variant)[1]))
colnames(frequency) <- "freq"
head(frequency) # Frequencies are under 0.999 

# Frequency of the genotype of the data I am ganna use for analyse 
hist(frequency$freq, main = "Allele frequency spectrum", 
     xlab="Frequency in Population", 
     ylab = "Proportion of Variants", 
     col = "grey")

# Check if the high freq 0s are not present in the data == YES all removed before.
which(lapply(apply(Variant, 2, table), length) == 1) # 0 

#clean the environment
rm(f, i, temp, num, frequency)

##*********************************************
## Another filter for HWE and MAF ----
##*********************************************
# Convert the Variant Matrix to snp.matrix class
snps <- new("SnpMatrix", Variant)
snps

# Check the final matrix statistics
summary(snps)

# Retrieve the column statistics
snpsum <- col.summary(snps)
summary(snpsum)
length(snpsum$z.HWE^2) ##3386

# Check the multiple allelle frequescies (MAF) and HWE 
par(mfrow = c(1, 2))
hist(snpsum$MAF)
hist(snpsum$z.HWE)
par(mfrow = c(1, 1))

#Convert z values of HWE to p values 
pvalue = pnorm(-abs(snpsum$z.HWE))
pvalue

# Filter the snps based on HWE and MFA 
use <- snpsum$MAF > 0.01 & pvalue < 0.05

sum(use) # 507 (TRUE), To see how many SNPs pass this filter

Variant <- Variant[,use]
dim(Variant) #461 507
length(colnames(Variant)) #507
length(unique(colnames(Variant))) #485, column IDs are not unique because we combine three vcf file for three different genes.
colnames(Variant) <- 1:507 # make them unique

# clean the environment
rm(pvalue, snps, use, snpsum)

##*********************************************

######### Retrieve the population names ----

##*********************************************
## Get the population names for each individual for all data in the 1000 genome
ped = read.table(url("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"),sep="\t",header=TRUE,row.names=2)
head(ped)
class(ped)
#write.csv(ped,'ped.csv')

# Extract the population names for the ids
ids2 <- ped[rownames(Variant), 6, drop=FALSE] #
ids2
length(ids2$Population) #461
unique(ids2)
summary(ids2)

# Check if the individuals names are in same order
ids2 <- as.data.frame(cbind(rownames(Variant), ids2)) 
head(ids2) #Yes 
ids2$Population # Use this when you need population names

##*********************************************

## PCA ----

##*********************************************

pca <- prcomp(Variant, center = TRUE, scale = TRUE)

# Quantile plot
#png("quantile_variant.png", units="in", width=5, height=5, res=300)
qqnorm(pca$x[,1]); qqline(pca$x[,1], col = "steelblue", lwd = 2)
#dev.off()

#png("hist_quantile_variant.png", units="in", width=5, height=5, res=300)
hist(pca$x[,1], main = "Histogram of quantiles", xlab= "Quantiles")
#dev.off()

# Importance of PCA components
#png("PCA_explained_variant.png", units="in", width=5, height=5, res=300)
plot(summary(pca)$importance[3,]*100, xlab = "Number of components", ylab = "Explained %")
#dev.off()

#Scree Plot
#png("Scree_plot_variant.png", units="in", width=5, height=5, res=300)
barplot((pca[[1]]^2)[1:8],
        main= "Scree Plot", 
        names.arg=1:8, 
        xlab="PCs", 
        ylab="variance")
#dev.off()

#Biplot
# make groups for the populations
unique(ids2$Population) # CHB ESN FIN MXL PJL
length(which(ids2$Population == "CHB")) #103
length(which(ids2$Population == "ESN")) #99
length(which(ids2$Population == "FIN")) #99
length(which(ids2$Population == "MXL"))# 64
length(which(ids2$Population == "PJL")) #96
pop <- c(rep("CHB", 103), rep("ESN", 99), rep("FIN", 99), rep("MXL", 64), rep("PJL", 96))
num_pop <- c(rep(1, 103), rep(2, 99), rep(3, 99), rep(4, 64), rep(5, 96))

rm(ped, ids2)

#png("Variant_biplot.png", units="in", width=5, height=5, res=300)
ggbiplot(pca, ellipse=TRUE, var.axes=FALSE, groups=pop) + 
  scale_colour_manual(name="Populations", , values= c("forest green", "red4", "dark blue", "orange", "purple")) +
  theme_minimal() +
  theme(legend.direction = 'horizontal', legend.position = 'top') 
#dev.off()

# Check if there are any relations between PCA components and populations for regression 
pcs <- as.data.frame(pca$x)
plot(pcs$PC1, num_pop)
plot(pcs$PC2, num_pop)
plot(pcs$PC3, num_pop)
plot(pcs$PC4, num_pop)

cor(pcs$PC1, num_pop)
cor(pcs$PC2, num_pop)
cor(pcs$PC3, num_pop)

# there are no mild or strong corelations between populations and pca components
# therefore regression analysis would not be result well 

rm(pcs)

##*********************************************

## Dis-similarity ----

##*********************************************
# Check the dis-similarity to understand if clustering is possible and there are dissimarity between the individuals
#png("distance_matrix.png", units="in", width=5, height=5, res=300)
fviz_dist(dist(Variant), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
#dev.off()

##*********************************************

## Check data tendency for clustering ----

##*********************************************
hopkins(Variant, 5, byrow = F, header = F) # 0.3639521
#Null Hypothesis (Ho) : Data points are generated by non-random, uniform distribution (implying no meaningful clusters)
#Alternate Hypothesis (Ha): Data points are generated by random data points (presence of clusters)
#If H>0.5, null hypothesis can be rejected and it is very much likely that data contains clusters. If H is more close to 0, then data set doesn’t have clustering tendency.
# 0.3639521 < 0.5 Therefore Data has no meaningful clusters

##*********************************************

## Clustering----

##*********************************************

Variant_pop <- Variant
rownames(Variant_pop) <- pop

# Hierarchical Clustering
hc <- hclust(dist(Variant_pop), method = "complete")
dend_hc <- as.dendrogram(hc)

# Agglomerative Nesting (Hierarchical Clustering)
ag <- agnes(dist(Variant_pop), method = "complete")
dend_ag <- as.dendrogram(ag)

### Divisive Clustering
dv <- diana(dist(Variant_pop))
dend_dv <- as.dendrogram(dv)

# obtain the dendogram leaf names as numbered group names to check the strenght of cluster measures 
hclust <- num_pop[order.dendrogram(dend_hc)]
hagnes <- num_pop[order.dendrogram(dend_ag)]
hdiana <- num_pop[order.dendrogram(dend_dv)]

hcm <- adjustedRand(hclust, num_pop)
ham <- adjustedRand(hagnes, num_pop)
hdm <- adjustedRand(hdiana, num_pop)

compare <- data.frame(hclust = hcm, hagnes=ham, hdiana=hdm)
compare

# hclust has the highest scores for all the measures so plot the hclust
dend <- dend_hc %>%
  #color_branches(k = length(unique(labels(dend)))) %>%
  set("branches_lwd", c(2,1,2)) %>%
  set("branches_lty", c(1,2,1))

# Labels
labels_colors(dend) <-
  rainbow_hcl(length(unique(labels(dend))))[sort_levels_values(
    as.numeric(num_pop)[order.dendrogram(dend)])]

par(mar = rep(0,4))

#png("circlize.png", units="in", width=10, height=10, res=300)
circlize_dendrogram(dend)
legend("topleft", legend = unique(labels(dend)), fill = rainbow_hcl(length(unique(labels(dend)))))
#dev.off()

par(mar = rep(5,4))

##*********************************************

## K-Means Clustering

##*********************************************
kmean <- kmeans(scale(Variant_pop[,1:461]), 5)
kpam <- pam(scale(Variant_pop[,1:461]), 5)	# develop medoids only for ability to plot kmeans

kpam$medoids <- kmean$centers   # this is done to make it possible to plot kmeans
kpam$clustering <- kmean$cluster

#png("kmean.png", units="in", width=10, height=10, res=300)
clusplot(kpam, labels=5, sub='', main = "5-means")
#dev.off()

med <- pam(scale(Variant_pop[,1:461]), 5)

#png("PAM.png", units="in", width=10, height=10, res=300)
clusplot(med, labels=5, sub='', main = "5-medoids")
#dev.off()

table(med$clustering)

# Aggrement of the clusters 
head(med$clustering)
head(num_pop)

PAM <- adjustedRand(med$clustering, num_pop)

kmean <- adjustedRand(kmean$cluster, num_pop)

compare <- cbind(compare, PAM, kmean)

compare

#write.csv(compare,'clustering_compare.csv')

# kmeans plotted againts the PCA scores 
rownames(Variant_pop) <- pop
k5 <- kmeans(Variant_pop, centers = 5)

Component1 <- pca$x[,1]
Component2 <- pca$x[,2]

#png("KMEANS_GGPLOT.png", units="in", width=10, height=10, res=300)
Variant %>%
  as_tibble() %>%
  mutate(cluster = k5$cluster, pop = pop) %>%
  ggplot(aes(Component1, Component2, color = factor(cluster), label = pop)) +
  geom_text()
#dev.off()


