#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Emine Ozsahin
# March 10, 2020
#Genomics Methods
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library("reshape2")
library("ggplot2")
library("splitstackshape")
library("data.table")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fig.2. Allele frequency spectrum freebayes ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Freebayes
freq_freebayes <- read.table ("freebayes_allele_freq.frq", sep=" ")
head(freq_freebayes)

freq_freebayes <- cSplit(data.table(freq_freebayes), "V1", "\t|:", fixed = FALSE)
head(freq_freebayes)
freq_freebayes <- freq_freebayes[-1, -(9:16)]
head(freq_freebayes)

colnames_freq <- c("CHROM", "POS", "N_ALLELES",  "N_CHR", "ALLELE1", "ALLELE1_FREQ", "ALLELE2", "ALLELE2_FREQ")

colnames(freq_freebayes) <- colnames_freq
freq_freebayes <- cbind(freq_freebayes, Variant_Caller= rep("freebayes", dim(freq_freebayes)[1]))
head(freq_freebayes)

freq_freebayes$ALLELE1_FREQ <- as.numeric(as.character(freq_freebayes$`ALLELE1_FREQ`))
head(freq_freebayes$ALLELE1_FREQ)

freq_freebayes$ALLELE2_FREQ <- as.numeric(as.character(freq_freebayes$`ALLELE2_FREQ`))
head(freq_freebayes$ALLELE2_FREQ)

tiff("hist1.tiff", units="in", width=5, height=5, res=300)
hist(freq_freebayes$`ALLELE1_FREQ`, main = "Allele frequency spectrum", xlab="Frequency in Population", ylab = "Proportion of SNPs", col = "dark blue")
hist(freq_freebayes$`ALLELE2_FREQ`, add=T, col = scales::alpha("pink", .5))
dev.off()

#bcf_mpileup
freq_bcf <- read.table ("bcftools_allele_freq.frq", sep=" ")
head(freq_bcf)

freq_bcf <- cSplit(data.table(freq_bcf), "V1", "\t|:", fixed = FALSE)
freq_bcf <- freq_bcf[-1, -(9:12)]
colnames(freq_bcf) <- colnames_freq
freq_bcf <- cbind(freq_bcf, Variant_Caller= rep("bcftools_mpileup", dim(freq_bcf)[1]))
head(freq_bcf)

freq_bcf$ALLELE1_FREQ <- as.numeric(as.character(freq_bcf$`ALLELE1_FREQ`))
head(freq_bcf)
dim(freq_bcf)

#merge
freq <- rbind(freq_bcf,freq_freebayes)
head(freq)
tail(freq)
dim(freq)

hist(freq2$ALLELE1_FREQ)

hist(freq$ALLELE1_FREQ)
#plot

tiff("histggplot.tiff", units="in", width=5, height=5, res=300)

ggplot(data = freq, aes(ALLELE1_FREQ)) + 
  geom_histogram(col="dark blue",
                 aes(col= "pink", fill=Variant_Caller), 
                 position = "dodge") + 
  labs(title="Allele frequency spectrum", 
       x="Frequency in Population", 
       y="Proportion of SNPs")
  

dev.off()

freq2 <- freq[!(which(freq$ALLELE1_FREQ==0)),]

tiff("histggplot_2.tiff", units="in", width=5, height=5, res=300)

ggplot(data = freq2, aes(ALLELE1_FREQ)) + 
  geom_histogram(col="dark blue",
                 aes(col= "pink", fill=Variant_Caller), 
                 position = "dodge") + 
  labs(title="Allele frequency spectrum", 
       x="Frequency in Population", 
       y="Proportion of SNPs")

dev.off()

#+++++++++++++++++++++++++++++++++++++++++++
#  Fig.3. Depth
#+++++++++++++++++++++++++++++++++++++++++++

depth_freebayes <- read.table("freebayes_depth.txt")[-1,]
colnames_depth <- c("Individual", "N_Sites", "Mean_Depth")
colnames(depth_freebayes) <- colnames_depth

depth_freebayes$Mean_Depth <- as.numeric(as.character(depth_freebayes$Mean_Depth))
depth_freebayes

depth_bcf <- read.table("bcftools_depth.txt") [-1,]
colnames(depth_bcf) <- colnames_depth

depth_bcf$Mean_Depth <- as.numeric(as.character(depth_bcf$Mean_Depth))
depth_bcf

depth <- data.frame(Individual= depth_bcf$Individual, bcf_mpileup= depth_bcf$Mean_Depth, Freebayes=depth_freebayes$Mean_Depth)

depth <- reshape2::melt(depth, id="Individual")  
head(depth)

colnames(depth) <- c("Individual", "Varient_Caller", "Depth")

tiff("depth.tiff", units="in", width=5, height=5, res=300)
ggplot() + geom_col(data = depth, aes(x = Individual, y = Depth, fill = Varient_Caller), position = "dodge") + coord_flip() + theme(legend.position = "top")

dev.off()

mean (depth_freebayes$Mean_Depth)
quantile(depth_freebayes$Mean_Depth)
median(depth_freebayes$Mean_Depth)

mean (depth_bcf$Mean_Depth)
quantile(depth_bcf$Mean_Depth)
median(depth_bcf$Mean_Depth)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fig.1. Total number of SNPs for per individual before filters. ----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SNPs_bcf <- read.table("bcftools_SNP_numbers.txt")

SNPs_freebayes <- read.table("freebayes_SNP_numbers.txt")

SNPs <- data.frame(Individuals = SNPs_bcf[-11,1], Samtools_mpileup = SNPs_bcf[-11,7], Freebayes = SNPs_freebayes[-11,7])

SNPs_count <- melt(SNPs, id="Individuals")  

colnames <- c("Individuals", "Varient_Caller", "Number_of_SNPs")

colnames(SNPs_count) <- colnames

SNPs_count

png("sample.png", 490, 350)

theme_update(plot.title = element_text(hjust = 0.5)) # This has to be run once to make the title centered for all of the ggplots
tiff("snp.tiff", units="in", width=5, height=5, res=300)
ggplot() + geom_col(data = SNPs_count, aes(x = Individuals, y = Number_of_SNPs, fill = Varient_Caller), position = "dodge") + coord_flip() + theme(legend.position = "top")

dev.off()

#+++++++++++++++++++++++++++++++++++++++++++
#  Fig.4. VennDiagram - overlapped SNPs
#+++++++++++++++++++++++++++++++++++++++++++

# Load library
#install.packages("VennDiagram")
library(VennDiagram)

#bcftools stat
# SN	0	number of SNPs:	9410 bcf only
# SN	1	number of SNPs:	62273 freebayes only
# SN	2	number of SNPs:	30548 intersect

# A more complicated diagram Demonstrating external area labels
tiff(filename = "Pairwise_Venn_diagram.tiff", compression = "lzw");
venn.plot <- draw.pairwise.venn(
  area1 = 92821,
  area2 = 39958,
  cross.area = 30548,
  category = c("Freebayes", "bcftools mpileup"),
  fill = c("yellow", "pink"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  margin = 0.1
  #cat.pos = c(285, 105),
  #cat.dist = 0.09,
  #cat.just = list(c(-1, -1), c(1, 1)),
  #ext.pos = 30,
  #ext.dist = -0.05,
  #ext.length = 0.85,
  #ext.line.lwd = 2,
  #ext.line.lty = "dashed"
);
dev.off();
grid.draw(venn.plot);
grid.newpage();


# SN	0	number of indels:	5056 5390
# SN	1	number of indels:	3576  3910
# SN  2 number of indels:	334

tiff(filename = "Pairwise_Venn_diagram_indels.tiff", compression = "lzw");
venn.plot <- draw.pairwise.venn(
  area1 = 5390,
  area2 = 3910,
  cross.area = 334,
  category = c("bcftools mpileup", "Freebayes"),
  fill = c("yellow", "pink"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  margin = 0.2
);
grid.draw(venn.plot);
dev.off();

#+++++++++++++++++++++++++++++++++++++++++++
#  Fig.5. Depth for loci
#+++++++++++++++++++++++++++++++++++++++++++
depth_site <- read.table("bcftools_vs_freebayes_depth.ldepth")
sum(is.na(depth_site$V3))

head(depth_site)

depth_site$V3 <- as.numeric(as.character(depth_site$V3))

tiff(filename = "depth_by_site.tiff", compression = "lzw")

hist(depth_site$V3, main = "Depth by site", ylab= "Number of locus", xlab="depth")

dev.off()



#+++++++++++++++++++++++++++++++++++++++++++
#  Fig.5. 
#+++++++++++++++++++++++++++++++++++++++++++
read.table("bcftools_vs_freebayes.diff.sites_in_files")

library(tidyverse)

vcf <- read.csv("compare_bcfstat", sep = " ")

head(vcf)

vcf <- vcf %>% 
       filter(str_detect(X., "QUAL"))

colnames(vcf)

vcf <- cSplit(data.table(vcf[,1]), "V1", "\t", fixed = FALSE)

colnames(vcf) <- c("QUAL", "id",	"Quality", "number of SNPs", "number of transitions (1st ALT)", "number of transversions (1st ALT)",	"number of indels")

vcf<- vcf[,c(2,3)]
vcf$Quality <- as.numeric(as.character(vcf$Quality))
vcf$id <- as.numeric(as.character(vcf$id))

vcf_bcf <- c(vcf[which(vcf$id==0),])
hist(vcf_bcf$Quality)

vcf_free <- c(vcf[which(vcf$id==1),])
hist(vcf_free$Quality)

vcf_com <- c(vcf[which(vcf$id==2),])
hist(vcf_com$Quality)


tiff("hist1.tiff", units="in", width=5, height=5, res=300)

hist(vcf_bcf$Quality, main = "", xlab="", ylab = "", col = "red")

hist(vcf_free$Quality, add=T, col = scales::alpha("pink", .5))

hist(vcf_com$Quality, add=T, col = scales::alpha("yellow", .5))

dev.off()


#+++++++++++++++++++++++++++++++++++++++++++
## Examples
#+++++++++++++++++++++++++++++++++++++++++++


# ggplot(data=SNPs_count,
#        aes(x=Individuals, y=Number_of_SNPs, colour=Varient_Caller)) + geom_line()
# 
# ggplot(data=SNPs_count,
#        aes(x=Individuals, y=Number_of_SNPs, colour=Varient_Caller)) + geom_col()
# 
# ggplot(data=SNPs_count,
#        aes(x=Individuals, y=Number_of_SNPs, colour=Varient_Caller)) + geom_count()
# 
# ggplot() + geom_col(data = SNPs_count, aes(x=Individuals, y=Number_of_SNPs, colour=Varient_Caller), position = "dodge")

#pdf("snp_plot.pdf", 7, 5)


# # Make long style 
# freq_h <- reshape2::melt(freq, id="Individuals")
# colnames(freq_h) <- c("Individuals", "Variant_Caller", "value" )
# head(freq_h)
# dim(freq_h)