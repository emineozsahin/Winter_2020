# 1 April 2020
# By: Emine Ozsahin

# Sheep breeds: Canadian Rideau Arcot and Spanish Churra sheep breeds

# Gene stable IDs
# File 1: RNA-Seq ̀Results of differential gene expression analysis performed with DESeq2 using abomasum tissue samples
# from resistant and parasite infected Canadian Rideau Arcot sheep using RNA- Sequencing technology.
File1 <- read.table("File1.txt")

head(File1)

DE_CA_Rideo <- rownames(File1)

head (DE_CA_Rideo)

length(DE_CA_Rideo) # 404

length(unique(DE_CA_Rideo)) # 404

# Gene names find the Gene stable IDs
# File 2: RNA-Seq ̀List of genes differentially expressed between resistant and parasite infected animals from liver tissue samples of Spanish Churra sheep.

#Importing dataset
file2 <-read.table("File2.txt",h=F, sep="\t", stringsAsFactors = F)
head(file2)
values_file2 <- unique(file2$V1)
length(values_file2) #100
length(unique(values_file2)) #100
head(values_file2) 
#Getting Entrez IDs
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("btaurus_gene_ensembl", mart)

filter_file2 <- "external_gene_name"
attributes_file2 <- c("ensembl_gene_id","external_gene_name","entrezgene_id")

File2_biomart <- getBM(attributes=attributes_file2, filters=filter_file2, values=values_file2, mart=mart)

DE_SP_Churra <- File2_biomart$ensembl_gene_id[-1] # there is header don't forget
head(DE_SP_Churra)

length(DE_SP_Churra) #97

length(unique(DE_SP_Churra)) #97

# Gene stable IDs
# File 3: Gene Networks   ̀Ensembl IDs for the genes associated
# with the intestinal immune network for IgA production pathway in sheep.

File3 <- read.table("File3.txt")

head(File3$V1)

Gene_Networks_IgA <- File3$V1

length(Gene_Networks_IgA) # 97
length(unique(Gene_Networks_IgA)) #85

# Cordinates
# File 4: Genomic Regions ̀Genomic coordinates for the genes identified in a genome-wide association study (GWAS)
# for fecal egg counts and immunoglobulin-A (IgA) in a population of Canadian Rideau Arcot sheep.
File4 <- read.table("File4.txt", header = TRUE, sep="\t", stringsAsFactors = F)
head(File4)

File4 <- cbind(File4, X=paste(File4$Chromosome.Name, File4$Gene.Start..bp., File4$Gene.End..bp., sep = ":"))
head(File4$X)

length(unique(File4$X)) #990 in total 941 unique

values_file4 <- File4$X
filter_file4 <-"chromosomal_region"
attributes_file4 <- c("ensembl_gene_id")

file4_biomart <- getBM(attributes=attributes_file4, filters=filter_file4, values=values_file4, mart=mart)

GWAS_CA_Rideau_IgA <- file4_biomart$ensembl_gene_id

length(unique(GWAS_CA_Rideau_IgA)) #832
length(GWAS_CA_Rideau_IgA) #832

# SNPs find rs and gene stable names
# File 5: GWAS (Genome-Wide Association Study) ̀Significant SNPs associated with immunoglobulin- A (IgA)
# identified in a population of Spanish Churra sheep.
file5 <- read.table("File5.txt", header = TRUE)
head(file5)

paste(as.character(file5$SNP), sep = ",",collapse = ",")

length(paste(as.character(file5$SNP), sep = ",",collapse = ","))

values_file5 <- paste(file5$CHR, ":", (file5$Coordinate)-500000, ":",
                      (file5$Coordinate)+500000, sep="")

head(values_file5)
length(values_file5) #130
length(unique(values_file5)) #130
filter_file5 <- "chromosomal_region"

attributes_file5 <- c("ensembl_gene_id")


file5_biomart <- getBM(attributes=attributes_file5, filters=filter_file5, values=values_file5, mart=mart)

GWAS_SP_Churra_IgA <- file5_biomart$ensembl_gene_id

head(GWAS_SP_Churra_IgA)

length(GWAS_SP_Churra_IgA)#111

length(unique(GWAS_SP_Churra_IgA)) #111

#Venn Diagram
library("VennDiagram")
png("Variant_biplot.png", units="in", width=5, height=5, res=300)
venn.diagram (x = list(x1 = DE_CA_Rideo, 
                       x2 = DE_SP_Churra,
                       x3 = Gene_Networks_IgA, 
                       x4 = GWAS_CA_Rideau_IgA,
                       x5 = GWAS_SP_Churra_IgA),
              category.names = c("DE_CA_Rideo", 
                                 "DE_SP_Churra",
                                 "Gene_Networks_IgA",
                                 "GWAS_CA_Rideau_IgA",
                                 "GWAS_SP_Churra_IgA"),
              fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), 
              main = "VennDiagram", main.pos = c(0.5, 1.05), main.just = c( 0.5, 1), 
              lwd = 2, lty = "solid", col = "black", 
              alpha = 0.5, rotation.degree = 0, 
              rotation.centre = c(0.5, 0.5), 
              label.col = "black", cex = 2, fontface = "plain", 
              fontfamily = "serif", 
              cat.dist = 0.10, cat.cex = 1, cat.col = "black", 
              cat.fontface = "plain", 
              cat.fontfamily = "serif", cat.prompts = FALSE,
              ext.text = FALSE, euler.d = TRUE, 
              scaled = TRUE, sep.dist = 0.05, 
              offset = 0, height = 6, width = 6, 
              resolution = 1000,
              units = "in", 
              filename = "a.tiff",
              sub.pos = c( 0.5, 1.05), # bu lazim olmayabilir bir bak
              ext.pos = "", ext.percent = "", ext.line.lwd = "", ext.line.lty = "", description = "",  ext.dist = "", ext.length = "")

dev.off()
#are oficial gene name IDs equal to Gene stable IDs?


# listMarts()
# ensembl <-  useMart("ENSEMBL_MART_SNP")
# datasets <- listDatasets(e)
# head(datasets)
# ensembl <- useDataset("btaurus_structvar", ensembl)
# searchFilters(mart = e, pattern = "var")
# searchAttributes(mart = e, pattern = "var")





