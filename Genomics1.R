
setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/BINF6110_Genomic_Methods")

# Triticum aestivum RefSeq: 22 chromosomes are present. There are no scaffolds in the headers only chromosomes but in the metadata they indicate 22 scaffolds are present and calculated Scaffold N50 is 709773743. Overal 
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/519/105/GCA_900519105.1_iwgsc_refseqv1.0/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna.gz
# gunzip GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna.gz
# Metadata : https://www.ncbi.nlm.nih.gov/assembly/GCA_900519105.1
#1-	How big is the genome?
# grep "^[^>]" GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna | tr -d "\n" | wc -c
# 14,547,261,565  # 14,547,261,565
#2-	How many scaffolds are included?  Are there both chromosomes and scaffolds?
# egrep "^[>+ ]" GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna | wc -l
# 22 # There are 22 chromosomes
# 3-	What is the average scaffold/chromosome length?  Plot scaffold length from largest to smallest.
# awk '/^>/ {if (seqlen) {print seqlen}; printf $0"\t"; seqlen=0;next; }{ seqlen += length($0)}END{print seqlen}' GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna > wheat_chromosome.txt

# 4-	What is the n50 given in the paper/metadata?  Can you recreate this result with your own calculations?

N50 <- function(filename){
  colnames = c ("name", "length")
  scaffold = read.table(filename, col.names = colnames, sep ="\t")
  scaffolds <- c(sort(scaffold$length, decreasing = T))
  for (i in cumsum(as.numeric(sort(scaffold$length, decreasing = T)))) {
    if (i < sum(scaffold$length)/2) {vector <- c(vector, i) } }
  result = scaffolds[length(vector)]
  result1 = plot(sort(scaffold$length), cumsum(as.numeric(sort(scaffold$length)), decreasing=T))
  result2= plot(sort(scaffold$length), cumsum(as.numeric(sort(scaffold$length))))
  result3= plot(1:length(scaffold$length), cumsum(as.numeric(sort(scaffold$length, decreasing=T))))
  return(result) 
  return(result1)
  }

N50("wheat_chromosome.txt") # 709,773,743

#5-	What percentage of the genome is N or n?
Nn_count=`grep -o "[Nn]" Org_Aviadenovirus_complete_genome_fetch.fasta | wc -l`
total_count=`grep "^[^>]" Org_Aviadenovirus_complete_genome_fetch.fasta | tr -d "\n" | wc -c` | 
echo "scale=5 ; $Nn_count / $total_count * 100" |bc

# 6-	If your genome has both uppercase and lowercase letters (denoting non-repetitive and repetitive regions, respectively), what is the ratio of repetitive to unique sequence?
lower_count=`grep -o [atgc] Org_Aviadenovirus_complete_genome_fetch.fasta | wc -l`
upper_count=`grep -o [ATGC] Org_Aviadenovirus_complete_genome_fetch.fasta | wc -l`
echo "scale = 5; $lower_count / $upper_count" |bc

# Drosophila melanogaster : if I grep scaffolds I cannot find the same result with metadata. if I use scaffolds and chromosomes together then I get the same result. They used scaffolds and chromosomes together to calculate N50
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
# gunzip GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
# Metadata: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/
# 1-	How big is the genome? 
# grep "^[^>]" GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna | tr -d "\n" | wc -c
# 143,726,002 
# 2-	How many scaffolds are included?  Are there both chromosomes and scaffolds?
# egrep "^[>+ ]" GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna | wc -l 
# 1870 
# egrep "^[>+ ]" GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna | grep "Scaffold" | wc -l
# 47 There are 47 Scaffolds 891 chromosomes
# 3-	What is the average scaffold/chromosome length?  Plot scaffold length from largest to smallest.
# awk '/^>/ {if (seqlen) {print seqlen}; printf $0"\t"; seqlen=0;next; }{ seqlen += length($0)}END{print seqlen}' GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna > drosophila_scaffold_sorted.txt
# 4-	What is the n50 given in the paper/metadata?  Can you recreate this result with your own calculations?
N50("drosophila_scaffold_sorted.txt") # 25,286,936
# 5-	What percentage of the genome is N or n?


colnames = c ("name", "length")
scaffolds = read.table("wheat_chromosome.txt", col.names = colnames, sep ="\t")
mean(scaffolds$length)

plot(sort(scaffolds$length), cumsum(as.numeric(sort(scaffolds$length))), xlab = "chromosome/scaffold length", ylab = "cumulative sum")

plot(1:length(scaffolds$length), cumsum (as.numeric(sort(scaffolds$length, decreasing=T))), , xlab = "chromosome/scaffold length", ylab = "cumulative sum")

hist(sort(scaffolds$length, decreasing = F), xlab = "Scaffold length", main = "Histogram of Scaffolds")

plot(sort(scaffolds$length, decreasing = F),  main = "Scaffolds")

# Triticum aestivum 3B chromosome WGS project CBUC00000000 : all scaffolds 
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/077/335/GCA_001077335.1_ASM107733v1/GCA_001077335.1_ASM107733v1_genomic.fna.gz
# gunzip GCA_001077335.1_ASM107733v1_genomic.fna.gz
#metadata : https://www.ncbi.nlm.nih.gov/assembly/GCA_001077335.1/
#1-	How big is the genome?
# grep "^[^>]" GCA_001077335.1_ASM107733v1_genomic.fna | tr -d "\n" | wc -c
# 58,502,153 #58,502,153
#2-	How many scaffolds are included?  Are there both chromosomes and scaffolds?
# 1450 # There are only 1450 scaffolds
# 3-	What is the average scaffold/chromosome length?  Plot scaffold length from largest to smallest.
# awk '/^>/ {if (seqlen) {print seqlen}; printf $0"\t"; seqlen=0;next; }{ seqlen += length($0)}END{print seqlen}' GCA_001077335.1_ASM107733v1_genomic.fna > wheat_scaffold_sorted.txt
# 4-	What is the n50 given in the paper/metadata?  Can you recreate this result with your own calculations?
N50("wheat_scaffold_sorted.txt") # 213,933
# 5-	What percentage of the genome is N or n?



