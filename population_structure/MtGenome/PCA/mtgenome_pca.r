# R environment in my DELL laptop
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")

library(gdsfmt)
library(SNPRelate)

setwd("C:/Users/wengz/Dropbox/Spiroplasma_project/MtGenome/PCA")

# convert vcf to gds
vcf.fn <- "mtgenome_snps_final_trimmed.vcf"
snpgdsVCF2GDS(vcf.fn, "mtgenome.gds", method="copy.num.of.ref")
snpgdsSummary("mtgenome.gds")
# The total number of samples: 375 
# The total number of SNPs: 199
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
# The number of valid samples: 375
# The number of biallelic unique SNPs: 191

genofile <- snpgdsOpen("mtgenome.gds")
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.9, autosome.only=FALSE)
# Chromosome tig00001416: 46.52%, 87/187
# Chromosome tig00009077: 6.44%, 2,396/37,194
# 2,483 markers are selected in total.

snpset.id <- unlist(unname(snpset))
samp_id <- read.gdsn(index.gdsn(genofile, "sample.id"))

pca <- snpgdsPCA(genofile, 
                 sample.id=samp_id, 
                 snp.id=snpset.id, 
                 autosome.only=FALSE, 
                 maf=0.05,
                 missing.rate=0.25,
                 num.thread=2)
pc.percent <- pca$varprop*100
# Excluding 170 SNPs
# number of SNPs: 2,313

head(round(pc.percent, 2))
# [1] 28.88 13.27  6.49  6.11  5.36  3.25

sample.info <- read.table("C:/Users/wengz/Dropbox/Spiroplasma_project/sample_data.txt", header=T, sep="\t")
pop.info=sample.info[,c(1,2,3)]
names(pop.info)[1] <- "sample.id"
names(pop.info)[3] <- "pop_code"
names(pop.info)[2] <- "lineage_code"
pop.info = data.frame(lapply(pop.info, function(x) {gsub("YMW17-0", "YMW17-", x)}))

tab <- data.frame(sample.id = pca$sample.id,
                  lineage = factor(pop.info$lineage_code)[match(pca$sample.id, pop.info$sample.id)],
                  pop = factor(pop.info$pop_code)[match(pca$sample.id, pop.info$sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)
x11()
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch=16, xlab="eigenvector 1", ylab="eigenvector 2")
legend("topleft", legend=levels(tab$pop), pch=16, col=1:nlevels(tab$pop))

x11()
plot(tab$EV1, tab$EV2, col=as.integer(tab$lineage), pch=16, xlab="eigenvector 1", ylab="eigenvector 2")
legend("topleft", legend=levels(tab$lineage), pch=16, col=1:nlevels(tab$lineage))


x11()
lbls <- paste("PC", 1:6, "\n", format(pc.percent[1:6], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:6], col=tab$lineage, labels=lbls)
