# 2021.08.04
**Population structure of the *Spiroplasma***

1. Merge the vcf files of the two contigs into one with Picard, and index the combinded vcf (spiroplasma_raw.vcf.gz)
```
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/GenomicsDBImport/
# combind the two vcf files
java -jar /media/Data1/Yiming/software/picard.jar MergeVcfs \
          I=tig00001416.vcf.gz \
          I=tig00009077.vcf.gz \
          O=spiroplasma_raw.vcf.gz

# index spiroplasma_raw.vcf.gz
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk IndexFeatureFile -I spiroplasma_raw.vcf.gz
```

2. To filter the SNPs in the vcf file, check the quality score distributions first:
  - extract the quality values from the raw vcf file
```
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/
mkdir vcf_filter
mv /media/Jade/YMW/spiroscan/vcf/GenomicsDBImport/spiroplasma_raw.vcf* ./vcf_filter

# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/vcf_filter
echo -e "QD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum" >> spiroplasma_quality
QD=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -o "QD=.*;" | cut -d ";" -f 1 | cut -d "=" -f 2)
FS=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -o "FS=.*;" | cut -d ";" -f 1 | cut -d "=" -f 2)
SOR=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -Po "SOR=.*?\t" | cut -d "=" -f 2)
MQ=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -o "MQ=.*;" | cut -d ";" -f 1 | cut -d "=" -f 2)
MQRankSum=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -o "MQRankSum=.*;" | cut -d ";" -f 1 | cut -d "=" -f 2)
ReadPosRankSum=$(zcat spiroplasma_raw.vcf.gz | grep -v "#" | grep -o "ReadPosRankSum=.*;" | cut -d ";" -f 1 | cut -d "=" -f 2)

paste <(echo "$QD") <(echo "$FS") <(echo "$SOR") <(echo "$MQ") <(echo "$MQRankSum") <(echo "$ReadPosRankSum") --delimiters $'\t' >> spiroplasma_quality
sed -i 's/\t\t/\t/' spiroplasma_quality
```
  - visualize the distribution of each quality score, and decide the cutoff
  - QD = 5.0
  - FS = 10.0
  - SOR = 3.0
  - MQ = 42.0
  - MQRankSum = -2.0
  - ReadPosRankSum = -2.0
  - **DP = 630 ** (=210\*3: the mean depth for site should be at least 3x)
![](@attachment/Clipboard_2021-08-04-09-46-13.png)
![](@attachment/Clipboard_2021-08-04-09-47-22.png)
![](@attachment/Clipboard_2021-08-04-09-48-24.png)
![](@attachment/Clipboard_2021-08-04-09-51-23.png)
![](@attachment/Clipboard_2021-08-04-09-52-27.png)

  - apply the filters to the vcf:
```
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/vcf_filter
# seperate the SNPs and indels
# extract SNPs and write to spiroplasma_snp_raw.vcf.gz
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk SelectVariants -V spiroplasma_raw.vcf.gz -select-type SNP --exclude-filtered -O spiroplasma_snp_raw.vcf.gz

# extract indels and write to spiroplasma_indels_raw.vcf.gz
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk SelectVariants -V spiroplasma_raw.vcf.gz -select-type INDEL -O spiroplasma_indels_raw.vcf.gz

# hard-filtering the SNPs file
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk VariantFiltration \
    -V spiroplasma_snp_raw.vcf.gz \
    -filter "QD < 5.0" --filter-name "QD5" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 10.0" --filter-name "FS10" \
    -filter "MQ < 42.0" --filter-name "MQ42" \
    -filter "MQRankSum < -2.0" --filter-name "MQRankSum-2" \
    -filter "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSum-2" \
    --filter "DP < 610" --filter-name "DP3" \
    -O spiroplasma_snps_filtered.vcf.gz

# remove the filtered site for snp vcf
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk SelectVariants -V spiroplasma_snps_filtered.vcf.gz -select-type SNP --exclude-filtered -O spiroplasma_snps_final.vcf.gz

# hard-filtering the indel file
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk VariantFiltration \
    -V spiroplasma_indels_raw.vcf.gz \
    -filter "QD < 5.0" --filter-name "QD5" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 10.0" --filter-name "FS10" \
    -filter "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSum-2" \
    --filter "DP < 610" --filter-name "DP3" \
    -O spiroplasma_indels_filtered.vcf.gz

# remove the filtered site for indel vcf
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk SelectVariants -V spiroplasma_indels_filtered.vcf.gz -select-type INDEL --exclude-filtered -O spiroplasma_indels_final.vcf.gz
    
```
  - The final snp vcf has 37,381 SNPs in the vcf file
  - The final indel vcf has 4,598 indels in the vcf file

3. PCA of the filtered SNPs
```
# R environment in my DELL laptop
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")

library(gdsfmt)
library(SNPRelate)

setwd("C:/Users/wengz/Box/Spiroplasma/population_structure/PCA/")

# convert vcf to gds
vcf.fn <- "C:/Users/wengz/Box/Spiroplasma/population_structure/PCA/spiroplasma_snps_final.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "spiro.gds", method="copy.num.of.ref")
snpgdsSummary("spiro.gds")
# The total number of samples: 210 
# The total number of SNPs: 37381 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
# The number of valid samples: 210 
# The number of biallelic unique SNPs: 35069

genofile <- snpgdsOpen("spiro.gds")
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

sample.info=read.csv("spiro_all_data.csv", fileEncoding="UTF-8-BOM")
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
```
![](@attachment/Clipboard_2021-08-04-11-34-24.png)

4. sNMF of the SNPs
  - prepare the input file by converting vcf to plink ped, than to geno in R
```
# working dir: working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/
# use vcftools to convert vcf to plink ped
vcftools --gzvcf spiroplasma_snps_final.vcf.gz --plink --out spiroplasma_snps
# note that this only takes biallelic loci
```