# working dir: ymweng@denali:/data/Yiming/GenomeProject/ingens_SNP/snmf/
# R envrionment
library(LEA) 
library(ggplot2)
library(randomcoloR)
library(vcfR)
library(adegenet)
library(maps)
library(mapplots)

setwd("C:/Users/wengz/Box/Spiroplasma/population_structure/sNMF")
# convert .ped file to .geno file
ped2geno("spiroplasma_snps.ped", "spiroplasma_snps.geno", force = TRUE)

# Create a project
spirosnmf = NULL

# Runs with K between 1 and 12 with cross-entropy and 10 repetitions
spirosnmf = snmf("spiroplasma_snps.geno", 
                  K=1:12, CPU = 28, entropy = TRUE, repetitions = 10, project = "new")

## plot cross-entropy criterion of all runs of the project
pdf("spiroplasma_snps_cross_entropy.pdf")
plot(spirosnmf, lwd = 1, col = "black", pch=1, type='o', main = "Cross-entropy \n versus number of clusters")
dev.off()


# to average the ten runs for each K by using for loop
ce <- NULL
for (k in 1:12) {
  # k <- 1
  ce[k] <- mean(cross.entropy(spirosnmf, K = k))
}

# select the run with the lowest cross-entropy
K1 <- which.min(ce)

# setup k=6 for downstream analyses
K1 <- 3

# Next step is to create ancestry coefficients barplot (Q matrix)
# We start with constructing specimen table with population information
# import geographical information
pop.info <- read.csv("spiro_all_data.csv", fileEncoding="UTF-8-BOM")

# Import SNP data as within vcf format
imputed <- read.vcfR("C:/Users/wengz/Box/Spiroplasma/population_structure/PCA/spiroplasma_snps_final.vcf.gz", verbose = FALSE)

# Convert SNP data from vcf to genlight format

imputed.genlight <- vcfR2genlight(imputed)
n <- nrow(as.matrix(imputed.genlight))
L <- ncol(as.matrix(imputed.genlight))

# Extract samples who have SNP data from the pop.info dataset
subpop <- pop.info[which(pop.info$sampleID %in% indNames(imputed.genlight)), ]

# Generate palettes Of optimally distinct clors for K1 (K1=6)
# https://www.rdocumentation.org/packages/randomcoloR/versions/1.1.0/topics/distinctColorPalette
palette <- distinctColorPalette(K1)
color <- c("yellowgreen", "bisque3", "cyan3", "springgreen3", "gray66", "lightpink", "tan4", "slateblue2", "dodgerblue4", "orchid3")
color <- c("yellowgreen", "bisque3", "cyan3", "springgreen3", "gray66", "lightpink", "tan4", "slateblue2", "dodgerblue4")
# color <- c("cyan3", "yellow","bisque3", "gray66", "tan4", "orchid3", "lightpink", "yellowgreen", "slateblue2", "springgreen3", "dodgerblue4"
# color <- c("cyan3", "yellow","bisque3", "gray66", "tan4", "orchid3", "lightpink", "yellowgreen")


# To generate Q matrix, the individual admixture coefficient matrix
# Use minimum ce in ten runs of K=6
qmatrix <- Q(spirosnmf, K = K1, run=which.min(ce[K1]))

# Now we create the table by combining population information and Q matrix
# https://www.rdocumentation.org/packages/base/versions/3.5.0/topics/gl
snmf.ret <- data.frame(sampleID = gl(n, K1, labels = indNames(imputed.genlight)),
                       ratio = as.numeric(t(qmatrix)),
                       cluster = gl(K1, 1, n * K1, labels = paste("cluster", seq(1:K1), sep = "")))

# merge subpop data with snmf.ret by specimenID
snmf.ret <- merge(subpop, snmf.ret, by = "sampleID")

# Creat a barplot to present the ancentry cofficents wihh population lable 
# http://ggplot2.tidyverse.org/reference/scale_manual.html
# http://ggplot2.tidyverse.org/reference/geom_bar.html
x11()
ggplot(snmf.ret, aes(x = sampleID, y = ratio, fill = cluster)) + 
  scale_fill_manual(values = color )+
  theme(axis.text.x=element_blank()) +
  geom_bar(stat = "identity") + 
  facet_grid(~siteID, switch = "x", scales = "free_x", space = "free_x")

pdf("/data/Yiming/GenomeProject/ingens_SNP/snmf/lcs_bqsr_snmf_11.pdf")
snmf.p <- ggplot(snmf.ret, aes(x = sampleID, y = ratio, fill = cluster)) + 
  scale_fill_manual(values = color )+
  theme(axis.text.x=element_blank()) +
  geom_bar(stat = "identity") + 
  facet_grid(~siteID, switch = "x", scales = "free_x", space = "free_x")
print(snmf.p)
dev.off()


# Next step is to combine q matrix to the population map
# Create new data contains sample ID, population location and Q
q.mat <- data.frame(sampleID = indNames(imputed.genlight), qmatrix)
pop.snmf <- merge(subpop, q.mat, by = "sampleID")

# dwar the map with the info from the table created above
# creat a region with longitude and latitude
x11()
plot(NA,NA, ylim = c(36.3, 38), xlim = c(-119.5, -117.8),
     xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90",  fill = TRUE)
pop.list <- unique(snmf.ret$siteID)
for(i in 1:length(pop.list)) {
  pop.ave <- apply(pop.snmf[which(pop.snmf$siteID == pop.list[i]), -(1:8)],
                   2, mean)
  lat.ave <- mean(pop.snmf[which(pop.snmf$siteID == pop.list[i]), 7])
  lon.ave <- mean(pop.snmf[which(pop.snmf$siteID == pop.list[i]), 8])
  add.pie(z = pop.ave, x = lon.ave, y = lat.ave,
          col = color, labels = NA, label.dist = 2, radius = 0.05)
}


pdf("/data/Yiming/GenomeProject/ingens_SNP/snmf/lcs_bqsr_snmf_k11map.pdf")
plot(NA,NA, ylim = c(36.3, 38), xlim = c(-119.5, -117.8),
     xlab = "Longitude", ylab = "Latitude", type = "n")


# apply map to the region defined above
map(add = T, col = "grey90",  fill = TRUE)

# Add the data to the map
# https://www.rdocumentation.org/packages/base/versions/3.5.0/topics/apply
pop.list <- unique(snmf.ret$PopulationCode)
for(i in 1:length(pop.list)) {
  pop.ave <- apply(pop.snmf[which(pop.snmf$PopulationCode == pop.list[i]), -(1:5)],
                   2, mean)
  lat.ave <- mean(pop.snmf[which(pop.snmf$PopulationCode == pop.list[i]), 3])
  lon.ave <- mean(pop.snmf[which(pop.snmf$PopulationCode == pop.list[i]), 4])
  add.pie(z = pop.ave, x = lon.ave, y = lat.ave,
          col = color, labels = NA, label.dist = 2, radius = 0.05)
}
dev.off()

