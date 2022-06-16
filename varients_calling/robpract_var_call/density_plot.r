# read variant quality scores from spiroplasma raw vcf file
mydat <- read.table("C:/Users/19207/Box/Spiroplasma/varients_calling/spiroplasma_quality", header=T, sep="\t", fill=T)

# density data for QD
QD <- na.omit(mydat$QD)
d.QD <- density(QD)
plot(d.QD)
polygon(d.QD, col="#FF000022")
abline(v=5,lwd=1, lty=2)

# density data for FS
FS <- na.omit(mydat$FS)
d.FS <- density(FS)
plot(d.FS, xlim=c(1, 50))
polygon(d.FS, col="#FF000022")
abline(v=10, lwd=1, lty=2)
hist(FS)
# density data for SOR
SOR <- na.omit(mydat$SOR)
d.SOR <- density(SOR)
plot(d.SOR, xlim=c(0, 6))
polygon(d.SOR, col="#FF000022")
abline(v=3, lwd=1, lty=2)

# density data for MQ
MQ <- na.omit(mydat$MQ)
d.MQ <- density(MQ)
plot(d.MQ, xlim=c(10, 70))
polygon(d.MQ, col="#FF000022")
abline(v=42.0, lwd=1, lty=2)

# density data for MQRankSum
MQRankSum <- na.omit(mydat$MQRankSum)
d.MQRankSum <- density(MQRankSum)
plot(d.MQRankSum, xlim=c(-3, 3))
polygon(d.MQRankSum, col="#FF000022")
abline(v=-2.0, lwd=1, lty=2)

# density data for ReadPosRankSum
ReadPosRankSum <- na.omit(mydat$ReadPosRankSum)
d.ReadPosRankSum <- density(ReadPosRankSum)
plot(d.ReadPosRankSum, xlim=c(-10, 10))
polygon(d.ReadPosRankSum, col="#FF000022")
abline(v=-2.0, lwd=1, lty=2)


# check match/mismatch coverages
# For vqsr match/mismatch out put of callerror3.sh 
library(ggplot2)
mismatch <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/error_rate/mismatch_rate/vqsr_mismatches", header=T, sep="/t")
match <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/error_rate/mismatch_rate/vqsr_matches", header=T, sep="/t")
DP <- rbind(match, mismatch)
mismatch.type <-rep("mismatch",times=nrow(mismatch))
match.type <-rep("match",times=nrow(match))
type <- c(match.type, mismatch.type)
DP.density <- cbind(DP, type)

ggplot(DP.density, aes(x = DP, fill = type)) + 
  geom_density(alpha = 0.5, bw=0.4) + 
  xlim(0,12) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# For vqsr mismatch out put of mismatch_spliter.sh
miss <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/error_rate/mismatch_rate/vqsr_missing_allele", header=T, sep="/t")
exceed <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/error_rate/mismatch_rate/vqsr_exceeding_allele", header=T, sep="/t")
error <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/error_rate/mismatch_rate/vqsr_calling_error", header=T, sep="/t")
DP.mis <- rbind(miss, exceed, error)
mis.type <-rep("missing",times=nrow(miss))
exceed.type <-rep("exceeding",times=nrow(exceed))
error.type <- rep("error",times=nrow(error))
type <- c(mis.type, exceed.type, error.type)
DP.density <- cbind(DP.mis, type)
ggplot(DP.density, aes(x = DP, fill = type)) + 
  geom_density(alpha = 0.5, bw=0.4) + 
  xlim(0,12) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))



# check final vcf mean coverage (of all sample for each site) distribution
library(ggplot2)
mydat <- read.table("C:/Users/wengz/Box/Chapter 3/Genetics/variant_calling/final_vqsr_DP", header=F, sep="=")
mydat <- transform(mydat, new = V2 / 382)
ggplot(mydat, aes(x=new)) + geom_density() + xlim(c(0, 10))

