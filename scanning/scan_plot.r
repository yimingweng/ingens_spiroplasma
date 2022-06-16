
setwd("C:/Users/wengz/Box/Spiroplasma/scanning")

coverage <- read.table("all_coverage", header=F, sep="\t")
hist(coverage$V2, main = "covered_percent", xlab= "percentage", ylab="count")

