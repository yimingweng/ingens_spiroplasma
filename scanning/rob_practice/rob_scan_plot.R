########## YiMing's code

#setwd("C:/Users/wengz/Box/Spiroplasma/scanning")

#coverage <- read.table("all_coverage", header=F, sep="\t")
#hist(coverage$V2, main = "covered_percent", xlab= "percentage", ylab="count")

######

####code copy cat#####

setwd("C:/Users/19207/Box/Spiroplasma/scanning")

coverage <- read.table("all_coverage", header = FALSE, sep = "\t")

str(coverage)
#


head(coverage)
tail(coverage)
length(coverage)
length(coverage$V1)


# ?hist               x = vector of values; main = title
# y axis too short; add height

hist(coverage$V2, main = "covered_percent", xlab = "percentage", ylab = "count")

#add height with ylim

hist(coverage$V2, main = "covered_percent", xlab = "percentage", ylab = "count",  ylim=c(0,250))

#add more bins

hist(coverage$V2, main = "covered_percent", xlab = "percentage", ylab = "count",  ylim=c(0,15), breaks = 1000)


#what is covered percent?