args <- commandArgs(trailingOnly=TRUE)
seed <- args[1]
prefix <- args[2]

#chrWide
tracts <- read.table(paste(prefix, ".tsv", sep=""), col.names = c("start","end","indiv","ancestry"))
chrTracts <- cbind(data.frame(chr="ScyDAA6-1196-HRSCAF-1406"), tracts)
tracts$length <- tracts$end - tracts$start
numtracts <- mean(table(tracts$indiv))
length_minor_tracts <- median(aggregate(length ~ indiv, subset(tracts, ancestry == 2), FUN = median)$length)
tracts$weightlen <- (tracts$length * tracts$ancestry/31064592)/2
avganc <- mean(1-aggregate(weightlen ~ indiv, data = tracts, FUN = sum)$weightlen)
varanc <- var(1-aggregate(weightlen ~ indiv, data = tracts, FUN = sum)$weightlen)
cvanc <- sqrt(varanc)/avganc

cat(paste(args[1], numtracts, length_minor_tracts, avganc, varanc, cvanc,"\n"))
###OUTPUT column order###
##simID numberOfTracts medianLengthOfMinorTracts chrWideAvgAnc chrWideVarAnc chrWideCoV
