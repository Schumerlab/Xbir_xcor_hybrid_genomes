options(stringAsFactors=FALSE)
corBirPops <- c("STAC", "HUEX")
regTypes <- c("Deserts", "Islands")
allPops <- c("STAC", "HUEX", "TLMC", "ACUA", "AGZC")

cMdata <- read.table("average_ancestry_allPops_codeConservedSynNonSyn_0.05cM_windows.txt", header=T)
cMdata$cMid <- paste(cMdata$chr, cMdata$start, cMdata$end, sep="_")
tlmc <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata$TLMC)
acua <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata$ACUA)
agzc <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata$AGZC)
stac <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata$STAC)
huex <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata$HUEX)

###STAC###
popName <- "STAC"
print(popName)
##Deserts##
regType <- "deserts"
print(regType)
popRegCombined <- data.frame()
regData <- read.table("STAC_crossTypeSharedDeserts.txt", header = T)
#crossRegs <- subset(regData, crossPops=="TRUE")
regRows <- cMdata[which(cMdata$cMid %in% regData$cMwin),]
focal<-regRows
#start simple null simulations
null_positives<-{}
#define the threshold for ancestry outlier
quant=0.05
#3 population requirement, can modify with | for two population requirement
for(x in 1:1000){
  tlmc_null<-tlmc
  acua_null<-acua
  agzc_null<-agzc
  tlmc_null$minParAnc<-sample(tlmc$minParAnc)
  acua_null$minParAnc<-sample(acua$minParAnc)
  agzc_null$minParAnc<-sample(agzc$minParAnc)
  pass<-subset(tlmc_null,tlmc_null$minParAnc<quantile(tlmc_null$minParAnc,quant,na.rm=TRUE) | acua_null$minParAnc<quantile(acua_null$minParAnc,quant,na.rm=TRUE) | agzc_null$minParAnc<quantile(agzc_null$minParAnc,quant,na.rm=TRUE))
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
crossDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop="crossPop", permutations=null_positives)
print(summary(crossDF))
popRegCombined <- data.frame()
for (otherPop in allPops) {
  if (otherPop!=popName) {
    popAncCol <- grep(otherPop, colnames(cMdata))
    testData <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata[,popAncCol])
    popShareCol <- grep(paste(otherPop, "mid10cM", sep=""), colnames(regData))
    crossRegs <- regData[which(regData[,popShareCol]=="TRUE"),]
    regRows <- cMdata[which(cMdata$cMid %in% crossRegs$cMwin),]
    focal<-regRows
    null_positives<-{}
    #define the threshold for ancestry outlier
    quant=0.05
    #3 population requirement, can modify with | for two population requirement
    for(x in 1:1000){
      pop_null<-testData
      pop_null$minParAnc<-sample(testData$minParAnc)
      pass<-subset(pop_null,pop_null$minParAnc<=quantile(pop_null$minParAnc,quant,na.rm=TRUE))
      permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
      null_positives<-c(null_positives,permuted_pos)
    }
    popDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop=otherPop, permutations=null_positives)
    print(otherPop)
    print(summary(popDF))
    popRegCombined <- rbind(popRegCombined, popDF)
  }
}
#write.table(null_positives, paste(popName, "_shared", regType, "_permutations.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
combined <- rbind(crossDF, popRegCombined)
write.table(combined, "STAC_shared_deserts_permutations.txt", quote=F, row.names=F, sep="\t")


###STAC###
popName <- "STAC"
print(popName)
##Deserts##
regType <- "islands"
print(regType)
popRegCombined <- data.frame()
regData <- read.table("STAC_crossTypeSharedIslands.txt", header = T)
#crossRegs <- subset(regData, crossPops=="TRUE")
regRows <- cMdata[which(cMdata$cMid %in% regData$cMwin),]
focal<-regRows
#start simple null simulations
null_positives<-{}
#define the threshold for ancestry outlier
quant=0.05
#3 population requirement, can modify with | for two population requirement
for(x in 1:1000){
  tlmc_null<-tlmc
  acua_null<-acua
  agzc_null<-agzc
  tlmc_null$minParAnc<-sample(tlmc$minParAnc)
  acua_null$minParAnc<-sample(acua$minParAnc)
  agzc_null$minParAnc<-sample(agzc$minParAnc)
  pass<-subset(tlmc_null,tlmc_null$minParAnc<quantile(tlmc_null$minParAnc,quant,na.rm=TRUE) | acua_null$minParAnc<quantile(acua_null$minParAnc,quant,na.rm=TRUE) | agzc_null$minParAnc<quantile(agzc_null$minParAnc,quant,na.rm=TRUE))
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
crossDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop="crossPop", permutations=null_positives)
print(summary(crossDF))
popRegCombined <- data.frame()
for (otherPop in allPops) {
  if (otherPop!=popName) {
    popAncCol <- grep(otherPop, colnames(cMdata))
    testData <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata[,popAncCol])
    popShareCol <- grep(paste(otherPop, "mid10cM", sep=""), colnames(regData))
    crossRegs <- regData[which(regData[,popShareCol]=="TRUE"),]
    regRows <- cMdata[which(cMdata$cMid %in% crossRegs$cMwin),]
    focal<-regRows
    null_positives<-{}
    #define the threshold for ancestry outlier
    quant=0.05
    #3 population requirement, can modify with | for two population requirement
    for(x in 1:1000){
      pop_null<-testData
      pop_null$minParAnc<-sample(testData$minParAnc)
      pass<-subset(pop_null,pop_null$minParAnc<=quantile(pop_null$minParAnc,quant,na.rm=TRUE))
      permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
      null_positives<-c(null_positives,permuted_pos)
    }
    popDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop=otherPop, permutations=null_positives)
    popRegCombined <- rbind(popRegCombined, popDF)
    print(otherPop)
    print(summary(popDF))
  }
}
#write.table(null_positives, paste(popName, "_shared", regType, "_permutations.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
combined <- rbind(crossDF, popRegCombined)
write.table(combined, "STAC_shared_islands_permutations.txt", quote=F, row.names=F, sep="\t")

###HUEX###
popName <- "HUEX"
print(popName)
##Deserts##
regType <- "deserts"
print(regType)
popRegCombined <- data.frame()
regData <- read.table("HUEX_crossTypeSharedDeserts.txt", header = T)
#crossRegs <- subset(regData, crossPops=="TRUE")
regRows <- cMdata[which(cMdata$cMid %in% regData$cMwin),]
focal<-regRows
#start simple null simulations
null_positives<-{}
#define the threshold for ancestry outlier
quant=0.05
#3 population requirement, can modify with | for two population requirement
for(x in 1:1000){
  tlmc_null<-tlmc
  acua_null<-acua
  agzc_null<-agzc
  tlmc_null$minParAnc<-sample(tlmc$minParAnc)
  acua_null$minParAnc<-sample(acua$minParAnc)
  agzc_null$minParAnc<-sample(agzc$minParAnc)
  pass<-subset(tlmc_null,tlmc_null$minParAnc<quantile(tlmc_null$minParAnc,quant,na.rm=TRUE) | acua_null$minParAnc<quantile(acua_null$minParAnc,quant,na.rm=TRUE) | agzc_null$minParAnc<quantile(agzc_null$minParAnc,quant,na.rm=TRUE))
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
crossDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop="crossPop", permutations=null_positives)
print(summary(crossDF))
popRegCombined <- data.frame()
for (otherPop in allPops) {
  if (otherPop!=popName) {
    popAncCol <- grep(otherPop, colnames(cMdata))
    testData <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata[,popAncCol])
    popShareCol <- grep(paste(otherPop, "mid10cM", sep=""), colnames(regData))
    crossRegs <- regData[which(regData[,popShareCol]=="TRUE"),]
    regRows <- cMdata[which(cMdata$cMid %in% crossRegs$cMwin),]
    focal<-regRows
    null_positives<-{}
    #define the threshold for ancestry outlier
    quant=0.05
    #3 population requirement, can modify with | for two population requirement
    for(x in 1:1000){
      pop_null<-testData
      pop_null$minParAnc<-sample(testData$minParAnc)
      pass<-subset(pop_null,pop_null$minParAnc<=quantile(pop_null$minParAnc,quant,na.rm=TRUE))
      permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
      null_positives<-c(null_positives,permuted_pos)
    }
    popDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop=otherPop, permutations=null_positives)
    popRegCombined <- rbind(popRegCombined, popDF)
    print(otherPop)
    print(summary(popDF))
  }
}
#write.table(null_positives, paste(popName, "_shared", regType, "_permutations.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
combined <- rbind(crossDF, popRegCombined)
write.table(combined, "HUEX_shared_deserts_permutations.txt", quote=F, row.names=F, sep="\t")

###HUEX###
popName <- "HUEX"
print(popName)
##Deserts##
regType <- "islands"
print(regType)
popRegCombined <- data.frame()
regData <- read.table("HUEX_crossTypeSharedIslands.txt", header = T)
#crossRegs <- subset(regData, crossPops=="TRUE")
regRows <- cMdata[which(cMdata$cMid %in% regData$cMwin),]
focal<-regRows
#start simple null simulations
null_positives<-{}
#define the threshold for ancestry outlier
quant=0.05
#3 population requirement, can modify with | for two population requirement
for(x in 1:1000){
  tlmc_null<-tlmc
  acua_null<-acua
  agzc_null<-agzc
  tlmc_null$minParAnc<-sample(tlmc$minParAnc)
  acua_null$minParAnc<-sample(acua$minParAnc)
  agzc_null$minParAnc<-sample(agzc$minParAnc)
  pass<-subset(tlmc_null,tlmc_null$minParAnc<quantile(tlmc_null$minParAnc,quant,na.rm=TRUE) | acua_null$minParAnc<quantile(acua_null$minParAnc,quant,na.rm=TRUE) | agzc_null$minParAnc<quantile(agzc_null$minParAnc,quant,na.rm=TRUE))
  permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
  null_positives<-c(null_positives,permuted_pos)
}
crossDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop="crossPop", permutations=null_positives)
print(summary(crossDF))
popRegCombined <- data.frame()
for (otherPop in allPops) {
  if (otherPop!=popName) {
    popAncCol <- grep(otherPop, colnames(cMdata))
    testData <- data.frame(cMchrPos=cMdata$cMid, minParAnc=cMdata[,popAncCol])
    popShareCol <- grep(paste(otherPop, "mid10cM", sep=""), colnames(regData))
    crossRegs <- regData[which(regData[,popShareCol]=="TRUE"),]
    regRows <- cMdata[which(cMdata$cMid %in% crossRegs$cMwin),]
    focal<-regRows
    null_positives<-{}
    #define the threshold for ancestry outlier
    quant=0.05
    #3 population requirement, can modify with | for two population requirement
    for(x in 1:1000){
      pop_null<-testData
      pop_null$minParAnc<-sample(testData$minParAnc)
      pass<-subset(pop_null,pop_null$minParAnc<=quantile(pop_null$minParAnc,quant,na.rm=TRUE))
      permuted_pos<-length(intersect(rownames(pass),rownames(focal)))
      null_positives<-c(null_positives,permuted_pos)
    }
    popDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, compPop=otherPop, permutations=null_positives)
    popRegCombined <- rbind(popRegCombined, popDF)
    print(otherPop)
    print(summary(popDF))
  }
}
#write.table(null_positives, paste(popName, "_shared", regType, "_permutations.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
combined <- rbind(crossDF, popRegCombined)
write.table(combined, "HUEX_shared_islands_permutations.txt", quote=F, row.names=F, sep="\t")


###Jackknife
chunkSize <- 200

uniChrs <- unique(cMdata$chr)

popName <- "STAC"
regTypes <- c("Deserts", "Islands")

for (regType in regTypes) {
  #STAC_deserts_crossPopsShared_Jun21.txt
  crossData <- read.table(paste(popName, "_crossTypeShared", regType, ".txt", sep=""), header = T)
  #totalSharedCount <- nrow(regData)
  #pairSharedCount <- sum(regData$inTypeShared, na.rm=T)
  #crossSharedCount <- nrow(subset(regData,  sharedPops!="none" & sharedPops!="Xcor/Xbir"))
  crossSharedCount <- nrow(crossData)
  
  jackknifeDF <- data.frame()
  for (chrName in uniChrs) {
    chrWins <- subset(cMdata, chr==chrName)
    chrStartRow <- as.numeric(row.names(chrWins)[1])
    chrEndRow <- as.numeric(row.names(chrWins)[nrow(chrWins)])
    start <- 1
    while(start<nrow(chrWins)) {
      chunkData <- chrWins[start:(start+chunkSize),]
      sharedOutside <- crossData[which(!(crossData$cMwin %in% chunkData$cMid)),]
      outsiteTotal <- nrow(sharedOutside)
      #outsidePair <- sum(sharedOutside$pairPops)
      #outsideCross <- nrow(subset(sharedOutside,  crossShared==TRUE))
      tempDF <- data.frame(popComp=paste(popName, regType, sep="_"), focalPop=popName, regType=regType, winName=chunkData$cMid[1], crossPopTotal=outsiteTotal)
      tempCols <- colnames(tempDF)
      for (otherPop in allPops) {
        popShareCol <- grep(paste(otherPop, "mid10cM", sep=""), colnames(crossData))
        outsidePop <- sum(sharedOutside[,popShareCol], na.rm=T)
        tempDF <- cbind(tempDF, count=outsidePop)
        tempCols <- c(tempCols, paste(otherPop, "total", sep=""))
      }
      colnames(tempDF) <- tempCols
      jackknifeDF <- rbind(jackknifeDF, tempDF)
      start = start+chunkSize
    }
  }
  #allPopsJack <- rbind(allPopsJack, jackknifeDF)
  #write.table(jackknifeDF,paste(popName, "_shared", regType, "_jackknife.txt", sep=""), row.names = F, quote=F, sep="\t")
  print(popName)
  print(regType)
  print(summary(jackknifeDF))
  write.table(jackknifeDF, paste(popName, "_", regType, "jackknife.txt", sep=""), row.names = F, quote=F, sep="\t")
}

allPopsJack <- data.frame()
#for (popName in corBirPops) {
  
#}
write.table(allPopsJack, "minorParent_cM_shared_jackknife.txt", row.names = F, quote=F, sep="\t")


#options(stringAsFactors=FALSE)
#corBirPops <- c("STAC", "HUEX")
#regTypes <- c("Deserts", "Islands")

#allPopsJack <- data.frame()
#for (popName in corBirPops) {
#  for (regType in regTypes) {
#    jackData <- read.table(paste(popName, "_shared", regType, "_jackknife.txt", sep=""), header=T)
#    jackData$popComp <- paste(popName, regType, sep="_")
#    allPopsJack <- rbind(allPopsJack, jackData)
#  }
#}

#write.table(allPopsJack, "minorParent_cM_shared_jackknife.txt", row.names = F, quote=F, sep="\t")
