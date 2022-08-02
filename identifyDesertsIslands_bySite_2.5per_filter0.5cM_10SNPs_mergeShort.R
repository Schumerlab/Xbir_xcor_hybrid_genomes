options(stringAsFactors=FALSE)
args <- commandArgs(TRUE)
popName <- args[1]

if (file.exists("average_ancestry_allPops_codeConservedSynNonSyn_0.05cM_windows.txt")) {
  print("0.05cM windowed data for all populations exists and full pipeline will run")
} else {
  print("Missing average_ancestry_allPops_codeConservedSynNonSyn_0.05cM_windows.txt will only identify raw deserts and islands, will not do the post filtering")
}

aimData <- read.table(paste("average_ancestry_by_site_", popName, ".txt", sep=""), header=T)
meanHI <- mean(aimData$hybrid_index, na.rm=T)
if (meanHI > 0.5) {
  aimData$majPar <- aimData$hybrid_index
  aimData$minPar <- 1-aimData$hybrid_index
} else {
  aimData$majPar <- 1-aimData$hybrid_index
  aimData$minPar <- aimData$hybrid_index
}

meanMajAnc <- mean(aimData$majPar, na.rm = T)
meanMinAnc <- mean(aimData$minPar, na.rm = T)
quant99 <- quantile(aimData$minPar, 0.99, na.rm = T)
quant97.5 <- quantile(aimData$minPar, 0.975, na.rm = T)
quant95 <- quantile(aimData$minPar, 0.95, na.rm = T)
quant1 <- quantile(aimData$minPar, 0.01, na.rm = T)
quant2.5 <- quantile(aimData$minPar, 0.025, na.rm = T)
quant5 <- quantile(aimData$minPar, 0.05, na.rm = T)

uniChrs <- unique(aimData$group)
minIsleAIMData <- data.frame()
minIsleAIMRegs <- data.frame()
for (chr in uniChrs) {
  chrData <- subset(aimData, group==chr)
  i = 1
  while(i < nrow(chrData)) {
    siteMinPar <- chrData$minPar[i]
    if (siteMinPar >= quant97.5 & !is.na(siteMinPar)) {
      #print(i)
      startRow <- i
      endRow <- i
      startMinPar <- siteMinPar
      endMinPar <- siteMinPar
      chrSur <- chrData[i,]
      while (startMinPar >= quant95 & !is.na(startMinPar)) {
        testRow <- startRow-1
        if (testRow<1) {
          startMinPar <- NA
        } else {
          nextRow <- chrData[testRow,]
          chrSur <- rbind(chrSur, nextRow)
          startMinPar <- nextRow$minPar
          startRow <- testRow
        }
      }
      while (endMinPar >= quant95 & !is.na(endMinPar)) {
        testRow <- endRow+1
        if (testRow > nrow(chrData)) {
          endMinPar <- NA
        } else {
          nextRow <- chrData[testRow,]
          chrSur <- rbind(chrSur, nextRow)
          endMinPar <- nextRow$minPar
          endRow <- testRow
        }
        i <- endRow
      }
      chrSurSort <- chrSur[order(chrSur$position),]
      regStart <- min(chrSurSort$position)
      regEnd <- max(chrSurSort$position)
      regName <- paste(chr, ":", regStart, "_", regEnd, sep="")
      regSum <- data.frame(chrName=chr, start=regStart, end=regEnd, length=regEnd-regStart, midpoint=regStart+((regEnd-regStart)/2), hybrid_index=mean(chrSurSort$hybrid_index, na.rm = T), minPar=mean(chrSurSort$minPar, na.rm = T), inds=mean(chrSurSort$indivs_cov), AIMs=nrow(chrSurSort), regName=regName)
      minIsleAIMRegs <- rbind(minIsleAIMRegs, regSum)
      chrSurSort$regName <- regName
      minIsleAIMData <- rbind(minIsleAIMData, chrSurSort)
    }
    else {
      i=i+1
    }
  }
}

#minIsleAIMRegs$midpoint <- minIsleAIMRegs$start+((minIsleAIMRegs$end-minIsleAIMRegs$start)/2)

write.table(minIsleAIMRegs, paste(popName, "_minorParentIslandsInfo_fromSiteData.txt", sep=""), quote=F, row.names = F, sep="\t")
#write.table(minIsleAIMData, paste(popName, "_minorParentIslands_ancestryBySites.txt", sep=""), quote=F, row.names = F, sep="\t")

minDesertAIMData <- data.frame()
minDesertAIMRegs <- data.frame()
for (chr in uniChrs) {
  #print(chr)
  chrData <- subset(aimData, group==chr)
  i = 1
  while(i < nrow(chrData)) {
    siteMinPar <- chrData$minPar[i]
    if (siteMinPar <= quant2.5 & !is.na(siteMinPar)) {
      #print(i)
      startRow <- i
      endRow <- i
      startMinPar <- siteMinPar
      endMinPar <- siteMinPar
      chrSur <- chrData[i,]
      while (startMinPar <= quant5 & !is.na(startMinPar)) {
        testRow <- startRow-1
        if (testRow<1) {
          startMinPar <- NA
        } else {
          nextRow <- chrData[testRow,]
          chrSur <- rbind(chrSur, nextRow)
          startMinPar <- nextRow$minPar
          startRow <- testRow
        }
      }
      while (endMinPar <= quant5 & !is.na(endMinPar)) {
        testRow <- endRow+1
        if (testRow > nrow(chrData)) {
          endMinPar <- NA
        } else {
          nextRow <- chrData[testRow,]
          chrSur <- rbind(chrSur, nextRow)
          endMinPar <- nextRow$minPar
          endRow <- testRow
        }
        i <- endRow
      }
      chrSurSort <- chrSur[order(chrSur$position),]
      regStart <- min(chrSurSort$position)
      regEnd <- max(chrSurSort$position)
      regName <- paste(chr, ":", regStart, "_", regEnd, sep="")
      regSum <- data.frame(chrName=chr, start=regStart, end=regEnd, length=regEnd-regStart, midpoint=regStart+((regEnd-regStart)/2), hybrid_index=mean(chrSurSort$hybrid_index, na.rm = T), minPar=mean(chrSurSort$minPar, na.rm = T), inds=mean(chrSurSort$indivs_cov), AIMs=nrow(chrSurSort), regName=regName)
      minDesertAIMRegs <- rbind(minDesertAIMRegs, regSum)
      chrSurSort$regName <- regName
      minDesertAIMData <- rbind(minDesertAIMData, chrSurSort)
      #print(i)
    }
    else {
      i=i+1
    }
  }
}

#minDesertAIMRegs$midpoint <- minDesertAIMRegs$start+((minDesertAIMRegs$end-minDesertAIMRegs$start)/2)

write.table(minDesertAIMRegs, paste(popName, "_minorParentDesertsInfo_fromSiteData.txt", sep=""), quote=F, row.names = F, sep="\t")
#write.table(minDesertAIMData, paste(popName, "_minorParentDesertsInfo_ancestryBySite.txt", sep=""), quote=F, row.names = F, sep="\t")

###Here's where the filtering for 0.5cM ancestry, removing regions without 10 SNPs, and merging close by regions starts
cMdata <- read.table("average_ancestry_allPops_codeConservedSynNonSyn_0.05cM_windows.txt", header=T)
if (file.exists("minorParentQuantiles_allPops_0.05cMwins.txt")) {
  print("Already have the quantile file!")
} else {
  print("Getting quantile info")
  popNames <- colnames(cMdata[10:ncol(cMdata)])
  cMquantiles <- data.frame()
  for (pop in popNames) {
    print(pop)
    compCol <- grep(paste("^", pop, "$", sep=""), colnames(cMdata))
    popAnc <- cMdata[,compCol]
    popQuants <- quantile(popAnc, c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99), na.rm=T)
    tmpRow <- data.frame(pop=pop, meanMinPar=mean(popAnc, na.rm = T), lower1=popQuants[[1]], lower2.5=popQuants[[2]], lower5=popQuants[[3]], lower10=popQuants[[4]], upper10=popQuants[[5]], upper5=popQuants[[6]], upper2.5=popQuants[[7]], upper1=popQuants[[8]])
    cMquantiles <- rbind(cMquantiles, tmpRow)
  }
  write.table(cMquantiles, "minorParentQuantiles_allPops_0.05cMwins.txt", row.names = F, quote=F, sep="\t")
}

cMquantiles <- read.table("minorParentQuantiles_allPops_0.05cMwins.txt", header=T)

regTypes <- c("Deserts", "Islands")
###Getting cM windows of interest
#popName <- "CHPL"
#regType <- "Deserts"
for (regType in regTypes) {
  compCol <- grep(paste("^", popName, "$", sep=""), colnames(cMdata))
  if (regType=="Deserts") {
    cutoff <- cMquantiles$lower10[which(cMquantiles$pop==popName)]
  } else if (regType=="Islands") {
    cutoff <- cMquantiles$upper10[which(cMquantiles$pop==popName)]
  }
  
  regionDataRaw <-read.table(paste(popName, "_minorParent", regType, "Info_fromSiteData.txt", sep=""), header=T)
  regionData <- data.frame(chrName=regionDataRaw$chrName, start=regionDataRaw$start, end=regionDataRaw$end, midpoint=regionDataRaw$midpoint, length=regionDataRaw$length, regName=as.character(regionDataRaw$regName)) 
  regionData$cMstart <- NA
  regionData$cMend <- NA
  regionData$cMwin <- NA
  regionData$cMsnps <- NA
  regionData$minAnc <- NA
  regionData$mid10cM <- NA
  
  for (i in 1:nrow(regionData)) {
    chrTest <- regionData$chrName[i]
    midTest <- regionData$midpoint[i]
    cMwin <- subset(cMdata, (as.character(chr)==chrTest & start<=midTest & end>=midTest))
    if (nrow(cMwin)>0) {
      winMinAnc <- cMwin[,compCol]
      regionData$minAnc[i] <- winMinAnc
      if (!is.na(winMinAnc)) {
        if (regType == "Deserts") {
          if (winMinAnc <= cutoff) {
            regionData$mid10cM[i] <- TRUE
          } else if (winMinAnc > cutoff) {
            regionData$mid10cM[i] <- FALSE
          }
        } else if (regType == "Islands") {
          if (winMinAnc >= cutoff) {
            regionData$mid10cM[i] <- TRUE
          } else if (winMinAnc < cutoff) {
            regionData$mid10cM[i] <- FALSE
          }
        }
        regionData$cMstart[i] <- cMwin$start
        regionData$cMend[i] <- cMwin$end
        regionData$cMwin[i] <- paste(chrTest,cMwin$start,cMwin$end,sep="_")
        regionData$cMsnps[i] <- cMwin$nSNPs
      }
    }
  }
  ##10 SNPs as a minimum
  regionTrue <- subset(regionData, mid10cM==TRUE)
  regionToKeep <- data.frame()
  for (cMwinName in unique(regionTrue$cMwin)) {
    cMwinData <- subset(regionTrue, cMwin==cMwinName)
    numSNPs <- cMwinData$cMsnps[1]
    if (numSNPs > 9) {
      if (nrow(cMwinData)>1) {
        newReg <- data.frame(chrName=cMwinData$chrName[1], start=min(cMwinData$start), end=max(cMwinData$end), midpoint=median(cMwinData$midpoint), length=max(cMwinData$end)-min(cMwinData$start), regName=paste(cMwinData$chrName[1], ":", min(cMwinData$start), "_", max(cMwinData$end), sep=""), cMstart=cMwinData$cMstart[1], cMend=cMwinData$cMend[1], cMwin=cMwinData$cMwin[1], cMsnps=cMwinData$cMsnps[1], minAnc=cMwinData$minAnc[1], mid10cM=TRUE)
        regionToKeep <- rbind(regionToKeep, newReg)
      } else {
        regionToKeep <- rbind(regionToKeep, cMwinData)
      }
    }
  }
  write.table(regionToKeep, paste(popName, "_minorParent", regType, "_cMpass.txt", sep=""), row.names=F, quote=F, sep="\t")
  print(popName)
  print(regType)
  print(dim(regionData))
  print(dim(regionTrue))
  print(dim(regionToKeep))
  print(summary(regionToKeep))
}


for (regType in regTypes) {
  compCol <- grep(paste("^", popName, "$", sep=""), colnames(cMdata))
  if (regType=="Deserts") {
    cutoff <- cMquantiles$lower10[which(cMquantiles$pop==popName)]
  } else if (regType=="Islands") {
    cutoff <- cMquantiles$upper10[which(cMquantiles$pop==popName)]
  }
  regions <- read.table(paste(popName, "_minorParent", regType, "_cMpass.txt", sep=""), header=T)
  regions$distToPrev <- NA
  regions$distToNext <- NA
  uniChrs <- unique(regions$chrName)
  for (chr in uniChrs) {
    chrData <- subset(regions, chrName==chr)
    preEnd <- NA
    nextStart <- NA
    for (i in 1:nrow(chrData)) {
      nextStart <- chrData$start[i+1]
      regInfo <- chrData[i,]
      regStart <- regInfo$start
      regEnd <- regInfo$end
      preDist <- regStart-preEnd
      nextDist <- nextStart-regEnd
      regions$distToPrev[which(regions$chrName==chr & regions$start==regStart & regions$end==regEnd)] <- preDist
      regions$distToNext[which(regions$chrName==chr & regions$start==regStart & regions$end==regEnd)] <- nextDist
      preEnd <- regEnd
    }
  }
  
  mergeDF <- data.frame()
  i = 1
  while(i < nrow(regions)) {
    regDist <- regions$distToNext[i]
    j <- i
    if (!is.na(regDist) & regDist < 50000) {
      #print(i)
      mergeRegs <- regions[i,]
      nextRegDist <- regions$distToNext[i+(j-i)]
      while (!is.na(nextRegDist) & nextRegDist < 50000) {
        #print(nextRegDist)
        j <- j + 1
        nextRegDist <- regions$distToNext[i+(j-i)]
        mergeRegs <- rbind(mergeRegs, regions[i+(j-i),])
      }
      #print(paste(as.character(regions$regName[i]), "to merge", (j-i), sep=" "))
      mergeChr <- mergeRegs$chrName[1]
      newStart <- min(mergeRegs$start)
      newEnd <- max(mergeRegs$end)
      newLen <- newEnd-newStart
      newMid <- as.integer(newStart+(newLen/2))
      newRegName <- paste(mergeChr, ":", newStart, "_", newEnd, sep="")
      cMwin <- subset(cMdata, (as.character(chr)==mergeChr & start<=newMid & end>=newMid-1))
      pass <- NA
      winMinAnc <- cMwin[,compCol]
      if (!is.na(winMinAnc)) {
        if (regType == "Deserts") {
          if (winMinAnc <= cutoff) {
            pass <- TRUE
          } else if (winMinAnc > cutoff) {
            pass <- FALSE
          }
        } else if (regType == "Islands") {
          if (winMinAnc >= cutoff) {
            pass <- TRUE
          } else if (winMinAnc < cutoff) {
            pass <- FALSE
          }
        }
      }
      cMstart <- cMwin$start
      cMend <- cMwin$end
      cMregName <- paste(mergeChr, cMstart, cMend, sep="_")
      cMsnps <- cMwin$nSNPs
      newPreDist <- NA
      if (i!=1) {
        if (mergeChr == regions$chrName[i-1]) {
          newPreDist <- newStart - regions$end[i-1]
          if (newPreDist==0) {
            #print(regions[i-1,])
            #print(regions[i,])
          }
        }
      }
      newNextDist <- NA
      if (j+1 <= nrow(regions)) {
        if (mergeChr == regions$chrName[j+1]) {
          newNextDist <- regions$start[j+1] - newEnd
        }
      }
      newRegDF <- data.frame(chrName=mergeChr, start=newStart, end=newEnd, midpoint=newMid, length=newLen, regName=newRegName, cMstart=cMstart, cMend=cMend, cMwin=cMregName, cMsnps=cMsnps, minAnc=winMinAnc, mid10cM=pass, distToPrev=newPreDist, distToNext=newNextDist)
      mergeDF <- rbind(mergeDF, newRegDF)
      #print(i)
      #print(paste("merge end = ", j, sep=""))
      i = j + 1
    } else {
      mergeDF <- rbind(mergeDF, regions[i,])
      unMergeDist <- regions$distToPrev[i]
      if (!is.na(unMergeDist) & unMergeDist==0) {
        print(i)
      }
      i = i+1
    }
  }
  toKeep <- subset(mergeDF, length>10000)
  toKeepReOrd <- toKeep[,c(1:10,13,14,11,12)]
  names(toKeepReOrd)[13] <- paste(popName,"minAnc",sep="")
  names(toKeepReOrd)[14] <- paste(popName,"mid10cM",sep="")
  write.table(toKeepReOrd, paste(popName, "_minorParent", regType, "_cMpass_shortMerged.txt", sep=""), row.names=F, quote=F, sep="\t")
  print(popName)
  print(regType)
  print(dim(regions))
  print(dim(mergeDF))
  print(dim(toKeep))
  print(summary(toKeepReOrd))
}


