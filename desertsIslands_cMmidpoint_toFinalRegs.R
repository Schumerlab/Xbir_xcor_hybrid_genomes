options(stringAsFactors=FALSE)

popNames <- c("STAC", "HUEX", "TLMC", "ACUA", "AGZC")
allPops <- c("STAC", "HUEX", "TLMC", "ACUA", "AGZC")
#corBirPops <- c("STAC", "HUEX")

cMdata <- read.table("average_ancestry_allPops_codeConservedSynNonSyn_0.05cM_windows.txt", header=T)

cMquantiles <- data.frame()
for (pop in popNames) {
  compCol <- grep(pop, colnames(cMdata))
  popAnc <- cMdata[,compCol]
  popQuants <- quantile(popAnc, c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99), na.rm=T)
  tmpRow <- data.frame(pop=pop, meanMinPar=mean(popAnc, na.rm = T), lower1=popQuants[[1]], lower2.5=popQuants[[2]], lower5=popQuants[[3]], lower10=popQuants[[4]], upper10=popQuants[[5]], upper5=popQuants[[6]], upper2.5=popQuants[[7]], upper1=popQuants[[8]])
  cMquantiles <- rbind(cMquantiles, tmpRow)
}
write.table(cMquantiles, "minorParentQuantiles_allPops_0.05cMwins.txt", row.names = F, quote=F, sep="\t")

regTypes <- c("Deserts", "Islands")
###Getting cM windows of interest
#popName <- "STAC"
#regType <- "Deserts"
for (popName in popNames) {
  for (regType in regTypes) {
    compCol <- grep(popName, colnames(cMdata))
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
          regionData$cMsnps[i] <- cMwin$SNPs
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
}

#popName <- "STAC"
#regType <- "Deserts"
for (popName in popNames) {
  for (regType in regTypes) {
    compCol <- grep(popName, colnames(cMdata))
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
        newMid <- newStart+newLen/2
        newRegName <- paste(mergeChr, ":", newStart, "_", newEnd, sep="")
        cMwin <- subset(cMdata, (as.character(chr)==mergeChr & start<=newMid & end>=newMid))
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
        cMsnps <- cMwin$SNPs
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
}

##And getting it from other pops
#otherPopsSTAC <- c("HUEX", "TLMC", "ACUA", "AGZC")
#otherPopsHUEX <- c("STAC", "TLMC", "ACUA", "AGZC")

for (popName in popNames) {
  for (regType in regTypes) {
    regionData <- read.table(paste(popName, "_minorParent", regType, "_cMpass_shortMerged.txt", sep=""), header=T)
    regionCols <- colnames(regionData)
    nextPopList <- popNames[!popNames %in% popName]
    for (nextPop in nextPopList) {
      regionCols <- c(regionCols, paste(nextPop,"minAnc",sep=""),paste(nextPop,"mid10cM", sep=""))
      regionData$nextMin <- NA
      regionData$next10cM <- NA
      nextPopCol <- grep(nextPop, colnames(cMdata))
      if (regType=="Deserts") {
        nextPopCutoff <- cMquantiles$lower10[which(cMquantiles$pop==nextPop)]
      } else if (regType=="Islands") {
        nextPopCutoff <- cMquantiles$upper10[which(cMquantiles$pop==nextPop)]
      }
      for (i in 1:nrow(regionData)) {
        chrTest <- regionData$chrName[i]
        startTest <- regionData$cMstart[i]
        cMwin <- subset(cMdata, (as.character(chr)==chrTest & start==startTest))
        winMinAnc <- cMwin[,nextPopCol]
        regionData$nextMin[i] <- winMinAnc
        if (!is.na(winMinAnc)) {
          if (regType == "Deserts") {
            if (winMinAnc <= nextPopCutoff) {
              regionData$next10cM[i] <- TRUE
            } else if (winMinAnc > nextPopCutoff) {
              regionData$next10cM[i] <- FALSE
            }
          } else if (regType == "Islands") {
            if (winMinAnc >= nextPopCutoff) {
              regionData$next10cM[i] <- TRUE
            } else if (winMinAnc < nextPopCutoff) {
              regionData$next10cM[i] <- FALSE
            }
          }
        }
      }
      colnames(regionData) <- regionCols
    }
    regionData$totalPops <- NA
    regionData$inTypeShared <- NA
    regionData$crossShared <- NA
    for(i in 1:nrow(regionData)) {
      regData <- regionData[i,]
      regionData$totalPops[i] <- sum(regData$STACmid10cM, regData$HUEXmid10cM, regData$TLMCmid10cM, regData$ACUAmid10cM, regData$AGZCmid10cM, na.rm = TRUE)
      if (popName=="ACUA" | popName=="TLMC" | popName=="AGZC") {
        if(sum(regData$TLMCmid10cM, regData$ACUAmid10cM, regData$AGZCmid10cM, na.rm = TRUE)>1){
          regionData$inTypeShared[i] <- TRUE
        } else {
          regionData$inTypeShared[i] <- FALSE
        }
        if ((regionData$totalPops[i]-sum(regData$TLMCmid10cM, regData$ACUAmid10cM, regData$AGZCmid10cM, na.rm = TRUE))>1) {
          regionData$crossShared[i] <- TRUE
        } else {
          regionData$crossShared[i] <- FALSE
        }
      } else {
        if(sum(regData$STACmid10cM, regData$HUEXmid10cM, na.rm = TRUE)>1){
          regionData$inTypeShared[i] <- TRUE
        } else {
          regionData$inTypeShared[i] <- FALSE
        }
        if ((regionData$totalPops[i]-sum(regData$STACmid10cM, regData$HUEXmid10cM, na.rm = TRUE))>0) {
          regionData$crossShared[i] <- TRUE
        } else {
          regionData$crossShared[i] <- FALSE
        }
      }
    }
    write.table(regionData, paste(popName, "_minorParent", regType, "_allPops_cMmidpoint.txt", sep=""), row.names=F, quote=F, sep="\t")
    print(popName)
    print(regType)
    print(dim(regionData))
    print(sum(regionData$inTypeShared, na.rm = T))
    print(sum(regionData$crossShared, na.rm = T))
    print(max(regionData$totalPops))
  }
}


cMdata <- read.table("average_ancestry_thinned_allPops_codeConservedSynNonSyn_0.05cM_windows.txt", header=T)
cMquantiles <- read.table("minorParentQuantiles_allPops_thinned_0.05cMwins.txt", header=T)

#popName <- "STAC"
#regType <- "Deserts"
for (popName in allPops) {
  for (regType in regTypes) {
    print(popName)
    print(regType)
    regionData <- read.table(paste(popName, "_minorParent", regType, "_allPops_cMmidpoint.txt", sep=""), header=T)
    for (nextPop in allPops) {
      regionCols <- colnames(regionData)
      regionCols <- c(regionCols, paste(nextPop,"minAncThinned",sep=""),paste(nextPop,"mid10cMthinned", sep=""))
      regionData$nextMin <- NA
      regionData$next10cM <- NA
      nextPopCol <- grep(nextPop, colnames(cMdata))
      if (regType=="Deserts") {
        nextPopCutoff <- cMquantiles$lower10[which(cMquantiles$pop==nextPop)]
      } else if (regType=="Islands") {
        nextPopCutoff <- cMquantiles$upper10[which(cMquantiles$pop==nextPop)]
      }
      for (i in 1:nrow(regionData)) {
        chrTest <- regionData$chrName[i]
        startTest <- regionData$cMstart[i]
        cMwin <- subset(cMdata, (as.character(chr)==chrTest & start==startTest))
        winMinAnc <- cMwin[,nextPopCol]
        regionData$nextMin[i] <- winMinAnc
        if (!is.na(winMinAnc)) {
          if (regType == "Deserts") {
            if (winMinAnc <= nextPopCutoff) {
              regionData$next10cM[i] <- TRUE
            } else if (winMinAnc > nextPopCutoff) {
              regionData$next10cM[i] <- FALSE
            }
          } else if (regType == "Islands") {
            if (winMinAnc >= nextPopCutoff) {
              regionData$next10cM[i] <- TRUE
            } else if (winMinAnc < nextPopCutoff) {
              regionData$next10cM[i] <- FALSE
            }
          }
        }
      }
      colnames(regionData) <- regionCols
    }
    regionData$totalPopsThinned <- NA
    regionData$pairPopsThinned <- NA
    regionData$crossPopsThinned <- NA
    for(i in 1:nrow(regionData)) {
      regData <- regionData[i,]
      regionData$totalPopsThinned[i] <- sum(regData$STACmid10cMthinned, regData$HUEXmid10cMthinned, regData$TLMCmid10cMthinned, regData$ACUAmid10cMthinned, regData$AGZCmid10cMthinned, na.rm = TRUE)
      if (!is.na(regData$HUEXmid10cMthinned) & !is.na(regData$STACmid10cMthinned)) {
        if (regData$STACmid10cMthinned==TRUE & regData$HUEXmid10cMthinned==TRUE) {
          regionData$pairPopsThinned[i] <- TRUE
        } else {
          regionData$pairPopsThinned[i] <- FALSE
        }
        if (regionData$pairPopsThinned[i]==TRUE) {
          if (isTRUE(regData$TLMCmid10cMthinned) | isTRUE(regData$ACUAmid10cMthinned) | isTRUE(regData$AGZCmid10cMthinned)) {
            regionData$crossPopsThinned[i] <- TRUE 
          } else {
            regionData$crossPopsThinned[i] <- FALSE
          }
        }
      }
    }
    for (nextPop in allPops) {
      popMean <- cMquantiles$meanMinPar[which(cMquantiles$pop==nextPop)]
      regionCols <- colnames(regionData)
      unColName <- paste(nextPop, "minAnc", sep="")
      thColName <- paste(nextPop, "minAncThinned", sep="")
      popUn <- min(grep(unColName, colnames(regionData)))
      popTh <- grep(thColName, colnames(regionData))
      regionCols <- c(regionCols, paste(nextPop,"minAncDiff",sep=""), paste(nextPop,"minAncComp",sep=""))
      regionData$popDiff <- regionData[,popUn]-regionData[,popTh]
      regionData$popComp <- regionData$popDiff/popMean
      print(nextPop)
      print(summary(regionData$popComp))
      colnames(regionData) <- regionCols
    }
    write.table(regionData, paste(popName, "_minorParent", regType, "_allPops_thinned_cMmidpoint.txt", sep=""), row.names=F, quote=F, sep="\t")
    print(dim(regionData))
    print(sum(regionData$pairPopsThinned, na.rm = T))
    print(sum(regionData$crossPopsThinned, na.rm = T))
    print(max(regionData$totalPopsThinned))
  }
}

for (popName in allPops) {
  for (regType in regTypes) {
    print(popName)
    print(regType)
    regionData <- read.table(paste(popName, "_minorParent", regType, "_allPops_thinned_cMmidpoint.txt", sep=""), header=T)
    nonThinCol <- min(grep(paste(popName, "minAnc", sep=""), colnames(regionData)))
    thinCol <- grep(paste(popName, "minAncThinned", sep=""), colnames(regionData))
    regionData$popPerDiff <- (regionData[,nonThinCol]-regionData[,thinCol])/((regionData[,nonThinCol]+regionData[,thinCol])/2)
    if (regType == "Deserts"){
      deserts <- regionData[,1:25]
      write.table(deserts, paste(popName, "_all", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
      inType <- subset(deserts, inTypeShared=="TRUE")
      write.table(inType, paste(popName, "_inTypeShared", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
      crossType <- subset(deserts, crossShared=="TRUE")
      write.table(crossType, paste(popName, "_crossTypeShared", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
    } else if (regType == "Islands") {
      thinnedIslands <- subset(regionData, popPerDiff<=0.1 & popPerDiff>=-0.1)
      islands <- thinnedIslands[,1:25]
      write.table(islands, paste(popName, "_all", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
      inType <- subset(islands, inTypeShared=="TRUE")
      write.table(inType, paste(popName, "_inTypeShared", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
      crossType <- subset(islands, crossShared=="TRUE")
      write.table(crossType, paste(popName, "_crossTypeShared", regType, ".txt", sep=""), row.names = F, quote = F, sep="\t")
    }
  }
}


# for (popName in popNames) {
#   for (regType in regTypes) {
#     regionData <- read.table(paste(popName, "_minorParent", regType, "_allPops_cMmidpoint.txt", sep=""), header=T)
#     if (popName=="STAC" | popName=="HUEX"){
#       shared <- subset(regionData, inTypeShared==TRUE & crossShared==TRUE)
#     } else {
#       shared <- subset(regionData, crossShared==TRUE)
#     }
#     write.table(shared, paste(popName, "_minorParent", regType, "_allPops_shared.txt", sep=""), row.names=F, quote=F, sep="\t")
#     print(popName)
#     print(regType)
#     print(dim(regionData))
#     print(dim(shared))
#     print(summary(shared))
#   }
# }

