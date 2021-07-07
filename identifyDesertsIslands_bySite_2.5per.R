options(stringAsFactors=FALSE)
args <- commandArgs(TRUE)
popName <- args[1]

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



