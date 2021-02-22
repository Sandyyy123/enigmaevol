#Script to find overlap between BICseq outputs
#It first defines deletions and duplications and splits them for each sample 
#Then it uses foverlaps to find any overlap
#Then it calculates this overlap and filters hits of foverlap that have a big difference in lenght and filters according to %based overlap (compared to the mean lenght of the two compared CNVs).

.libPaths(c(.libPaths(), "D:/R_packages"))
library(data.table)

SYD12_BICseq <- read.table("P:/lg-ngs/working/SYD_analysis_Ivo/BICseq/CNVs/SYD12/SYD12.CNVprofile_L2.out", header = TRUE)
SYD12_BICseq$length <- (SYD12_BICseq$end - SYD12_BICseq$start)
SYD13_BICseq <- read.table("P:/lg-ngs/working/SYD_analysis_Ivo/BICseq/CNVs/SYD13/SYD13.CNVprofile_L2.out", header = TRUE)
SYD13_BICseq$length <- (SYD13_BICseq$end - SYD13_BICseq$start)

SYD12deletions <- data.table(SYD12_BICseq[which(SYD12_BICseq$log2.copyRatio < -0.5 & SYD12_BICseq$log2.copyRatio > -2), c(1:9)])
SYD12duplications <- data.table(SYD12_BICseq[which(SYD12_BICseq$log2.copyRatio < 2 & SYD12_BICseq$log2.copyRatio > 0.4), c(1:9)])

SYD13deletions <- data.table(SYD13_BICseq[which(SYD13_BICseq$log2.copyRatio < -0.5 & SYD13_BICseq$log2.copyRatio > -2), c(1:9)])
SYD13duplications <- data.table(SYD13_BICseq[which(SYD13_BICseq$log2.copyRatio < 2 & SYD13_BICseq$log2.copyRatio > 0.4), c(1:9)])

write.table(SYD12deletions, "P:/lg-ngs/working/SYD_analysis_Ivo/bedtools/BICseq.deletions/SYD12deletions-lambda2.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE, dec = ".", quote = FALSE, append = FALSE)
write.table(SYD12duplications, "P:/lg-ngs/working/SYD_analysis_Ivo/bedtools/BICseq.duplications/SYD12duplications-lambda2.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE, dec = ".", quote = FALSE, append = FALSE)
write.table(SYD13deletions, "P:/lg-ngs/working/SYD_analysis_Ivo/bedtools/BICseq.deletions/SYD13deletions-lambda2.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE, dec = ".", quote = FALSE, append = FALSE)
write.table(SYD13duplications, "P:/lg-ngs/working/SYD_analysis_Ivo/bedtools/BICseq.duplications/SYD13duplications-lambda2.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE, dec = ".", quote = FALSE, append = FALSE)

setkey(SYD06deletions, chrom, start, end)
setkey(SYD06duplications, chrom, start, end)
deletionoverlap <- foverlaps(SYD01deletions, SYD06deletions, by.x=c("chrom", "start", "end"), type="any", mult="all", which=TRUE, nomatch = 0)
duplicationoverlap <- foverlaps(SYD01duplications, SYD06duplications, by.x=c("chrom", "start", "end"), type="any", mult="all", which=TRUE, nomatch = 0)

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
# combine all regions with low mappability that are close to each other (within a boundary)
# do this per chromosome to save time

#CHR <- "chr10"
CHR <- unique(map$chr)
boundary <- 25

map3 <- c()
for (chr in CHR) {
  map_chr <- map[which(map$chr == chr),]
  Sys.sleep(0.1)
  print(chr)
  flush.console()
  cond <- c(map_chr$pos2[1:nrow(map_chr)-1] > as.numeric(map_chr$pos1[2:nrow(map_chr)]) - boundary, FALSE)
  toRemove <- rep("keep", nrow(map_chr))
  # first remove all rows that are not necessary
  for (i in 2:(length(cond)-1)) {
    if (cond[i-1] && cond[i+1] && cond[i]) {
      toRemove[i] <- "remove"
    }
  }
  map2_chr <- map_chr[which(toRemove == "keep"),]
  cond2 <- cond[which(toRemove == "keep")]
  Sys.sleep(0.1)
  print("first step is ready")
  flush.console()
  #now merge the boundaries
  j <- 1
  toRemove2 <- rep("keep", nrow(map2_chr))
  for (j in nrow(map2_chr):1){
    if (cond2[j]) {
      map2_chr$pos2[j] <- map2_chr$pos2[j+1]
      toRemove2[j+1] <- "remove"
    }              
  }
  map3_chr <- map2_chr[which(toRemove2 == "keep"),]
  map3 <- rbind(map3, map3_chr)
}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

deletions.goodoverlap <- c()
duplications.goodoverlap <- c()
#filter the deletions by % overlap and difference in total lenght
filtered.deletion.overlap <- c()
deletions.combined <- c()
for (i in 1:nrow(deletionoverlap)) {
  SYD01line <- SYD01deletions[as.numeric(deletionoverlap[i,1])]                                #use the data from foverlaps to select the lines that need to be compared.
  SYD06line <- SYD06deletions[as.numeric(deletionoverlap[i,2])]
  combinedlines <- cbind(SYD01line, SYD06line)
  if (((SYD01line$length/SYD06line$length) > 0.3) & ((SYD01line$length/SYD06line$length) < 3 )) {   #filter deletions that have a big difference in length
    #Calculate %based overlap according to average length of the two deletions that are being compared. 
    #There are four differend possible ways of overlap.
    if ((SYD06line$start >= SYD01line$start)&(SYD06line$end <= SYD01line$end)){                #check if SYDxx falls totaly within SYD01
      if ((SYD06line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){  #calculate wether the overlap is more than 50%
        combinedlines$overlap <- (SYD06line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD06line$end - SYD06line$start)
        deletions.goodoverlap <- rbind(deletions.goodoverlap, combinedlines)                                 #if that is the case, save the line in a new data frame.
      } else {print(paste0("omitted due to less than 50% overlap:", deletionoverlap[i]))}
    } else if ((SYD06line$start >= SYD01line$start)&(SYD06line$end > SYD01line$end)){                #check if SYDxx overlaps with the tail of SYD01
      if ((SYD01line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){
        combinedlines$overlap <-(SYD01line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD01line$end - SYD06line$start)
        deletions.goodoverlap <- rbind(deletions.goodoverlap, combinedlines)
      } else {print(paste0("omitted due to less than 50% overlap:", deletionoverlap[i]))}
    } else if ((SYD06line$start < SYD01line$start)&(SYD06line$end <= SYD01line$end)){                #check if SYDxx overlaps with its tail with SYD01
      if ((SYD06line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){
        combinedlines$overlap <- (SYD06line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD06line$end - SYD01line$start)
        deletions.goodoverlap <- rbind(deletions.goodoverlap, combinedlines)
      }else {print(paste0("omitted due to less than 50% overlap:", deletionoverlap[i]))}
    } else if ((SYD06line$start < SYD01line$start)&(SYD06line$end > SYD01line$end)){                #check if SYD01 falls within SYDxx
      if ((SYD01line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){  
        combinedlines$overlap <- (SYD01line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD01line$end - SYD01line$start)
        deletions.goodoverlap <- rbind(deletions.goodoverlap, combinedlines)
      }else {print(paste0("omitted due to less than 0.5 overlap:", deletionoverlap[i]))}
    }
  } else {
    print(paste0("omitted:", deletionoverlap[i]))
  }
}
colnames(deletions.goodoverlap) <- c("SYD01chrom", "SYD01start", "SYD01end", "SYD01log2.copyRatio", "SYD01pvalue", "SYD01length", "SYD06chrom", "SYD06start", "SYD06end", "SYD06log2.copyRatio", "SYD06pvalue", "SYD06length", "overlap", "overlap length")
#print(deletions.goodoverlap)
SYD01vsSYD06deletionoverlap <- as.data.frame(deletions.goodoverlap)
#write.table(SYD01vsSYD06deletionoverlap, "P:/lg-ngs/working/SYD_analysis_Ivo/BICseq/SYD01vsSYD06deletions.txt", sep = " ", na = "NA", row.names = FALSE, col.names = TRUE, dec = ".")

#Now the same for duplications
filtered.duplication.overlap <- c()
for (i in 1:nrow(duplicationoverlap)) {
  SYD01line <- SYD01duplications[as.numeric(duplicationoverlap[i,1])]                                #use the data from foverlaps to select the lines that need to be compared.
  SYD06line <- SYD06duplications[as.numeric(duplicationoverlap[i,2])]
  combinedlines <- cbind(SYD01line, SYD06line)
  if (((SYD01line$length/SYD06line$length) > 0.3) & ((SYD01line$length/SYD06line$length) < 3 )) {   #filter duplications that have a big difference in length
    #Calculate %based overlap according to average length of the two duplications that are being compared. 
    #There are four differend possible ways of overlap.
    #print(duplicationoverlap[i])
    if ((SYD06line$start >= SYD01line$start)&(SYD06line$end <= SYD01line$end)){                #check if SYDxx falls totaly within SYD01
      if ((SYD06line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){  #calculate wether the overlap is more than 50%
        combinedlines$overlap <- (SYD06line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD06line$end - SYD06line$start)
        duplications.goodoverlap <- rbind(duplications.goodoverlap, combinedlines)                                 #if that is the case, save the line in a new data frame.
      } else {print(paste0("omitted due to less than 50% overlap:", duplicationoverlap[i]))}
    } else if ((SYD06line$start >= SYD01line$start)&(SYD06line$end > SYD01line$end)){                #check if SYDxx overlaps with the tail of SYD01
      if ((SYD01line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){
        combinedlines$overlap <-(SYD01line$end - SYD06line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD01line$end - SYD06line$start)
        duplications.goodoverlap <- rbind(duplications.goodoverlap, combinedlines)
      } else {print(paste0("omitted due to less than 50% overlap:", duplicationoverlap[i]))}
    } else if ((SYD06line$start < SYD01line$start)&(SYD06line$end <= SYD01line$end)){                #check if SYDxx overlaps with its tail with SYD01
      if ((SYD06line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){
        combinedlines$overlap <- (SYD06line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD06line$end - SYD01line$start)
        duplications.goodoverlap <- rbind(duplications.goodoverlap, combinedlines)
      }else {print(paste0("omitted due to less than 50% overlap:", duplicationoverlap[i]))}
    } else if ((SYD06line$start < SYD01line$start)&(SYD06line$end > SYD01line$end)){                #check if SYD01 falls within SYDxx
      if ((SYD01line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2) > 0.5){  
        combinedlines$overlap <- (SYD01line$end - SYD01line$start)/((SYD06line$length + SYD01line$length)/2)
        combinedlines$overlaplenght <- (SYD01line$end - SYD01line$start)
        duplications.goodoverlap <- rbind(duplications.goodoverlap, combinedlines)
      }else {print(paste0("omitted due to less than 0.5 overlap:", duplicationoverlap[i]))}
    }
  } else {
    print(paste0("omitted:", duplicationoverlap[i]))
  }
}
colnames(duplications.goodoverlap) <- c("SYD01chrom", "SYD01start", "SYD01end", "SYD01log2.copyRatio", "SYD01pvalue", "SYD01length", "SYD06chrom", "SYD06start", "SYD06end", "SYD06log2.copyRatio", "SYD06pvalue", "SYD06length", "overlap", "overlap length")
SYD01vsSYD06duplicationoverlap <- as.data.frame(duplications.goodoverlap)
#write.table(SYD01vsSYD06duplicationoverlap, "P:/lg-ngs/working/SYD_analysis_Ivo/BICseq/SYD01vsSYD06duplications.txt", sep = " ", na = "NA", row.names = FALSE, col.names = TRUE, dec = ".")
#print(duplications.goodoverlap)
