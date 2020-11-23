# Jason's ancestry regressed summary statistics are in Rdata format
# and I need to move them to .txt files so that I can run the LDSC
# munge_sumstats.py script on them, prior to actually running LDSC 
# or partitioned heritability.

# ran on R 3.3.3, using GenomicRanges for bioconductor v3.4
# note: ran as an R session within the USC INI server

# These files do NOT have genomic control correction applied. They were made by Jason on October 25, 2019

library(tools)
library(GenomicRanges)

files <-  as.vector(read.delim("AncestryRegressionData_noGC_Rdata_files.txt", header = FALSE, stringsAsFactors = FALSE))
# files <- files[[1]][grepl(pattern = "Mean_", x = files)] #so we don't include Jason's list of the files.


for (i in 1:length(files[[1]])) {
  fname = file_path_sans_ext(basename(files[[1]][i]))
  if (file.exists(paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/wo_1KGP3_ancreg/",fname,"_1KGP3_nonancreg_noGC.txt"))) {
    message("File exists, moving on.")
    } else {
      load(file = files[[1]][i])
      message(paste("Working on file: ",fname))

      sumstats <- data.frame(CHR = seqnames(mergedGR),
                             BP = start(mergedGR),
                             SNP = mergedGR$SNP,
                             A1 = mergedGR$A1.x,
                             A2 = mergedGR$A2.x,
                             FREQ1 = mergedGR$FREQ1,
                             SE = mergedGR$ancSE,
                             P = mergedGR$ancP,
                             N = mergedGR$N,
                             MARKER = mergedGR$MARKER,
                             BETA = mergedGR$ancBETA)
      write.table(sumstats, file = paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/w_1KGP3_ancreg/",fname,"_1KGP3_ancreg_noGC.txt"), 
                  quote = F, sep = "\t", row.names = F, col.names = T)

      no_ancreg_sumstats <- data.frame(CHR = seqnames(mergedGR),
                             BP = start(mergedGR),
                             SNP = mergedGR$SNP,
                             A1 = mergedGR$A1.x,
                             A2 = mergedGR$A2.x,
                             FREQ1 = mergedGR$FREQ1,
                             SE = mergedGR$SE,
                             P = mergedGR$P,
                             N = mergedGR$N,
                             MARKER = mergedGR$MARKER,
                             BETA = mergedGR$BETA)
      write.table(no_ancreg_sumstats, file = paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/wo_1KGP3_ancreg/",fname,"_1KGP3_nonancreg_noGC.txt"), 
                  quote = F, sep = "\t", row.names = F, col.names = T)
      }
    }

# For Full Surface Area, where we have to scale the betas
files2 <- files[grepl(pattern = "Full_Surf", x = files[[1]])] #so we don't include Height

for (i in 1:length(files2[[1]])) {
  load(file = files2[[1]][i])
  fname = file_path_sans_ext(basename(files2[[1]][i]))
  message(paste("Working on the scaled version of file: ",fname))
  sumstats <- data.frame(CHR = seqnames(mergedGR),
                         BP = start(mergedGR),
                         SNP = mergedGR$SNP,
                         A1 = mergedGR$A1.x,
                         A2 = mergedGR$A2.x,
                         FREQ1 = mergedGR$FREQ1,
                         SE = mergedGR$ancSE/1000,
                         P = mergedGR$ancP,
                         N = mergedGR$N,
                         MARKER = mergedGR$MARKER,
                         BETA = mergedGR$ancBETA/1000)
  write.table(sumstats, file = paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/w_1KGP3_ancreg/",fname,"_1KGP3_ancreg_noGC.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = T)
  no_ancreg_sumstats <- data.frame(CHR = seqnames(mergedGR),
                                   BP = start(mergedGR),
                                   SNP = mergedGR$SNP,
                                   A1 = mergedGR$A1.x,
                                   A2 = mergedGR$A2.x,
                                   FREQ1 = mergedGR$FREQ1,
                                   SE = mergedGR$SE/1000,
                                   P = mergedGR$P,
                                   N = mergedGR$N,
                                   MARKER = mergedGR$MARKER,
                                   BETA = mergedGR$BETA/1000)
  write.table(no_ancreg_sumstats, file = paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/wo_1KGP3_ancreg/",fname,"_1KGP3_nonancreg_noGC.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = T)
}

# Make lists of files for use with LDSC script
ancreg_noGC_files <- paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/w_1KGP3_ancreg/",file_path_sans_ext(basename(files[[1]])),"_1KGP3_ancreg_noGC.txt")
write(ancreg_noGC_files, file = "/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/w_1KGP3_ancreg/1KGPhase3_ancreg_noGC_premunge_sumstats_filenames.txt")

nonancreg_noGC_files <- paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/wo_1KGP3_ancreg/",file_path_sans_ext(basename(files[[1]])),"_1KGP3_nonancreg_noGC.txt")
write(nonancreg_noGC_files, file = "/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/premunge_sumstats/wo_1KGP3_ancreg/1KGPhase3_nonancreg_noGC_premunge_sumstats_filenames.txt")
