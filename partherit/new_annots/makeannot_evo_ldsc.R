##Script to make annotation files for LD score regression
## 
library(GenomicRanges);
# library(GenomicFeatures);
library(tools)
options(stringsAsFactors=FALSE);

##Output directory for annotation files
outputdir = "/data/clusterfs/lag/users/gokala/enigma-evol/partherit/new_annot/"

##Get a vector of bed file paths
bedFiles <- Sys.glob("/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/supplemental_table1/beds/*.bed")

##Function
bed_to_granges <- function(bedfile){
  bedfile <- read.table(file = bedfile, header = FALSE)
  # my bed files are all chr# start stop, no headers
  evo_granges <- GRanges(bedfile$V1, IRanges(bedfile$V2, bedfile$V3))
  newseqlevelnames = unique(substring(text = bedfile$V1,first = 4,last = nchar(bedfile$V1)));
  evo_granges <- renameSeqlevels(x = evo_granges,value = newseqlevelnames);
  
  return(evo_granges)
}


a <- bed_to_granges(bedFiles[[1]])

##loop over each of the tissue types
for (i in 1:length(bedFiles)) {

  evo_name <- file_path_sans_ext(basename(bedFiles[[i]]))
  message("*******")
  message(evo_name)
  message("*******")

  evo_granges <- bed_to_granges(bedFiles[[i]])
  print(head(evo_granges))
  
    ##Loop over chromosomes
    for (chr in 1:22) {
    	##Read in locations of SNPs from the annotation files
    	snps.table = read.table(gzfile(paste0("/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/baselineLD_v2.2/baselineLD.",chr,".annot.gz")),sep="\t",header=TRUE);
    	snps = GRanges(snps.table$CHR,IRanges(snps.table$BP,snps.table$BP),SNP=snps.table$SNP,CM=snps.table$CM);
    	## make a blank column to fill with 1s where evo SNP is also in baseline list
    	indicator = matrix(data = 0, nrow = length(snps), ncol = 1);
    	##Find the evo snps that overlap with the baseline
      olap.cat = findOverlaps(snps,evo_granges);     
      ##Replace 0s with 1s where there's overlap
      indicator[queryHits(olap.cat),1] = 1;
      
      print(sum(indicator))

	    ##Write out the new annotation file 
      snps.cat = snps;
      mcols(snps.cat) = cbind(as.data.frame(mcols(snps)),indicator);
      
      print(head(snps.cat))
      
      outframe = data.frame(CHR=as.character(seqnames(snps.cat)),
                            BP=start(snps.cat),
                            SNP=mcols(snps.cat)$SNP,
                            CM=mcols(snps.cat)$CM,
                            indicator=mcols(snps.cat)$indicator);
      
      # outframe = cbind(outframe,as.data.frame(mcols(snps.cat))[,3:(length(categories)+2)]);
      colnames(outframe)[5] = evo_name;
      dir.create(paste0(outputdir,evo_name),showWarnings=FALSE);
      gz1 = gzfile(paste0(outputdir,evo_name,"/",evo_name, ".",chr,".annot.gz"), "w")
      write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
      close(gz1);
    }
}
