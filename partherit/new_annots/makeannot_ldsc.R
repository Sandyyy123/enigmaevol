##Script to make annotation file for LD score regression
library(GenomicRanges);
library(GenomicFeatures);
library(readxl);
options(stringsAsFactors=FALSE);

##Peak annotation file
fpeakannot = "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/EpiRoadmap/jul2013.roadmapData.qc.xlsx";
##Directory containing the epigenetic annotations
fpeakdir = "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/EpiRoadmap/";
##List all the epigenetic marks to get (these are the ones where brain has pre-dominated)
activemarks = c("1_TssA","2_TssAFlnk","7_Enh","6_EnhG");
repressivemarks = c("8_ZNF/Rpts", "13_ReprPC", "14_ReprPCWk", "9_Het");
##GZ/CP differences from human fetal cortex (de la Torre-Ubieta, Stein et al., Cell)
fGZCP = "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/differentialbinding_CQN_delatorrestein_cell.csv";
##Output directory for annotation file
outputdir = "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/EpiRoadmap/";

##Read in the GZCP differences file
GZCP = read.csv(fGZCP,row.names=1);
GZCPranges = GRanges(GZCP$seqnames,IRanges(GZCP$start,GZCP$end));
mcols(GZCPranges) = GZCP[,c(4:10)];
GZCPranges = renameSeqlevels(GZCPranges,paste0("chr",seqlevels(GZCPranges)));
##Find the ranges for GZ>CP or CP>GZ
GZgreaterCP = GZCPranges[which(mcols(GZCPranges)$padj < 0.05 & mcols(GZCPranges)$log2FoldChange > 0)];
CPgreaterGZ = GZCPranges[which(mcols(GZCPranges)$padj < 0.05 & mcols(GZCPranges)$log2FoldChange < 0)];

##Read in the annotation file
peakannot = read_excel(fpeakannot,sheet=1);
peakannot = peakannot[3:nrow(peakannot),];
EID = peakannot$"Epigenome ID (EID)";

tissueepi = vector("list",length(EID));
names(tissueepi) = EID;

##Create a GRList of features with which to annotate that holds all the epi marks for each tissue
#for (i in 1:length(EID)) {
#    ##Read in the epi locations	
#    epifile = paste0(fpeakdir,"/",EID[i],"_15_coreMarks_mnemonics.bed.gz");
#    bedfile = read.table(epifile,header=FALSE,sep="\t",quote="",fill=TRUE);
#    colnames(bedfile) = c("chrom","chromStart","chromEnd","ChromState");
#    bedGR = GRanges(bedfile$chrom, IRanges(bedfile$chromStart,bedfile$chromEnd),ChromState=bedfile$ChromState);
#    activeind = match(mcols(bedGR)$ChromState,activemarks);
#    activeind = which(!is.na(activeind));
#    activeGR = bedGR[activeind];
#    repressiveind = match(mcols(bedGR)$ChromState,repressivemarks); 
#    repressiveind = which(!is.na(repressiveind));
#    repressiveGR = bedGR[repressiveind];
#    epimarksgranges = GRangesList(activeGR,repressiveGR);    
#    names(epimarksgranges) = c("active","repressive");
#    tissueepi[[i]] = epimarksgranges;
#}
#save(tissueepi,file="tissueepi.Rdata");
load('tissueepi.Rdata');

##loop over each of the tissue types
for (i in 1:length(EID)) {
    ##Concatenate all categories of interest into a Genomic Ranges List
    categories = tissueepi[[i]];
    ##Change the seqnames
    newseqlevelnames = substring(seqlevels(categories),4,nchar(seqlevels(categories)));
    categories = renameSeqlevels(categories,newseqlevelnames);
    ##Loop over chromosomes
    for (chr in 1:22) {
    	##Read in locations of SNPs from the annotation files
    	snps.table = read.table(gzfile(paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/baseline/baseline.",chr,".annot.gz")),sep="\t",header=TRUE);
    	snps = GRanges(snps.table$CHR,IRanges(snps.table$BP,snps.table$BP),SNP=snps.table$SNP,CM=snps.table$CM);
	##Loop over categories
    	indicator = matrix(0,length(snps),length(categories));
    	for (j in 1:length(categories)) {
            ##Find the snps that overlap with the category
            category = categories[[j]];
            olap.cat = findOverlaps(snps,category);        
            indicator[queryHits(olap.cat),j] = 1;
    	}    
	##Write out the new annotation file        
        snps.cat = snps;
        mcols(snps.cat) = cbind(as.data.frame(mcols(snps)),indicator);
    	outframe = data.frame(CHR=as.character(seqnames(snps.cat)),BP=start(snps.cat),SNP=mcols(snps.cat)$SNP,CM=mcols(snps.cat)$CM);
   	outframe = cbind(outframe,as.data.frame(mcols(snps.cat))[,3:(length(categories)+2)]);
    	colnames(outframe)[5:ncol(outframe)] = names(categories);
	dir.create(paste0(outputdir,EID[i]),showWarnings=FALSE);
    	gz1 = gzfile(paste0(outputdir,EID[i],"/",chr,".annot.gz"), "w")
    	write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    	close(gz1);
    }
}


categories = GRangesList(GZgreaterCP=GZgreaterCP,CPgreaterGZ=CPgreaterGZ);
##Change the seqnames
newseqlevelnames = substring(seqlevels(categories),4,nchar(seqlevels(categories)));
categories = renameSeqlevels(categories,newseqlevelnames);
##Loop over chromosomes
for (chr in 1:22) {
    ##Read in locations of SNPs from the annotation files
    snps.table = read.table(gzfile(paste0("/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/baseline/baseline.",chr,".annot.gz")),sep="\t",header=TRUE);
    snps = GRanges(snps.table$CHR,IRanges(snps.table$BP,snps.table$BP),SNP=snps.table$SNP,CM=snps.table$CM);
    ##Loop over categories
    indicator = matrix(0,length(snps),length(categories));
    for (j in 1:length(categories)) {
    	##Find the snps that overlap with the category
        category = categories[[j]];
        olap.cat = findOverlaps(snps,category);        
        indicator[queryHits(olap.cat),j] = 1;
    }    
    ##Write out the new annotation file        
    snps.cat = snps;
    mcols(snps.cat) = cbind(as.data.frame(mcols(snps)),indicator);
    outframe = data.frame(CHR=as.character(seqnames(snps.cat)),BP=start(snps.cat),SNP=mcols(snps.cat)$SNP,CM=mcols(snps.cat)$CM);
    outframe = cbind(outframe,as.data.frame(mcols(snps.cat))[,3:(length(categories)+2)]);
    colnames(outframe)[5:ncol(outframe)] = names(categories);
    dir.create(paste0(outputdir,"GZCP"),showWarnings=FALSE);
    gz1 = gzfile(paste0(outputdir,"GZCP/",chr,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    close(gz1);
}
        
    
    
