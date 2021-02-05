##Use R from /usr/local/R-3.2.3/bin/R
##Run gene-ontology enrichments on variants within HGEs that also impact brain structure
##mapping to genes using PsychENCODE eQTL
options(stringsAsFactors=FALSE)
library(GenomicRanges);
library(biomaRt);

##BED file containing HGE_7PCW, as well as other annotations
fannot = "allAnnots_MAe3ukw3_shortnames.bed";
##eQTL data downloaded from PsychENCODE
feqtl = "DER-08a_hg19_eQTL.significant.txt";
##SNP information corresponding to eQTL data
fsnp = "SNP_Information_Table_with_Alleles.txt";
##ancestry regressed gwas summary stats
fGWASsumstats = "/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/1000Gphase3_PC_cor/AncestryRegressionData_noGC/GWASfiles.txt"
##Clumped ancestry regressed data
fclumpdir = "/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/E3ancreg1KGP3_noGC_CLUMPED/CLUMPED5e-8";

##Match the Rdata file locations of sumstats, text file sumstats
GWASsumstats=read.table(fGWASsumstats, header=FALSE)$V1
##Parse to get trait name
tmpname = sapply(GWASsumstats,function (x) {unlist(strsplit(x,"/",fixed=TRUE))[11]});
phenoname = substr(tmpname,1,nchar(tmpname)-13);
allfileloc = data.frame(rdatafile=GWASsumstats,clumpfile=paste0(fclumpdir,"/",phenoname,".clumped"));

####################
##Full surface area
##Start by getting all the clumped SNPs only for full surface area
fullsurfind = which(phenoname=="Mean_Full_SurfArea");
clump = read.table(allfileloc$clumpfile[fullsurfind],header=TRUE);
##Loop over all SNPs
for (i in 1:nrow(clump)) {
    ##Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("/ifshome/smedland/bin/plink1.9 --bfile /ifshome/smedland/ARCHIEVE/refs/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes --r2  --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld"));
    ##Read in the LD calculated from plink
    LD = read.table('tmpld.ld',header=TRUE);
    ##Turn into a genomic ranges object
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A);
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs);
    }
}

##Read in the BED file containing the HGE_7PCW file
annot = read.table(fannot,header=FALSE);
annot = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4);
##Restrict to only HGE_7PCW
##HGE = annot[which(annot$type=="HGE.7")];
##Use all HGEs for this analysis
HGE = annot[which(annot$type=="HGE.12F" | annot$type=="HGE.12O" | annot$type=="HGE.7" | annot$type=="HGE.8")]
HGE = reduce(HGE);
##Find overlaps
olap = findOverlaps(HGE,LDSNPs);
globalSASNPsolapHGE = unique(LDSNPs$indexSNP[subjectHits(olap)]);
globalHGE = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,globalSASNPsolapHGE)))];
##Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with an HGE is: ',length(unique(globalHGE$indexSNP)),'\n');


##Overlap with eQTL data from psychENCODE
eqtl = read.table(feqtl,header=TRUE);
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end));
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)];
olap2 = findOverlaps(eqtl.GR,globalHGE);
cat('Number of HGE overlapping loci that also have an eQTL: ',length(unique(globalHGE$indexSNP[subjectHits(olap2)])),'\n');
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)]);
##Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]});
##Convert these genes to hgnc_id
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org");
geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart);
cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n');
write.csv(geneannot,file="GlobalSAHGE7PCW.csv",row.names=FALSE,quote=FALSE);

####################
##All regional surface areas
##Start by getting all the clumped SNPs only for any regional surface area
surfareaind = grep("surfavg",phenoname);
##Loop over all regions
for (j in 1:length(surfareaind)) {
    if (j==1) {
       clump = read.table(allfileloc$clumpfile[surfareaind[j]],header=TRUE);
    } else {
       if (file.exists(allfileloc$clumpfile[surfareaind[j]])) {
       	  clump = rbind(clump,read.table(allfileloc$clumpfile[surfareaind[j]],header=TRUE));
       }
    }
}

##Loop over all SNPs
for (i in 1:nrow(clump)) {
    ##Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("/ifshome/smedland/bin/plink1.9 --bfile /ifshome/smedland/ARCHIEVE/refs/ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes --r2  --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld"));
    ##Read in the LD calculated from plink
    LD = read.table('tmpld.ld',header=TRUE);
    ##Turn into a genomic ranges object
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A);
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs);
    }
}

##Make sure that all index SNPs are in this list
LDSNPs = c(GRanges(clump$CHR,IRanges(clump$BP,clump$BP),SNP=clump$SNP,indexSNP=clump$SNP),LDSNPs);

save(LDSNPs,file="RegionalLDSNPs.Rdata");
##Read in the BED file containing the HGE_7PCW file
annot = read.table(fannot,header=FALSE);
annot = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4);
##Restrict to only HGE_7PCW
##HGE = annot[which(annot$type=="HGE.7")];
##Use all HGEs for this analysis
HGE = annot[which(annot$type=="HGE.12F" | annot$type=="HGE.12O" | annot$type=="HGE.7" | annot$type=="HGE.8")]
HGE = reduce(HGE);
##Find overlaps
olap = findOverlaps(HGE,LDSNPs);
regionalSASNPsolapHGE = unique(LDSNPs$indexSNP[subjectHits(olap)]);
regionalHGE = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,regionalSASNPsolapHGE)))];
##Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with an HGE is: ',length(unique(regionalHGE$indexSNP)),'\n');

##Overlap with eQTL data from psychENCODE
eqtl = read.table(feqtl,header=TRUE);
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end));
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)];
olap2 = findOverlaps(eqtl.GR,regionalHGE);
cat('Number of HGE overlapping loci that also have an eQTL: ',length(unique(regionalHGE$indexSNP[subjectHits(olap2)])),'\n');
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)]);
##Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]});
##Convert these genes to hgnc_id
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org");
geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart);
cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n');
write.csv(geneannot,file="RegionalSAHGE7PCW.csv",row.names=FALSE,quote=FALSE);

##GO enrichment done with ggprolifer2, but can't get it to install here so doing it locally
