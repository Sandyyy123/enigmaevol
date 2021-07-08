#!/usr/bin/env Rscript
#
# This script was built on the GO_enrichnment.R
# script from Jason Stein and Amanda Tilot.
#
# Gokberk Alagoz
# Created on: 30.06.2021
#--------------------------------------------------#

# This script finds overlapping SNPs in clumped
# summary statistics and bed files.
# Runs gene-ontology enrichment on variants within 
# provided annotations that also impact cortical
# structure by mapping variants to genes using the
# PsychENCODE eQTL list.

#--------------------------------------------------#
# Provide the arguments below and submit it to the
# grid.

args = commandArgs(trailingOnly=TRUE)
bedfile = args[1] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/adult_hge_hg19.merged.sorted.bed"
sumstatsList = args[2] # "/data/clusterfs/lag/users/gokala/enigma-evol/data/european_lr/sumstats_txt_list.txt"
clumpedSumstatsDir = args[3] # "/data/clusterfs/lag/users/gokala/enigma-evol/eqtl/clumped_sumstats/european_lr"
outDir = args[4] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_lr/results"

#--------------------------------------------------#

options(stringsAsFactors=FALSE)
library(GenomicRanges)
library(biomaRt)

#--------------------------------------------------#
# PARSE & SET PATHS

# eQTL data downloaded from PsychENCODE
feqtl = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/old_results/eqtl/DER-08a_hg19_eQTL.significant.txt"

# Read sumstats paths
GWASsumstats = read.table(sumstatsList, header=FALSE)$V1

# Parse to get trait name
tmpname = sapply(GWASsumstats,function (x) {unlist(strsplit(x,"/",fixed=TRUE))[10]})
phenoname = paste0(sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[4]}),"_",sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[5]}),"_",sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[6]}))
clumpfileloc = dir(clumpedSumstatsDir,pattern = ".clumped",full.names = T)
annot_name = unlist(strsplit(unlist(strsplit(bedfile, "/", fixed=T))[15], ".", fixed=T))[1]

#--------------------------------------------------#
# FUNCTIONS

# Full surface area (LEFT HEM.)
# Start by getting all the clumped SNPs only for full surface area

fullsurfind = which(phenoname=="surface_le_Full")
clump = read.table(clumpfileloc[fullsurfind],header=TRUE)

# Loop over all clumped SNPs

for (i in 1:nrow(clump)) {
    
    # Find all SNPs in LD (r2>0.6) with each clumped SNP
  
    system(paste0("module load plink/1.9b6 \
                   plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld_",annot_name))
    
    # Read in the LD calculated from plink
    
    LD = read.table(paste0("tmpld_",annot_name,".ld"),header=TRUE)
    
    # Turn into a genomic ranges object
    
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A)
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs)
    }
}
  
# Read in the eQTL data
eqtl = read.table(feqtl,header=TRUE)
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)]

# Read in the BED file containing the annotation file
annot = read.table(bedfile,header=FALSE)
annotGR = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4)
        
# Find overlaps
olap = findOverlaps(annotGR,LDSNPs)
globalSA_SNPs_olap_annot = unique(LDSNPs$indexSNP[subjectHits(olap)])
globalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,globalSA_SNPs_olap_annot)))]
        
# Output the number of loci that overlap with an HGE
cat('Number of loci that overlap with ',annot_name,' is: ',length(unique(globalAnnot$indexSNP)),'\n')
        
##Overlap with eQTL data from psychENCODE
olap2 = findOverlaps(eqtl.GR,globalAnnot)
cat('Number of ',annot_name,' overlapping loci that also have an eQTL: ',length(unique(globalAnnot$indexSNP[subjectHits(olap2)])),'\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
        
# Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]})

if (length(unique(globalAnnot$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart)
  cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n')
  write.csv(geneannot,file=paste0(outDir,"/leftHem_SA_",annot_name,"european_lr.csv"),row.names=FALSE,quote=FALSE)
  
} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

#--------------------------------------------------#
# FUNCTIONS

# Full surface area (RIGHT HEM.)
# Start by getting all the clumped SNPs only for full surface area

fullsurfind = which(phenoname=="surface_re_Full")
clump = read.table(clumpfileloc[fullsurfind],header=TRUE)

# Loop over all clumped SNPs

for (i in 1:nrow(clump)) {
  
  # Find all SNPs in LD (r2>0.6) with each clumped SNP
  
  system(paste0("module load plink/1.9b6 \
                   plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out tmpld_",annot_name))
  
  # Read in the LD calculated from plink
  
  LD = read.table(paste0("tmpld_",annot_name,".ld"),header=TRUE)
  
  # Turn into a genomic ranges object
  
  if (i==1) { 
    LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A)
  } else {
    LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs)
  }
}

# Read in the eQTL data
eqtl = read.table(feqtl,header=TRUE)
eqtl.GR = GRanges(gsub("chr","",eqtl$SNP_chr),IRanges(eqtl$SNP_start,eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[,c(1:8,12:15)]

# Read in the BED file containing the annotation file
annot = read.table(bedfile,header=FALSE)
annotGR = GRanges(gsub("chr","",annot$V1),IRanges(annot$V2,annot$V3),type=annot$V4)

# Find overlaps
olap = findOverlaps(annotGR,LDSNPs)
globalSA_SNPs_olap_annot = unique(LDSNPs$indexSNP[subjectHits(olap)])
globalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,globalSA_SNPs_olap_annot)))]

# Output the number of loci that overlap with an HGE
cat('Number of loci that overlap with ',annot_name,' is: ',length(unique(globalAnnot$indexSNP)),'\n')

##Overlap with eQTL data from psychENCODE
olap2 = findOverlaps(eqtl.GR,globalAnnot)
cat('Number of ',annot_name,' overlapping loci that also have an eQTL: ',length(unique(globalAnnot$indexSNP[subjectHits(olap2)])),'\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])

# Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]})

if (length(unique(globalAnnot$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart)
  cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n')
  write.csv(geneannot,file=paste0(outDir,"/rightHem_",annot_name,"european_lr.csv"),row.names=FALSE,quote=FALSE)
  
} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

#--------------------------------------------------#

# All regional surface areas
# Start by getting all the clumped SNPs only for any regional surface area


surfareaind = grep("surface",phenoname)
surfareaind = surfareaind[-7] #TODO Find a proper way to exclude Full Surface area index here!

# Loop over all regions
for (j in 1:length(surfareaind)) {
  if (j==1) {
    clump = read.table(clumpfileloc[surfareaind[7]],header=TRUE)
  } else {
    if (file.exists(clumpfileloc[surfareaind[j]])) {
      clump = rbind(clump,read.table(clumpfileloc[surfareaind[j]],header=TRUE))
    }
  }
}

  # Loop over all SNPs
  for (i in 1:nrow(clump)) {
    
    # Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("module load plink/1.9b6 \
                    plink --bfile /data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",clump$SNP[i]," --out reg_tmpld_",annot_name))
    
    # Read in the LD calculated from plink
    LD = read.table(paste0("reg_tmpld_",annot_name,".ld"),header=TRUE)
    
    # Turn into a genomic ranges object
    if (i==1) { 
      LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A)
    } else {
      LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP=LD$SNP_B,indexSNP=LD$SNP_A),LDSNPs)
    }
  }
  
# Make sure that all index SNPs are in this list
LDSNPs = c(GRanges(clump$CHR,IRanges(clump$BP,clump$BP),SNP=clump$SNP,indexSNP=clump$SNP),LDSNPs)
  
#save(LDSNPs,file="RegionalLDSNPs.Rdata")

# Find overlaps
olap = findOverlaps(annotGR,LDSNPs)
regionalSASNPsolapHGE = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP,regionalSASNPsolapHGE)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ',annot_name,' is: ',length(unique(regionalAnnot$indexSNP)),'\n')

olap2 = findOverlaps(eqtl.GR,regionalAnnot)
cat('Number of ',annot_name,' overlapping loci that also have an eQTL: ',length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])),'\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
  
# Remove the . annotation
eqtlgenes = sapply(eqtlgenes,function (x) {unlist(strsplit(x,".",fixed=TRUE))[1]})

if (length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"),filters="ensembl_gene_id",values=eqtlgenes,mart=mart)
  cat('These eQTLs impact ',length(which(geneannot=="protein_coding")),' protein-coding eGenes\n')
  write.csv(geneannot,file=paste0(outDir,"/regionalSA_",annot_name,"_european_lr.csv"),row.names=FALSE,quote=FALSE)

} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

#--------------------------------------------------#
print("#--------------------------------------------------#",quote = F)
sessionInfo()