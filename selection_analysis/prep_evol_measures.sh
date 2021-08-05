#!/bin/bash

# This script explains how to get phyloP, phastCons and GERP scores from UCSC servers and how to match them
# with your variants.

# Gokberk Alagoz
# Created on: 31/07/2021

# Download phyloP46way and phastCons46way data sets from UCSC server. They are in bigWig format.
# Download the bigWig toolkit from UCSC: rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./
# Add those functions to your bash path so you can call them directly from the command line. Actually, do it for any
# tool that you install yourself :).

module load bedtools/2.29.2

# Convert bigWig files to wig
bigWigToWig primates.phyloP46way.bw primates.phyloP46way.wig
bigWigToWig primates.phastCons46way.bw primates.phastCons46way.wig

# Convert wig files to bed using wig2bed from BEDOPS
wig2bed < primates.phyloP46way.wig > primates.phyloP46way.bed
wig2bed < primates.phastCons46way.wig > primates.phastCons46way.bed

###########################
# Full bigWig files which should contain all chromosomes seem to have only
# chr9, chrX and chrY. So I'll try independent chr files version of phastCons
# and phyloP scores.

# Download data
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/primates /path/to/dir
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons46way/primates /path/to/dir
