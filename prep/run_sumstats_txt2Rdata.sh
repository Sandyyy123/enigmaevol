#!/bin/bash
#$ -N run_sumstats_txt2Rdata
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

#-----Convert sumstats.txt to sumstats.Rdata-----

# Convert txt formatted summary statisctics files to Rdata for downstream analysis.
# Rdata files will be used by PC_corBJK_noGC.R script.

#-----
echo "Starting to convert sumstats from txt to Rdata"
/usr/local/bin/Rscript sumstats_txt2Rdata.R
echo "Done!"
#-----
