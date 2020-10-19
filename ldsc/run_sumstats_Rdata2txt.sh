#!/bin/bash
#$ -N run_sumstats_Rdata2txt
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Convert sumstats.Rdata to sumstats.txt-----

# Convert Rdata formatted summary statistics files to txt for downstream analysis.

#-----
echo "Starting to convert sumstats from Rdata to txt"
Rscript sumstats_Rdata2txt.R
echo "Done!"
#-----
