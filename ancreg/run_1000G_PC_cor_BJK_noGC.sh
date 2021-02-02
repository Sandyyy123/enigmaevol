#!/bin/bash
#$ -N run_PC_cor_BJK_noGC
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----run_1000G_PC_cor_BJK_noGC.sh-----

# Run correlation analysis between 1000G phase 3 PC loadings and sumstats BETAs for all phenotypes
# /usr/local/lib64/R/bin/
#-----Variables-----
# $rdataDir - Directory containing GWAS summary statistics (in Rdata format)
# $outDir - Directory to write Spearman's correlation test results
# $rdataList - A txt file containing the list of Rdata files with their paths
rdataDir="/data/clusterfs/lag/users/gokala/enigma-evol/dti"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/dti/corvals/"
rdataList="/data/clusterfs/lag/users/gokala/enigma-evol/dti/sumstats_rdata_list.txt"

#-----
# Read rdataList
mkdir ${outDir}scripts
cd $rdataDir

while read line; do 
   echo $line
   LINE=$line
   tmp_file_name=$(basename "$line")
   echo $tmp_file_name
   pheno_name="$(cut -d'_' -f3,6,7 <<<"$tmp_file_name")"
   echo $pheno_name
   tmp_run_file="${outDir}scripts/${pheno_name}.sh"
   echo '#!/bin/sh
#$ -N PC_cor_BJK_noGC
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

Rscript /data/clusterfs/lag/users/gokala/enigma-evol/1000G_PC_cor_BJK_noGC.R' $LINE $pheno_name $outDir > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd $outDir"scripts" $tmp_run_file
done < $rdataList
