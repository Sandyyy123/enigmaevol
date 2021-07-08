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
rdataDir="/data/clusterfs/lag/users/gokala/enigma-evol/Rdata/replication"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/corvals/replication/"
rdataList="/data/clusterfs/lag/users/gokala/enigma-evol/Rdata/replication/munged_sumstats_list.txt"

#-----
# Read rdataList
mkdir ${outDir}scripts
cd $rdataDir

while read line; do 
   echo $line
   LINE=$line
   tmp_file_name=$(basename "$line")
   echo $tmp_file_name
   pheno_name="$(cut -d'_' -f4,5,6 <<<"$tmp_file_name")"
   echo $pheno_name
   tmp_run_file="${outDir}scripts/${pheno_name}.sh"
   echo '#!/bin/sh
#$ -N PC_cor_BJK_noGC
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

Rscript /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/ancreg/1000G_PC_cor_BJK_noGC.R' $LINE $pheno_name $outDir > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd $outDir"scripts" $tmp_run_file
done < $rdataList
