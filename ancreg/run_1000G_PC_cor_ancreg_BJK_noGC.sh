#!/bin/bash
# Run it on bash -> bash run_1000G_PC_cor_ancreg_BJK_noGC_par.sh

#-----1000G_PC_cor_ancreg_BJK_noGC_par.sh-----

# Run correlation analysis between 1000G phase 3 PC loadings and ancestry regressed GenLang sumstats for all phenotypes

#-----Variables-----
# $rdataDir - Directory containing GWAS summary statistics (in Rdata format)
# $outDir - Directory to write Spearman's correlation test results
# $rdataList - A txt file containing the list of "/path/to/dir/summary_statisctics.Rdata" of all phenotypes
rdataDir="/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC/"
outDir="/data/clusterfs/lag/users/gokala/beta/"
rdataList="/data/clusterfs/lag/users/gokala/beta/ancreg/ancreg_sumstats_rdata_list.txt"

#-----
cd /data/clusterfs/lag/users/gokala
while read line; do
   echo $line
   LINE=$line
   tmp_file_name=$(basename "$line")
   echo $tmp_file_name
   pheno_name="$(cut -d'_' -f1,2 <<<"$tmp_file_name")"
   echo $pheno_name
   tmp_run_file="${outDir}scripts/${pheno_name}_ancreg.sh"
   echo '#!/bin/sh
#$ -N PC_cor_ancreg_BJK_noGC
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

Rscript /data/clusterfs/lag/users/gokala/1000G_PC_cor_ancreg_BJK_noGC.R' $LINE $pheno_name $outDir > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "/data/clusterfs/lag/users/gokala/beta/scripts" $tmp_run_file
done < $rdataList
