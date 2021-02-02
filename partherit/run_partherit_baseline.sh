#-----Run Partitioned Heritability-----

# Estimating heritability, genetic correlation and the LD Score Regression Intercept
# ldsc function is from github.com/bulik/ldsc
# For interpreting results, check the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

#-----Variables-----
# $inDir - ancestry regressed + munged summary statistics directory
# $iv - LD Score files to use as the independent variable in the LD Score regression
# $rw - LD Scores to use for the regression weights
# For this analysis, same file is used as iv and rw -> eur_w_ld_chr

inDir="/data/clusterfs/lag/users/gokala/enigma-evol/sumstats/munged/replication_v1/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/nonancreg_intercepts/replication_v1/"
iv_rw="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/"
mungedSumstatsList="/data/clusterfs/lag/users/gokala/enigma-evol/sumstats/munged/replication_v1/munged_sumstat_list.txt"

#-----
mkdir "${inDir}scripts"

echo $inDir
echo $outDir
echo $iv_rw
echo "Starting to compute partitioned heritability"

#module load python/3.7.2 \

while read line; do
   echo $line
   tmp_base_name=$(basename "$line")
   echo $tmp_base_name
   trait="$(cut -d'_' -f1,2 <<<"$tmp_base_name")"
   echo $pheno_name
   output="${outDir}${trait}"
   tmp_run_file="${inDir}scripts/${trait}.sh"
   echo '#!/bin/sh
   #$ -N partherit
   #$ -cwd
   #$ -q multi15.q
   #$ -S /bin/bash

   /data/clusterfs/lag/users/gokala/enigma-evol/partherit/run_partitioned_heritability_1KG_Phase3_baseline_template.sh '$tmp_base_name' '$annot' '$outDir'' > tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "${inDir}scripts" $tmp_run_file
done < $mungedSumstatsList

echo "Done!"

#-----

#!/bin/sh \n"+"#$ -N partherit \n"+"#$ -cwd \n"+"#$ -q multi15.q \n"+"#$ -S /bin/bash \n"+mainDir+"run_partitioned_heritability_1KG_Phase3_baseline_plus_extra_template.sh "+E3MA+" "+annot+" "+outDir+baseE3MA+" "+baseline+"\" > "+shellFile

