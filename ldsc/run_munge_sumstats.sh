#!/bin/sh
# Gokberk Alagoz, 2020

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

#-----Variables-----
# $input - summary statistic file
# $output - outfile_name

inDir="/data/clusterfs/lag/users/gokala/enigma-evol/data/replication/thickness/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/data/replication/thickness/munged/"
sumstatsList="/data/clusterfs/lag/users/gokala/enigma-evol/data/replication/thickness/sumstats_txt_list.txt"

#-----

echo $inDir
echo $outDir
echo "Munging sumstats"
mkdir "${inDir}scripts"
mkdir "${inDir}munged"

while read line; do
   echo $line
   LINE=$line
   tmp_base_name=$(basename "$line")
   echo $tmp_base_name
   pheno_name="$(cut -d'_' -f5 <<<"$tmp_base_name")"
   echo $pheno_name
   tmp_run_file="${inDir}scripts/${pheno_name}.sh"
   output="${outDir}${tmp_base_name%.txt}_munged.txt"
   echo $line
   echo '#!/bin/sh
#$ -N munge_sumstats
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

echo '$LINE'
echo '$output'
#module load python/3.7.2
#module load ldsc/v1.0.1 

python /home/gokala/ldsc/munge_sumstats.py --sumstats '$LINE' --out '$output' --merge-alleles /data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/resources/w_hm3.snplist' > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "${inDir}scripts" $tmp_run_file
done < $sumstatsList

echo "Finished!"

#-----
