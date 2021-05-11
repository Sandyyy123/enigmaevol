#!/bin/bash

# Make .aanot and l2.ldscore files for all chromosomes and annotations.
# Run this on grid as "bash run_annot_and_l2.sh"

# Variables
inDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds"
#annotList="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/amanda_annotations/annot_list.txt"

# Run LD Score Estimation

mkdir $inDir/scripts

for annot in ${inDir}/*.bed; do
   echo $annot
   tmp_annot=$(basename "$annot")
   echo $tmp_annot
   tmp_run_file="${inDir}/scripts/${tmp_annot}.sh"
   echo '
#$ -N make_l2_ldscore
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

mkdir '$inDir'/'$tmp_annot'
bash /data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds/make_annot_l2.sh ' $annot > $tmp_run_file

   chmod a+x $tmp_run_file
   echo "Created the script for cluster -> submitting ${annot} to the Grid"
   qsub -wd "/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds/scripts" $tmp_run_file

done
