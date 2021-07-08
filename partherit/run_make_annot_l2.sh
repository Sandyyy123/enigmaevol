#!/bin/bash

# Make .annot and l2.ldscore files for all chromosomes and annotations.
# Run this on grid as "bash run_make_annot_l2.sh"

# Variables
inDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37"

# Run LD Score Estimation

mkdir $inDir/scripts

for annot in ${inDir}/nean_RA_hg19-rinker_et_al.sorted.bed; do
   echo $annot
   tmp_annot=$(basename "$annot")
   echo $tmp_annot
   tmp_run_file="${inDir}/scripts/${tmp_annot}.sh"
   echo '#$ -N make_annot_and_ldscore
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

export PATH="/home/gokala/bedtools2/bin:$PATH"
mkdir '$inDir'/'$tmp_annot'
bash /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/partherit/make_annot_l2.sh ' $annot > $tmp_run_file

   chmod a+x $tmp_run_file
   echo "Created the script for cluster -> submitting ${annot} to the Grid"
   qsub -wd "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotation/beds/GRCh37/scripts" $tmp_run_file

done
