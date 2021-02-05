#--------------------------
#---- Clumping w/PLINK ----
#--------------------------

# 

#-----Variables-----

sumstatsList="/data/clusterfs/lag/users/gokala/enigma-evol/ancreg/replication_v1/txt_sumstats/SA_sumstats_txt_list.txt"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/eqtl/clumped_ancestry_regressed_sumstats/"
genotypeFile="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"

#-----

echo "Starting to clump..."

while read i; do

	echo $i
	base_name=$(basename "$i")
	echo "#!/bin/sh
#$ -N plink_clump
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

cd '$outDir'
module load plink/1.9b6

plink --bfile '$genotypeFile' --clump '$i' --clump-r2 0.6 --clump-kb 100000 --clump-p1 5e-8 --out '${outDir}${base_name}'" >> $outDir$base_name.sh
	chmod a+x $outDir$base_name.sh
	qsub -o $base_name.out -j y "$outDir$base_name.sh";

done < $sumstatsList

echo "Done!"

# Now type "qstat" on bash and follow the status of your jobs.
