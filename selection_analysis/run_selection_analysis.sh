#!/bin/sh
#$ -N plink_ld
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

inDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/selection_analysis"

while read leadSNP; do
	
	parsed_varID="$(cut -d':' -f2 <<<"${leadSNP}")"
	mkdir ${inDir}/scripts
	mkdir ${inDir}/shelloutput
	shellFile="${inDir}/scripts/${parsed_varID}.sh"
	logFile="${inDir}/shelloutput/${parsed_varID}.out"

	echo '#$ -N get_ld_buddies
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

Rscript --vanilla /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/selection_analysis/selection_analysis.R '${leadSNP}'' > $shellFile

	chmod a+x $shellFile
	echo "Created the script for the cluster - submitting ${leadSNP} to grid"
	qsub -o ${logFile} -j y ${shellFile}

done < ${inDir}/lead_snps.txt
