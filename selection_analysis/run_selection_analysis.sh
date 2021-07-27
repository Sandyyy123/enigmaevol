#!/bin/sh
#
# This script will make sub-scripts for each lead SNP
# and will run selection_analysis.R script in grid.
# selection_analysis.R script will run plink and
# get LD-buddies of control variants.
#

inDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/european_hemave"

mkdir ${inDir}/scripts
mkdir ${inDir}/shelloutput

while read leadSNP; do
	
	parsed_varID="$(cut -d':' -f2 <<<"${leadSNP}")"
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

done < ${inDir}/SNPsnap_european_hemave/working_snps.txt #working_snps list is the list of lead SNPs survived SNPsnap filtering.
