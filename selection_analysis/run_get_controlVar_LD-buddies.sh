#!/bin/bash
#
# Gokberk Alagoz - 22.10.21
#
# This script will parallelize over each control SNP
# and will run PLINK --ld-snp.
# Gridmaster does not work if you have a colon in output
# file name. Replace colons with underscores.
# Run this via gridportal1, as grid cannot find PLINK
# for some reason.

# SET PATHS
input_snps="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication/SNPsnap_preterm_birth/input_snps_identifer_mapping.txt"
genotypeF="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"
controlVars="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication/SNPsnap_preterm_birth/matched_snps_annotated_subset.filtered.txt"
outDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication"

# MAIN

echo "******************************************"
echo "Creating output files for each clumped SNP"
echo "******************************************"

while read snpID rsID; do

        snpID=$(echo "${snpID}" | sed -r 's/[:]+/_/g')
	mkdir "${outDir}/${snpID}"
	mkdir "${outDir}/${snpID}/scripts"
	mkdir "${outDir}/${snpID}/results"
	mkdir "${outDir}/${snpID}/shelloutput";

done < ${input_snps}

echo "*************************************"
echo "Starting to submit PLINK jobs to grid"
echo "*************************************"

while read set_col input_snp snpID rsID junk; do

	echo "Trait-associated SNP: $input_snp"
	echo "Control SNP: $rsID"

	input_snp=$(echo "${input_snp}" | sed -r 's/[:]+/_/g')
	snpID=$(echo "${snpID}" | sed -r 's/[:]+/_/g')

	shellFile="${outDir}/${input_snp}/scripts/chr${input_snp}_${set_col}_${rsID}.sh"
        logFile="${outDir}/${input_snp}/shelloutput/${input_snp}_${set_col}_${rsID}.out"

	echo '#!/bin/sh
#$ -N get_LD_buds
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load plink/1.9b6
	
plink --bfile '${genotypeF}' --r2 dprime --ld-snp '${rsID}' \
	--ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.9 \
	--out '${outDir}'/'${input_snp}'/'${input_snp}'_'${set_col}'_'${rsID}'' > ${shellFile}

        chmod a+x ${shellFile}
        echo "Created the script for cluster. Submitting job to the grid."
        qsub -o ${logFile} -j y ${shellFile};

done < ${controlVars}
