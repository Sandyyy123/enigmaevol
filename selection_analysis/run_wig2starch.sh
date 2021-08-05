#!/bin/bash

# This script will generate shell files to unzip wig files and 
# convert wig2starch, and will submit all to the grid.
# Wrote this script because it takes ages to convert wig2starch,
# so parallelized it over chrs.

# Gokberk Alagoz
# Created on: 08.08.2021

phylop="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/primate_phylop"
phastcons="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/primate_phastcons"

mkdir ${phylop}/scripts
mkdir ${phastcons}/scripts
mkdir ${phylop}/logs
mkdir ${phastcons}/logs

for i in `seq 1 22` X Y;
	
	do
	echo "working on chr " ${i}
	phylopSF="${phylop}/scripts/chr${i}_phylop.sh"
	phastconsSF="${phastcons}/scripts/chr${i}_phastcons.sh"
	phylopLog="${phylop}/logs/chr${i}_phylop.out"
	phastconsLog="${phastcons}/logs/chr${i}_phastcons.out"
	phylopWigFn="${phylop}/chr${i}.phyloP46way.primate.wigFix"
	phastconsWigFn="${phastcons}/chr${i}.phastCons46way.primates.wigFix"

	# make a script to convert phyloP and submit
	echo '#$ -N phyloP_wig2starch
#$ -cwd
#$ -q multi.q
#$ -S /bin/bash

gunzip -c '${phylopWigFn}' | wig2starch > '${phylopWigFn%.wigFix}'.starch' > ${phylopSF}
	
	chmod a+x ${phylopSF}
	echo "Submitting chr${i} phyloP to the Grid"
	qsub -o ${phylopLog} -j y ${phylopSF}

	# make a script to convert phastCons and submit
	echo '#$ -N phastCons_wig2starch
#$ -cwd
#$ -q multi.q
#$ -S /bin/bash

gunzip -c '${phastconsWigFn}' | wig2starch > '${phastconsWigFn%.wigFix}'.starch' > ${phastconsSF}

        chmod a+x ${phastconsSF}
        echo "Submitting chr${i} phastCons to the Grid"
        qsub -o ${phastconsLog} -j y ${phastconsSF};

done

echo "All submitted!"
