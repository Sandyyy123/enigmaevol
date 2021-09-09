#!/bin/bash
# Gokberk Alagoz
# 01.06.2021
#
# This script derives CSA genomic regions from GWAS summary statistics
# and creates matched control regions for CSA-associated regions.
# All pipeline is based on LaBella et al. (2020).
#---------------------------------------------------------------------- 
sumstatsDir="/data/clusterfs/lag/users/gokala/enigma-evol/data/european_lr"
sortedDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/data/european_lr/pvalue_sorted"
genotypeFile="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/selection_analysis/clumped_sumstats/european_lr/clumped_sumstats_p10e-5"

######################################################
# Sort summary statistics based on pvalue (descending)
######################################################

#echo "Sorting sumstats based on pvalues"
#for i in ${sumstatsDir}/*allChr.txt; do
#	tmp_name=$(basename "$i")
#	sort -k7 -n ${i} > ${sortedDir}/${tmp_name%.txt}_Psorted.txt
#done
#echo "Done!"

#####################
# Take first 10k SNPs
#####################

#echo "Subsetting sumstats - getting top 10k"
#for i in ${sortedDir}/*.txt; do
#	tmp_name=$(basename "$i")
#	head -n 10001 ${i} > ${sortedDir}/${tmp_name%.txt}_top10k.txt
#done
#echo "Done!"

#########################################################################
# Clump 10k SNPs baed on LD: p <= e10-4, r^2 > 0.9
# "A liberal p value threshold is used to increase the # of CSA-associated
# variants evaluated. We anticipate false positives will not have sig.
# evolutionary signals." (LaBella et al. 2020)
#########################################################################

mkdir ${outDir}/scripts
mkdir ${outDir}/shell_logs
echo "Starting to clump..."
for i in ${sortedDir}/*_top10k.txt; do
	base_name=$(basename "${i}")
	echo "Working on ${base_name}"
	echo "#!/bin/sh
#$ -N plink_clump
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
cd '${outDir}'
module load plink/1.9b6
plink --bfile '${genotypeFile}' --clump '${i}' --clump-r2 0.9 --clump-p1 10e-5 --out '${outDir}/${base_name%.txt}'" >> ${outDir}/scripts/${base_name%.txt}.sh
	chmod a+x ${outDir}/scripts/${base_name%.txt}.sh
	qsub -o ${outDir}/shell_logs/${base_name%.txt}.out -j y "${outDir}/scripts/${base_name%.txt}.sh"
done

echo "Done!"

# Continue with creating matched control regions for CSA-associated regions
# with SNPsnap website.
