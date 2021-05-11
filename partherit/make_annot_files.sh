#!/bin/bash
########################################################
# Make .annot and .l2 files for LDSC part. heritability#
########################################################

# Set variables
ldscDir="/home/gokala/ldsc/"
hapMap="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_EUR_Phase3_baseline/print_snps.txt"
oneKg="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/annotations/1000G_EUR_Phase3_plink/"
mainDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds/new_annots/"

########################################################
# Make .annot files
# make_annot.py function from LDSC is used
########################################################

module load python/2.7.15

for annot in ${mainDir}*.bed; do

	tmp_base=$(basename "$annot")
	tmp_annot="${tmp_base%.*}"
	echo $tmp_annot
	mkdir ${outDir}${tmp_annot}
	
	for chr in {1..22}; do
		python ${ldscDir}make_annot.py \
			--bed-file $annot \
			--bimfile ${oneKg}1000G.EUR.QC.${chr}.bim \
			--annot-file ${outDir}${tmp_annot}/${tmp_annot}.annot.gz
	done
done

########################################################
# Make l2.ldscore files
# ldsc.py from LDSC is used
########################################################

#awk '{if ($1!="SNP") {print $1} }' w_hm3.snplist > listHM3.txt for CHR in {1..22} do python ${ldscDir}/ldsc.py
#--l2
#--bfile 1000G.EUR.QC.$CHR
#--ld-wind-cm 1
#--print-snps $hapMap
#--annot yourannot.$CHR.annot.gz
#--out yourannot.$CHR done
