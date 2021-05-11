#!/bin/bash

# another way to count SNPs in
# your annotations

# paths
bedDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds"
vcfDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/1000G_Phase3_genotypes_GRCh38"

module load bedtools/2.29.2

for bed in ${bedDir}/*.bed; do
	echo $bed
	tmp_bed=$(basename "$bed")
	for vcf in ${vcfDir}/*.vcf; do
		echo $vcf
		tmp_vcf=$(basename "$vcf")
		bed_vcf=$(vcf2bed < ${vcf})
		bedmap --echo --count --delim '\t' ${bed} ${bed_vcf} > ${bed%.bed}${tmp_vcf%.vcf}_SNPcount.bed;
	done;
done
