#!/bin/bash
#
# Gokberk Alagoz
# Created on: 08.06.2021
#
#-----------------------
# This script will...
# 1) Extract set, input_snp, snpID and rsID
# columns from the matched_snps_annotated.txt
# 2) Remove control variants wo an rsID

#-----------------------
# PATHS

matched_snps_annotatedDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication/SNPsnap_preterm_birth"

#-----------------------
# MAIN

awk '{print $1"\t"$2"\t"$3"\t"$4}' $matched_snps_annotatedDir/matched_snps_annotated.txt > $matched_snps_annotatedDir/matched_snps_annotated_subset.txt

awk '$4 ~ /^rs/ {print}' $matched_snps_annotatedDir/matched_snps_annotated_subset.txt > $matched_snps_annotatedDir/matched_snps_annotated_subset.filtered.txt
