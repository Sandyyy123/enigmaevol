#!/bin/bash

# Gokberk Alagoz
# 04.05.2021
# This script merges all chromosomes for PLINK .bed files

# directories
bed_file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.1"
merge_list="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/1000G_EUR_Phase3_plink/merge_list.txt"
outDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/1000G_EUR_Phase3_plink/allchr/"

# load modules
module load plink/1.9b6

plink --bfile ${bed_file} \
      --make-bed \
      --merge-list ${merge_list} \
      --out ${outDir}1000G_Phase3_GRCh37_allchr
