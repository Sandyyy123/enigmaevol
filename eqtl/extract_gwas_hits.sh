#!/bin/bash
# Gokberk Alagoz - 05/10/21
# This script will extract info of genome-wide significant SNPs
# from the relevant summary stats file.

var_list='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_lr/results/results/regLeft_list_clean.txt'
outDir='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_lr/results/results'

#set -xv

while read line; do
      	echo $line 
	rsID=$(cut -d'	' -f2 <<< $line)
	region=$(cut -d'	' -f1 <<< $line)

#	sumstat="/data/clusterfs/lag/users/gokala/enigma-evol/data/european_lr/sumstats_ukb43760_regionalDK_surface_le_${region}_withGlob_european_allChr.txt"
	sumstat="sumstats_ukb43760_regionalDK_surface_le_${region}_withGlob_european_allChr.txt"
	
	grep -Fw $rsID $sumstat >> test_file.txt

done < $var_list
