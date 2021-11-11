# Gokberk Alagoz - 30.09.2021
# This script will run PLINK on 5k control variants of each
# SA-associated lead SNP and get LD-buddies for them.

#-----Variables-----
# $inDir - ancestry regressed + munged summary statistics directory

genotypeFile = "/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"
controlVars = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication/SNPsnap_preterm_birth/matched_snps_annotated_subset.filtered.txt"
outDir= "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/preterm_birth_replication"

#-----

library(data.table)
options(stringsAsFactors=FALSE)

# runs plink --ld on control snps

controlVarsTable = read.table(controlVars, header = T)

for (i in 1:nrow(controlVarsTable)) {

      	system(paste0("module load plink/1.9b6 \
      	               plink --bfile ", genotypeFile, " --r2 dprime --ld-snp ",
                       controlVarsTable$rsID[i], " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.9 --out ",
                       outDir, "/", 
		       controlVarsTable$input_snp[i], "_", 
		       controlVarsTable$set[i], "_",
                       controlVarsTable$rsID[i]))
}
