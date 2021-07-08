#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#
# This script will...
# 1- Merge clumped files.
# 2- Read control variants generated
# by SNPsnap.
# 3- Run plink on control variants to
# expand control sets.
# 4- Pick random variants in LD with
# controls.
#
# Run this script with R 4.0.3
#
# Gokberk Alagoz
# created on: 02.06.2021

#####

#library(biomaRt)
options(stringsAsFactors=FALSE)

#--------------------------------------------
# PATHS

clumpedDir="/data/clusterfs/lag/users/gokala/enigma-evol/selection_analysis/clumped_sumstats/european_hemave"
genotypeFile = "/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"
controlVars = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/snpsnap_lists/european_hemave/matched_snps_annotated_filtered.txt"
outDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/selection_analysis/snpsnap_lists/european_hemave"

#--------------------------------------------
# FUNCTIONS

#assure_dir_exists #TODO

#check_for_plink_error #TODO

#make_CSA_region_list = function(clumpedDir) {
#  # gather rsID and TOTAL info lead SNPs from
#  # each summary stats
#  clumpedSNPs = data.frame()
#  for (i in dir(clumpedDir, full.names = T, pattern = "clumped")) {
#      region = sapply(i, function (x) {unlist(strsplit(as.character(x),"_",fixed=TRUE))[9]})
#      tmp_clump = read.table(i, header = TRUE)
#      tmp_df = data.frame(SNP = tmp_clump$SNP, TOTAL = tmp_clump$TOTAL)
#      clumpedSNPs = rbind(clumpedSNPs, tmp_df)
#  }
#  write.table(clumpedSNPs,paste0(clumpedDir,"/snpsnap_list/surface_area_all.txt"),row.names = F, col.names = F, quote = F)
#}

make_CSA_region_list_only_rsIDs = function(clumpedDir) {
  # gather rsID info lead SNPs from
  # each summary stats, use for SNPsnap
  # control variant generation
  clumpedSNPs = data.frame()
  for (i in dir(clumpedDir,full.names = T,pattern = "clumped")) {
    region = sapply(i, function (x) {unlist(strsplit(as.character(x),"_",fixed=TRUE))[9]})
    tmp_clump = read.table(i,header=TRUE)
    tmp_df = data.frame(SNP=tmp_clump$SNP)
    clumpedSNPs=rbind(clumpedSNPs,tmp_df)
  }
  write.table(clumpedSNPs,paste0(outDir,"/surface_area_all.txt"),row.names = F, col.names = F, quote = F)
}

#get_num_LD_buddies = function(clumpedDir, outDir) {
#  # get the rsIDs and LD buddy number from clumped file
#  for (i in dir(clumpedDir,full.names = F, pattern = ".clumped")) {
#    file_name = sapply(i, function (x) {unlist(strsplit(as.character(x),".",fixed=TRUE))[1]})
#    system(paste0("awk '{print $3\"\t\"$6}' ", paste0(clumpedDir,"/",i), " > ", paste0(outDir,"/",i,"_LDbuddy_counts.txt"))) #TODO use column names instead of field numbers
#  }
#}

#get_rsIDs_from_pos = function() { # instead of this, download annotations from SNPsnap
#                                  # when generating the control variants.
#  # read control variants, convert variants IDs
#  # (chr:pos) to rsIDs
#  controlVarIDs_table = read.table(controlVarIDs,header = T)
#  controlVar_rsIDs_table = data.frame(matrix(NA,
#                                             nrow = nrow(controlVarIDs_table),
#                                             ncol = ncol(controlVarIDs_table)))
#  colnames(controlVar_rsIDs_table) = colnames(controlVarIDs_table)
#  controlVar_rsIDs_table[,1] = controlVarIDs_table$Input_SNP
#  for (i in 1:nrow(controlVarIDs_table)) {
#    for (j in 1:ncol(controlVarIDs_table)) {
#      tmp_chr = as.numeric(strsplit(controlVarIDs_table[i,j],":",fixed=TRUE)[[1]][1])
#      tmp_pos = as.numeric(strsplit(controlVarIDs_table[i,j],":",fixed=TRUE)[[1]][2])
#      snpmart = useMart("ENSEMBL_MART_SNP", 
#                        dataset = "hsapiens_snp", 
#                        path="/biomart/martservice", 
#                        host="https://grch37.ensembl.org")
#      tmp_rsID = getBM("refsnp_id", 
#                       filters = c("chr_name","start","end"), 
#                       values = list(tmp_chr,tmp_pos,tmp_pos), 
#                       mart = snpmart)
#      controlVar_rsIDs_table[i,j] = tmp_rsID[1,1]
#    }
#  }
#}

#run_plink_ld = function(genotypeFile, controlVars, outDir) {
#  # runs plink --ld in control snps
#  
#  controlVarsTable = read.table(controlVars, header = T)
#  controlVarsTable = controlVarsTable[controlVarsTable$input_snp == args[1],]
#  
#  for (i in 1:nrow(controlVarsTable)) {
#    system(paste0("module load plink/1.9b6 \
#                   plink --bfile ",genotypeFile," --r2 dprime --ld-snp ",controlVarsTable$rsID[i]," --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.9 --out ",outDir,"/",controlVarsTable$input_snp[i],"_",controlVarsTable$set[i],"_",controlVarsTable$rsID[i]))
#    }
#}

#--------------------------------------------
# MAIN

###
# Read in clumped SNPs from all summary stats
# and merge them.
###

make_CSA_region_list_only_rsIDs(clumpedDir)

###
# Get number of SNPs in LD with each tag SNP
###

# get_num_LD_buddies(clumpedDir,outDir)

###
# Generate control variants with SNPsnap using
# guidelines in the paper.
###

#run_plink_ld(genotypeFile,controlVars,outDir)

#--------------------------------------------