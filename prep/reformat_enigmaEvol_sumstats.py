# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
# create the input and output path variables
file_name ='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/QC_5MAD/bgenie/bgenie_raw_output/bgenie_enigmaEvo_globalValues_averaged_chr{}.out.gz'
output_name='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/raw_sumstats/bgenie_enigmaEvo_globalValues_SA_averaged.out.gz'
#empty list to store the seperate dataframes (per chromosom)
df_list = []
# lets only import what we need as files are already massive
#loop over all the chromosomes
for i in range(0, 22):
	#please check this again - format is outline on https://jmarchini.org/bgenie-usage/
	#this is now just ONE phenotype, 0-4 are the details you always need, the other 3 are phenotype specific
	# Beta, SE, and P value; so do same thing again for second phenotype (0-4 + columns 11,12,14)
	df_list.append(pd.read_csv(file_name.format(i+1),sep='\s+', usecols =[*range(0,6),7,8,10]))
#concatenate the list into one single data frame
df_list = pd.concat(df_list)
# convert all -log10p values to normal p values again - easier to handle later on
#check if p value is really in column 10 --- it's in col. 8
df_list.iloc[:,8]=10**(-df_list.iloc[:,8])
#merge chr and pos columns and use ':' as separator to get the MARKER column
df_list['MARKER'] = df_list['chr'].astype(str) + ":" + df_list['pos'].astype(str)
df_list['N'] = 18901
# reorder columns
df_list_reorder = df_list[['rsid', 'a_0', 'a_1', 'af', 'mean_total_surface_beta', 'mean_total_surface_se', 'mean_total_surface-log10p', 'N', 'MARKER', 'chr', 'pos']]
# rename columns
df_list_reorder.columns = ['SNP','A1','A2','FREQ1','BETA','SE','P', 'N','MARKER','CHR','BP']
# save the dataframe and empty the list again
df_list_reorder.to_csv(output_name,sep=" ",index = None)