#Singleton Density Score (SDS) Analysis

We asked whether alleles under recent selective pressure (over the past ~2kya) impact the surface area of the human cortex. This analysis is based on measures of selection for SNPs acquired from this dataset: https://www.ncbi.nlm.nih.gov/pubmed/27738015. SDS values for each SNP were downloaded from: https://datadryad.org/resource/doi:10.5061/dryad.kd58f.

##Correlating SDS to ancestry regressed brain GWAS summary statistics

We correlated the SDS values downloaded as described above to the ancestry regressed summary statistics (described in a different folder) using a block jackknife spearman correlation. This is implemented in `Blocked_Jackknife_clumpedSNPs_gr_forqsub.R`.

##Plotting the SDS values on the brain surface

We FDR corrected across all 34 regions tested plus global surface area, and then plotted the significant (FDR < 0.05) correlations on the brain surface. This was implemented via `plotSDS_ancreg1KGPphase3_noGC_MA6.R` which uses the surface of the Bert brain from Freesurfer, where the surface files are found in `FreesurferRegionalObjs.Rdata`.