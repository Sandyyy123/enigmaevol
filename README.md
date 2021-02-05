# ENIGMA-evolution

All original scripts are available at Prof. Jason Stein's repository (https://bitbucket.org/jasonlouisstein/enigmaevolma6/src/master/). Original codes are adapted to our own data and research purposes.

## 1. Pre-processing

**a. Reformatting your GWAS summary statistics:**  
Reformat your GWAS summary statistics using `reformat_sumstats.sh`.  
Your summary stats must be formatted as the following:  
  
```
SNP A1 A2 FREQ1 BETA SE P N MARKER CHR BP
rs2977670 c g 0.9608 758.3807 485.5590 0.1183 10332 1:723891 1 723891
rs143225517 t c 0.8454 718.8055 232.3162 0.001974 15846 1:751756 1 751756
rs146277091 a g 0.0344 -1186.6514 501.2970 0.01793 10501 1:752478 1 752478
```
  
**b. Save your summary statistics files in Rdata format:**  
Save summary statistics files as Rdata files using `sumstats_txt2Rdata.R` and `run_sumstats_txt2Rdata.sh`.  
Rdata files are smaller in size compared to txt files. Thus, working with them will speed up your analysis in the upcoming steps.  
  
**c. List summary statistics file names:**  
Compile Rdata files names ("/path/to/dir/summary_stats.txt") to a single txt file using `list_rdata_files.R` and `run_list_rdata_files.sh`.

**d. Munge non-ancsetry regressed summary statistics:**  
Munge non-ancestry regressed `.txt` files to reformat and prepare files for the LDSC. Run `run_munge_sumstats.sh` to call `munge_sumstats.py`.

## 2. Assessing the impact of population stratification on summary statistics and Ancestry Regression

**a. Assessment Prior to ancestry regression:**  
       **a1.** Run the first correlation test between your summary statistics and 1000G phase 3 PC loadings (first 20 PCs).  
        Use `1000G_PC_cor_BJK_noGC.R` and `run_1000G_PC_cor_BJK_noGC.sh`.  
        Plot your results using `plot1000G_PC_cor_noGC_BJK.R` and `run_plot1000G_PC_cor_noGC_BJK.sh` (Fig.1a from Tilot et al., 2019).  
       **a2.** Compute LDSC intercepts prior to ancestry regression for each brain region.
        Run `ldsc.sh`.
    
**b. Ancestry Regression - correcting for population stratification:**  
"Using a standard multivariable regression implemented in R with the lm() function, we regress the GWAS summary statistics prior to ancestry regression (Beta_strat) with the ancestry PCs (Beta_PCs). The residuals of this model are then ancestry corrected effect sizes (Beta_r) and also have ancestry corrected standard errors and P-values." (from original repo.)  
Use `AncestryRegression_noGC.R` and `run_AncestryRegression_noGC.sh`.  
  
**c. Assessment After ancestry regression:**  
    **c1.** Run the second correlation test between your summary statistics and 1000G phase 3 PC loadings (first 20 PCs).  
        Use `1000G_PC_cor_ancreg_BJK_noGC.R` and `run_1000G_PC_cor_ancreg_BJK_noGC.sh`.  
        Plot your results using `plot1000G_PC_cor_ancreg_noGC_BJK.R` and `run_plot1000G_PC_cor_ancreg_noGC_BJK.sh`  (Fig.1b from Tilot et al., 2019).  
    **c2.** Convert ancestry regressed summary statistics to `.txt` format using `run_sumstats_Rdata2txt.sh` and `sumstats_Rdata2txt.R`.
    **c3.** Compute LDSC intercepts after ancestry regression for each brain region.
        Change `inDir` and `outDir` variables to the folder containing ancestry regressed in `ldsc.sh` and run.
        Extract LDSC intercept values from LDSC output with `Extract_values_LDSC_logs.R`.
        Plot LDSC intercept comparison (before/after ancestry regression) using `LDSC_ancreg_before_after.R`.

## Singleton Density Score (SDS)

Download SDS values for each SNP from: https://datadryad.org/resource/doi:10.5061/dryad.kd58f.
Run a correalation analysis between SDS values and GWAS effect sizes using `run_Blocked_Jackknife_clumpedSNPs_gr_forqsub.sh` and `Blocked_Jackknife_clumpedSNPs_gr_forqsub.R`.
Plot results with `plotSDS_ancreg1KGPphase3_noGC_MA6.R`.

## Partitioned Heritability
If you will run partitioned heritability by using annotation, make a list of your annotations with 2 columns (`annots` and `baseline`). with your annotations and save as annots_and_baselines.csv
