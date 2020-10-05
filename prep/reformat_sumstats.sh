#!/bin/bash
#$ -N reformat_sumstats
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

#-----Reformat Sumstats-----

# Reformat your non-GC (genomic control) corrected GWAS summary stats at the beginning of your analysis
# according to the guideline at https://bitbucket.org/jasonlouisstein/enigmaevolma6/src/master/1000Gphase3_PC_cor/

# The new format should be as the following:
# SNP A1 A2 FREQ1 BETA SE P N MARKER CHR BP
# rs2977670 c g 0.9608 758.3807 485.5590 0.1183 10332 1:723891 1 723891
# rs143225517 t c 0.8454 718.8055 232.3162 0.001974 15846 1:751756 1 751756
# rs146277091 a g 0.0344 -1186.6514 501.2970 0.01793 10501 1:752478 1 752478

#-----Variables-----
dir="/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/sumstats"

#-----
echo "Starting to reformat sumstats at $dir"
for file in $dir/*.txt; do
    echo "basename $file is being reformatted..."
    awk '{print $1, $2, $3, $4, $8, $9, $10, $16, $19, $17, $18}' $file > ${file%.txt}_formatted.txt
    sed -Ei 's/MarkerName/SNP/;s/Allele1/A1/;s/Allele2/A2/;s/Freq1/FREQ1/;s/Effect/BETA/;s/StdErr/SE/;s/P-value/P/;s/TotalSampleSize/N/;s/chrposID/MARKER/;s/POS/BP/' ${file%.txt}_formatted.txt
    chmod 777 ${file%.txt}_formatted.txt
    echo "Done!"
done
echo "All done."
#-----
