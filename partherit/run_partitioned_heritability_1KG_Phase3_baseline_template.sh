#!/bin/sh
#$ -N part_herit
#$ -cwd
#$ -q multi.q
#$ -S /bin/bash

set ftrait = $1
set fannot = $2
set foutput = $3
set fbaseline = $4

set trait = `basename $ftrait`
set annot = `basename $fannot`

echo "*** Beginning LDSC partitioned heritability..."
echo "Trait: $trait"
echo "Annot: $annot"
echo "Output: $foutput"
echo "Baseline: $fbaseline"

python /home/gokala/ldsc/ldsc.py  \
--h2 ${ftrait} \
--out ${foutput} \
--frqfile-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--ref-ld-chr /data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MAe3ukw3_ancreg/annotations/${fannot}/${annot}.,/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/baselineLD_v2.2/baselineLD.  \
--w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients

echo "Finished!"
