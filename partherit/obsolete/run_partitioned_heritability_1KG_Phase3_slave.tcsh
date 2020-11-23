#!/bin/tcsh

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

/ifshome/smedland/bin/anaconda2/bin/python /ifshome/smedland/bin/ldsc/ldsc.py  \
--h2 /ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/sumstats/w_1KGP3_ancreg/${ftrait} \
--out ${foutput} \
--frqfile-chr /ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/resources/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--ref-ld-chr /ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/partherit/annots_1KGPhase3/${fannot}/${annot}.,/ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/partherit/annots_1KGPhase3/${fbaseline}/${fbaseline}. \
--w-ld-chr /ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients

echo "Finished!"
