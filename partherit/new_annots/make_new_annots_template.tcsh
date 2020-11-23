#!/bin/tcsh

set annot = $1


echo "*** Beginning LDSC partitioned heritability..."
echo "Annot: $annot"

/ifshome/smedland/bin/anaconda2/bin/python \
/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/partherit/Step2_make_evo_based_annots_1KGPhase3.py ${annot} 


echo "Finished!"