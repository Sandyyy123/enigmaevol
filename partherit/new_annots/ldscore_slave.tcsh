#!/bin/tcsh

set chr = $1
set fannot = $2

/ifshome/smedland/bin/anaconda2/bin/python /ifshome/smedland/bin/ldsc/ldsc.py --l2 --bfile /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/1000G_plinkfiles/1000G.mac5eur.$chr --ld-wind-cm 1 --annot $fannot/$chr.annot.gz --out $fannot/$chr --print-snps /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/hapmap3_snps/hm.$chr.snp



