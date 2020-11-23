#!/bin/tcsh

set annot = /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/EpiRoadmap_FetalscRNASeq
set baseannot = "FetalscRNASeq"
foreach chr (`seq 1 22`) 
    echo "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/ldscore_slave.tcsh $chr ${annot}" > partheritscripts/${baseannot}_${chr}.sh
    chmod a+x partheritscripts/${baseannot}_${chr}.sh
    qsub -o `pwd`/shelloutput/${baseannot}_${chr}.out -j y `pwd`/partheritscripts/${baseannot}_${chr}.sh
end




##Loop over each annotation for a tissue
##foreach baseannot (`cat /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/consolidated.txt`) 
##foreach baseannot (`cat /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/fetalbrainadultspecific.txt`) 
##    set annot = /ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/EpiRoadmap_FetalAdultSpecific/${baseannot}/
##    foreach chr (`seq 1 22`) 
##	echo "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3UKBB/ldscorescripts/partherit/ldscore_slave.tcsh $chr ${annot}" > partheritscripts/${baseannot}_${chr}.sh
##	chmod a+x partheritscripts/${baseannot}_${chr}.sh
##	qsub -o `pwd`/shelloutput/${baseannot}_${chr}.out -j y `pwd`/partheritscripts/${baseannot}_${chr}.sh
##    end
##end

