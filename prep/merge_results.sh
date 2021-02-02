#!/bin/bash

## Put this scipt to the directory where your .results files are, run it from there!

# Make the initial file with header to fill with data later
file=$(shuf -n1 -e *)
dir=${PWD}
annot=$(echo $dir | cut -d "/" -f9)
head -1 $file > merged_results_${annot}.txt

for i in *.results; do
   line=$(awk 'NR==1' $i)
   echo $line >> merged_results_${annot}.txt;
done
