#!/bin/sh
#$ -N list_files
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

echo "Starting"
/usr/local/bin/Rscript list_files.R
echo "Done!"
