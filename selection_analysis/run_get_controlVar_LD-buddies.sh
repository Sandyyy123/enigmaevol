#!/bin/sh
#$ -N get_LDbuds_rightHem
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

echo "Starting"
Rscript get_controlVar_LD-buddies.R
echo "Done!"
