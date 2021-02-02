#!/bin/bash

cd /data/clusterfs/lag/users/gokala/enigma-evol/sumstats/munged/replication_v2

for fname in *.log; do
     mv "$fname" "$(echo "$fname" | sed -r 's/[_]{2}/_/')";
done
