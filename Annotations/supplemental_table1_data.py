# import pandas as pd 
import numpy as np 
import os, sys, glob

bedDir = "/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/supplemental_table1/beds/"

plinkDir = "/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/"


os.chdir(bedDir)
annots = glob.glob("*.bed")


for i in annots:
	annot = i[:-4]
	resultDir = "/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/supplemental_table1/"
	if not os.path.exists(resultDir):
		os.makedirs(resultDir)
	print("***********Running PLINK to caclulate the MAF for "+annot+" SNPs")
	os.system("plink \
		--bfile "+plinkDir+"/1KG_phase3_GRCh37_allchr \
		--keep "+resultDir+"unrelated_EUR_1KG_phase1_PLINKformat_v2.txt \
		--freq \
		--extract range "+i+" \
		--out "+resultDir+"1000G.Phase3."+annot)
print("*******All done with  "+annot+"!")
