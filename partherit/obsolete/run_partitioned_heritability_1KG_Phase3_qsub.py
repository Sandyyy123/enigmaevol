# Updated file to work with 1000G Phase 3 data!
# Notes: 
# This is designed to take a single argument that should be a csv file with two columns:
# "Annotation" and "Baseline"
# e.g. "/ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/partherit/annots_and_baselines.csv"

# The way to run this on the USC servers is:
#  /ifshome/smedland/bin/anaconda2/bin/python run_partitioned_heritability_1KG_Phase3_qsub.py\
#  annots_and_baselines.csv \
#  /ifs/loni/faculty/dhibar/ENIGMA3/MAe3ukw3/evolution/LDSC/1000G_phase3_ancreg_Rdata_sumstats_noGC.txt

# This way all of the information on what we ran and how is contained in a single file, and there aren't
# 800 separate qsub files to keep track of.

import os, sys
import pandas as pd

#----------------------------------------------------------------------
# Set up the list of annotations and other folders
#----------------------------------------------------------------------
annotFile=sys.argv[1]
ancregSumstats=sys.argv[2]

mainDir="/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/partherit/"

annots = pd.read_csv(annotFile, header = 0)
print (annots.head(n=5))

#----------------------------------------------------------------------
# Loop over the annotations and ancestry regressed GWAS sumstats
#----------------------------------------------------------------------
for index, row in annots.iterrows():
	annot = row['Annotation']
	baseline = row['Baseline']
	print ("***********************************")
	print (annot)
	print (baseline)
	print ("***********************************")
	with open(ancregSumstats) as sumstats:
		sumstat_list = [line.rstrip('\n') for line in sumstats]
		for E3MA in sumstat_list:
			baseE3MA = os.path.basename(E3MA)
			outDir = mainDir+"ancreg_Phase3_results/"+annot+"/"
			shellFile = "evo_partheritscripts/"+annot+"_vs_"+baseline+"_"+baseE3MA+".sh"
			logFile = "shelloutput/"+annot+"_vs_"+baseline+"_"+baseE3MA+".out"
			if not os.path.isdir(outDir):
				os.mkdir(outDir)

			# First check if there are already results for this annotation
			# and if not, continue running the analysis.
			if not os.path.exists(outDir+baseE3MA+".results"):
				if not os.path.exists(shellFile):
					echoLine = "echo \""+mainDir+"run_partitioned_heritability_1KG_Phase3_slave.tcsh "+baseE3MA+" "+annot+" "+outDir+baseE3MA+" "+baseline+"\" > "+shellFile
				else: 
					print ("There are already .sh files for this annotation, please go delete them first. :) ")
					exit() # we'll exit the script here so that the user can go delete the shell files
				# os.system("echo \""+mainDir+"scripts/run_partitioned_heritability_1KG_Phase3_slave.tcsh "+baseE3MA+" "+annot+" "+outDir+baseE3MA+" "+baseline+"\" > "+shellFile)
				os.system(echoLine) # making the shell file that will be submitted to the queue
				os.system("chmod a+x "+shellFile) # everyone can execute it: rwxr-xr-x
				qsubLine = "qsub -o `pwd`/"+logFile+" -j y `pwd`/"+shellFile
				# os.system("qsub -o `pwd`/"+logFile+" -j y `pwd`/"+shellFile)
				os.system(qsubLine)
			# If there *are* results files for that annotation, we want to 
			# move on to the next sumstats file, and then on to the next annotation 
			# if every set of sumstats already has results.
			else:
				continue
	sumstats.close()
