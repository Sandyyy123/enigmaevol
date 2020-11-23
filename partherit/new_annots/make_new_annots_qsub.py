import os, sys

#----------------------------------------------------------------------
# Set up the list of annotations
#----------------------------------------------------------------------
annotFile=sys.argv[1]

mainDir="/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/partherit/"

with open('/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/partherit/annots_list.txt', 'r') as f:
    annots = [line.strip() for line in f]
    for annot in annots:
        print ("***********************************")
        print (i)
        print ("***********************************")
        shellFile = "evo_partheritscripts/make_new_"+annot+"annot.sh"
        logFile = "shelloutput/make_new_"+annot+"_annot.out"
        echoLine = "echo \""+mainDir+"make_new_annots_template.tcsh "+annot+"\" > "+shellFile
        os.system(echoLine) # making the shell file that will be submitted to the queue
        os.system("chmod a+x "+shellFile) # everyone can execute it: rwxr-xr-x
        qsubLine = "qsub -o `pwd`/"+logFile+" -j y `pwd`/"+shellFile
        # os.system("qsub -o `pwd`/"+logFile+" -j y `pwd`/"+shellFile)
        os.system(qsubLine)
    f.close()
