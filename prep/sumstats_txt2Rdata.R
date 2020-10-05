## This script converts txt files to Rdata format.
## PC_cor_BJK_noGC.R script reads the output of this script.

inputDir="/data/clusterfs/lag/users/gokala/enigma-evol/"
outputDir="/data/clusterfs/lag/users/gokala/enigma-evol/"
#P:/workspaces/lg-genlang/Working/Evolution/sumstatsRdata/

for (i in dir(inputDir, pattern="txt", all.files=F, full.names=F)) {
  tmp_dir=gsub(" ","", paste(inputDir,i), fixed=T)
  tmp_ss_table=read.table(tmp_dir, header=T, fill=T, stringsAsFactors=F)
  RdataDir=paste(outputDir,gsub("txt", "Rdata", i))
  save(tmp_ss_table,file=gsub(" ","", RdataDir))
}
