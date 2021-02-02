## This script converts txt files to Rdata format

inputDir="/data/clusterfs/lag/users/gokala/enigma-evol/dti/"
outputDir="/data/clusterfs/lag/users/gokala/enigma-evol/dti/"

for (i in dir(inputDir, pattern="txt", all.files=F, full.names=F)) {
  tmp_dir=gsub(" ","", paste(inputDir,i), fixed=T)
  mergedGR=read.table(tmp_dir, header=T, fill=T, stringsAsFactors=F)
  RdataDir=paste(outputDir,gsub("txt", "Rdata", i))
  save(mergedGR,file=gsub(" ","", RdataDir))
}
