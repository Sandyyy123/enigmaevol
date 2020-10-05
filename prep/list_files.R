## The folder containing GWAS summary statistics
fileloc = "/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC"
outputdir = "/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC/"
filename = "sumstats_txt_list.txt"
curPattern = "txt"
#P:/workspaces/lg-genlang/Working/Evolution/

## read in gwas statistics file (compiled for all traits)
fGWASsumstats=gsub(" ","",paste(outputdir,filename))
write.table(dir(fileloc, pattern=curPattern, all.files=T,full.names=T),file=fGWASsumstats,row.names=F,col.names=F,quote=F)
