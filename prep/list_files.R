## The folder containing GWAS summary statistics
fileloc = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/data/ukb_eur_sumstats"
outputdir = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/data/ukb_eur_sumstats/"
filename = "SA_sumstats_list.txt"
curPattern = ".txt"

## read in gwas statistics file (compiled for all traits)
fGWASsumstats=gsub(" ","",paste(outputdir,filename))
write.table(dir(fileloc, pattern=curPattern, all.files=T,full.names=T),file=fGWASsumstats,row.names=F,col.names=F,quote=F)
