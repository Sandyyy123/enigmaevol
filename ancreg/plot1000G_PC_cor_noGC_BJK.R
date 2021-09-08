##/usr/local/R-3.2.3/bin/R
##Run Correlation of effect sizes from GWAS of 1000G PC components with effect sizes from GIANT height after ancestry regression
##The effects sizes from GWAS of 1000G PC components were acquired from Katya and are from phase 3
##They were calculated by deriving PCs from 1000G (all populations) and correlating that with SNPs
##The goal here is to see if population stratification is driving the results
 
##Install hexbin package for making hexagon scatterplots
##install.packages("hexbin")
##library(hexbin)
options(stringsAsFactors=FALSE)
#library(GenomicRanges);

dircorvals = "/data/clusterfs/lag/users/gokala/enigma-evol/corvals/replication/surface/"
#"P:/workspaces/lg-genlang/Working/Evolution/all_corvals/"
##Output file
foutput = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/plots/ancestry_correlations/raw/anc_cor_SA_wGlobCov_onePerPage.pdf"

fcorvals = dir(dircorvals,full.names=TRUE,pattern="csv");

##Make a pdf file of correlation estimates
pdf(foutput);
par(las=2);

for (i in 1:length(fcorvals)) {

    ##Read in correlation values
    corvals = read.csv(fcorvals[i]);

    ##Make barplots of correlation estimates for each
    corind = grep("BJK_cor",colnames(corvals));
    pind = grep("BJK_P",colnames(corvals));
    seind = grep("BJK_SE",colnames(corvals));
    region = strsplit(corvals$X[1], split = "_")[[1]][2]
    
    x = barplot(as.matrix(corvals[1,corind]),main=region,ylab="correlation coefficient",xlab="Ancestry PC",names.arg=paste0("PC",seq(1,20)),ylim=c(-0.05,0.05));
    y0 = as.numeric(corvals[1,corind]-corvals[1,seind]);
    y1 = as.numeric(corvals[1,corind]+corvals[1,seind]);
    arrows(x,y0,x,y1,angle=90,length=0)
    ##Add asterisk when significant (Bonferroni corrected)
    bonfsigind = which(corvals[1,pind] < 0.05/20);
    ## draw asterisk for significant ones
    ##text(x[sigind]+((x[sigind+1]-x[sigind])/2),0.095,"*");
    text(x[bonfsigind],0.095,"*");
    ##Add o when nominally significant
    nomsigind = which(corvals[1,pind] >= 0.05/20 & corvals[1,pind] < 0.05);
    text(x[nomsigind],0.095,"o");

}
dev.off();

##Output file
foutput = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/plots/ancestry_correlations/raw/anc_cor_SA_wGlobCov.pdf"

fcorvals = dir(dircorvals,full.names=TRUE,pattern="csv");
ind = grep("corvalues",fcorvals);
fcorvals = fcorvals[ind];

##Make a pdf file of correlation estimates
pdf(foutput,width=8.5,height=11);
par(las=2,mfrow=c(3,2));

for (i in 1:length(fcorvals)) {

    ##Read in correlation values
    corvals = read.csv(fcorvals[i])

    ##Make barplots of correlation estimates for each
    corind = grep("BJK_cor",colnames(corvals));
    pind = grep("BJK_P",colnames(corvals));
    seind = grep("BJK_SE",colnames(corvals));
    region = strsplit(corvals$X[1], split = "_")[[1]][2]
    
    x = barplot(as.matrix(corvals[1,corind]),main=region,ylab="correlation coefficient",xlab="Ancestry PC",names.arg=paste0("PC",seq(1,20)),ylim=c(-0.05,0.05));
    y0 = as.numeric(corvals[1,corind]-corvals[1,seind]);
    y1 = as.numeric(corvals[1,corind]+corvals[1,seind]);
    arrows(x,y0,x,y1,angle=90,length=0)
    ##Add asterisk when significant (Bonferroni corrected)
    bonfsigind = which(corvals[1,pind] < 0.05/20);
    ## draw asterisk for significant ones
    ##text(x[sigind]+((x[sigind+1]-x[sigind])/2),0.095,"*");
    text(x[bonfsigind],0.095,"*");
    ##Add o when nominally significant
    nomsigind = which(corvals[1,pind] >= 0.05/20 & corvals[1,pind] < 0.05);
    text(x[nomsigind],0.095,"o");

}
dev.off()
