#fGWASsumstats = "P:/workspaces/lg-genlang/Working/Evolution/sumstats_rdata_list.txt"
#fLDSCint = "P:/workspaces/lg-genlang/Working/Evolution/BETA.ancBETA_WR_Eur.csv";

##/usr/local/R-3.2.3/bin/R
##Correlate the effect sizes before and after ancestry regression
options(stringsAsFactors=FALSE)
library(GenomicRanges);

##Rdata files containing GWAS summary statistics
rdatafileloc = "P:/workspaces/lg-genlang/Working/Evolution/AncestryRegressionData_noGC/"
#"/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC/ancreg_Rdata"
##read in gwas statistics file (compiled for all traits)
fGWASsumstats = "P:/workspaces/lg-genlang/Working/Evolution/AncestryRegressionData_noGC/sumstats_rdata_list.txt"
#"/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC/ancreg_Rdata/sumstats_rdata_list.txt"
##LDSC intercept values from original ENIGMA (not ancestry regressed)
fLDSCint = "P:/workspaces/lg-genlang/Working/Evolution/LDSC_intercepts_w_and_wo_ancreg.csv"
#"/data/clusterfs/lag/users/gokala/LDSC/LDSC_intercepts_w_and_wo_ancreg.csv"

##Block jackknife function
cor.test.jack = function(x, y, blocks=1000) {
  
  keep = is.finite(x) & is.finite(y)
  x = x[keep]
  y = y[keep]
  mat = cbind(x, y, floor(seq(1, blocks+1-1e-3, len=length(x))))
  part_ests = sapply(1:blocks, function(i) cor(mat[mat[,3] != i, 1], mat[mat[,3] != i, 2]))
  jack_est = mean(part_ests)
  jack_se = sqrt((blocks-1)/blocks * sum((jack_est - part_ests)^2))
  jack_z = jack_est/jack_se;
  jack_p = 2*pnorm( -abs(jack_z));
  list(estimate=jack_est, se=jack_se, p=jack_p, estimate_normal=cor(x, y), se_normal=cor.test.plus(x, y)$se)
  
}

cor.test.plus <- function(x, y, ...) {
  
  # like cor.test, but also returns se of correlation
  corr = cor.test(x, y, ...)
  corr$se = unname(sqrt((1 - corr$estimate^2)/corr$parameter))
  corr
  
}

##Match the Rdata file locations of sumstats, text file sumstats, and clumped files
GWASsumstats=read.table(fGWASsumstats, header=FALSE)$V1;
##Parse to get trait name
tmpname = sapply(GWASsumstats,function (x) {unlist(strsplit(x,"/",fixed=TRUE))[7]});
tmpnamepheno = paste(sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[1]}),sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]}),sep="_")
phenoname = tmpnamepheno

##Make output matrices for the correlations
output = matrix(NA,nrow=length(phenoname),ncol=3);
rownames(output) = phenoname;
colnames(output) = c("corBETA","corSE","corP");

##Perform correlations pre and post ancestry regression
for (i in 1:length(GWASsumstats)) {
  load(GWASsumstats[i]);
  print(GWASsumstats[i])
  coroutput=cor.test.jack(mergedGR$BETA,mergedGR$ancBETA);
  output[i,1] = coroutput$estimate;
  output[i,2] = coroutput$se;
  output[i,3] = coroutput$p;
}

write.csv(output,file="P:/workspaces/lg-genlang/Working/Evolution/BETA.ancBETA_noGC_BJK.csv",quote=FALSE);
#"/data/clusterfs/lag/users/gokala/LDSC/BETA.ancBETA_noGC_BJK.csv"
##output = read.csv('P:/workspaces/lg-genlang/Working/Evolution/BETA.ancBETA_noGC_BJK.csv',row.names=1);
stop();

##Plot the results
pdf("P:/workspaces/lg-genlang/Working/Evolution/BETA.ancBETA_noGC_BJK.pdf");
#"/data/clusterfs/lag/users/gokala/LDSC/BETA.ancBETA_noGC_BJK.pdf"
plot(1:11,output[,1],ylab="cor(BETA,ancestry adj BETA)",xlab="GWAS",pch=19,xaxt="n",main="Effect sizes");
axis(1,at=1:11,labels=phenoname,las=2,cex.axis=0.5);

#pdf("/data/clusterfs/lag/users/gokala/LDSC/BETA.ancBETA_noGC_BJK_v2.pdf");
#plot(1:11,output[,1],ylim=c(0.9,1.1),ylab="cor(BETA,ancestry adj BETA)",xlab="GWAS",pch=19,xaxt="n",main="Effect sizes")
#axis(1,at=1:11,labels=phenoname,las=2,cex.axis=0.5)
#dev.off()

##Correlate these to the LDSC intercept values
##Take only the LDSC intercepts from ancestry regressed, non-GC corrected, surface area
LDSCint = read.csv(fLDSCint);
ind = which(LDSCint$Anc_reg==TRUE);
LDSCint = LDSCint[ind,];
##Rearrange in same order
#LDSCint = LDSCint[c(7,1:6,8:35),];

#matchedout = output[c(grep("Full_SurfArea",rownames(output)),grep("surfavg",rownames(output))),];
output_df=as.data.frame(output)
tmpcor = cor.test(output_df$corBETA,LDSCint$LDSC_intercept);
corval = tmpcor$estimate;
pval = tmpcor$p.value;
model = lm(LDSCint$LDSC_intercept ~ output_df$corBETA);
plot(output_df$corBETA,LDSCint$LDSC_intercept,xlab="cor(BETA, ancestry adj BETA)",ylab="LDSC Intercept unadjusted",main=paste0("LDSC intercept vs BETA correlation\ncor=",signif(corval,3),"; pval=",signif(pval,3)),pch=19);
x = range(output_df$corBETA);
y = model$coefficients[1] + model$coefficients[2]*x;
lines(x,y,lty=2,col="grey");
##Label some of the big ones with text
labelind = which(LDSCint$LDSC_intercept>1.02 | output_df$corBETA<0.9980);
text(output_df$corBETA[labelind],LDSCint$LDSC_intercept[labelind],rownames(output)[labelind],pos=4);

dev.off()