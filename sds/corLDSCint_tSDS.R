##/usr/local/R-3.2.3/bin/R
##Correlate the LDSC intercepts (measure of pop stratification) with the measure of selection
##cor(tSDS,beta_r) to determine if residual ancestry is impacting 
options(stringsAsFactors=FALSE);
library(GenomicRanges);

##The relationships between a measure of selection and genetic influences on cortical structure (cor(tSDS),beta_r)
fSDS = "/data/clusterfs/lag/users/gokala/genlang-evol/sds/SDS_bjk_ancreg_1kblocks.csv";
##LDSC intercept values from original ENIGMA (not ancestry regressed)
fLDSCint = "/data/clusterfs/lag/users/gokala/genlang-evol/LDSC/LDSC_intercepts_w_and_wo_ancreg.csv";

##Read in the SDS files
SDS = read.csv(fSDS,row.names=1);

##Correlate these to the LDSC intercept values
##Take only the LDSC intercepts from ancestry regressed, non-GC corrected, surface area
LDSCint = read.csv(fLDSCint);
ind = which(LDSCint$Surf_Thic=="Surface Area" & LDSCint$Anc_reg==FALSE);
LDSCint = LDSCint[ind,];
##Rearrange in same order
LDSCint = LDSCint[c(7,1:6,8:35),];

##Organize the two in the same order, only take surface area
matchedLDSC = LDSCint;
matchedSDS = SDS[c(grep("Full_SurfArea",rownames(SDS)),grep("surfavg",rownames(SDS))),];

##Determine if there is a relationship between the two
corout = cor.test(abs(matchedSDS$BJK_ESTIM_Z),matchedLDSC$LDSC_intercept);
corval = corout$estimate;
corpval = corout$p.value;
model = lm(matchedLDSC$LDSC_intercept ~ abs(matchedSDS$BJK_ESTIM_Z));

##Plot the relationship
pdf("SDS_E3_VS_LDSCint_wo_anc_reg.pdf");
plot(abs(matchedSDS$BJK_ESTIM_Z),matchedLDSC$LDSC_intercept,xlab="Association between tSDS and Cortical SA GWAS prior to ancestry regression\nabs(Z-score)",ylab="LDSC Intercept prior to ancestry regression",main=paste0("The impact of subtle ancestry regression on selection measured with tSDS\ncor=",signif(corval,3),"; P=",signif(corpval,3)),pch=19);
x = range(abs(matchedSDS$BJK_ESTIM_Z));
y = model$coefficients[1] + model$coefficients[2]*x;
lines(x,y,lty=2,col="grey");
##Label some of the big ones with text
##labelind = which(matchedLDSC$LDSC_intercept_wo_anc_reg>1.02);
##text(abs(matchedSDS$BJK_ESTIM_Z)[labelind],matchedLDSC$LDSC_intercept_wo_anc_reg[labelind],rownames(matchedSDS)[labelind],pos=4);

dev.off();





