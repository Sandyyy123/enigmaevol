library(kableExtra)
library(ggplot2)
##

mi=read.csv("/data/clusterfs/lag/users/gokala/enigma-evol/sds_test/sds/SDS_bjk_milk_intake_ancreg_1kblocks.csv",header=T)
lf=read.csv("/data/clusterfs/lag/users/gokala/enigma-evol/sds_test/sds/SDS_bjk_lactose_free_diet_ancreg_1kblocks.csv",header=T)
sm=read.csv("/data/clusterfs/lag/users/gokala/enigma-evol/sds_test/sds/SDS_bjk_spontaneousMiscarriage_stillBirth_termination_ancreg_1kblocks.csv",header=T)
bw=read.csv("/data/clusterfs/lag/users/gokala/enigma-evol/sds_test/sds/SDS_bjk_birth_weight_ancreg_1kblocks.csv",header=T)

df=rbind(mi,lf,sm,bw)
colnames(df)[1]=''
rownames(df)=df[,1]
df=df[,-1]

caption = "SDS results for test sumstats"
kable_classic(full_width = F, html_font = "Cambria")

df %>%
  kbl(caption = "SDS results of test sumstats from the Neale group - UKB phenotypes") %>%
  kable_classic(full_width = F, html_font = "Cambria")

  save_kable("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/sds_test.pdf")