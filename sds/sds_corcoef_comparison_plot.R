
options(stringsAsFactors=FALSE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)

fcorvals = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/sds/SDS_BJK_ancreg_replication_1kblocks.txt"
#fcorvals_eur = "/data/clusterfs/lag/users/gokala/enigma-evol/sds/SDS_bjk_ancreg_1kblocks_replication2.csv"
fjason_corvals = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/sds/SDS_bjk_ancreg_1kblocks_fromJason.csv"
corvals = read.table(fcorvals,header = T)
#corvals_eur = read.csv(fcorvals_eur)
jason_corvals = read.csv(fjason_corvals)

hist(jason_corvals$BJK_ESTIM_AVE)
hist(as.numeric(corvals$BJK_ESTIM_AVE))

## Pick relevant rows and columns from raw SDS outputs

# Prep E3 results
jason_corvals = jason_corvals[grepl("(?=.*_surfavg)|(?=.*Full)(?=.*Surf)",jason_corvals$X,perl=T),]
jason_corvals=jason_corvals[,c(1,3,6)]
jason_corvals=jason_corvals[jason_corvals$X!="Mean_temporalpole_surfavg",] # remove temporal lobe, because it's not included in the replication
jason_corvals$mergeCol = sapply(jason_corvals$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

# Prep White British replication results

thisSA1 = corvals[grep("surf",corvals$X),]
thisSA = corvals[grep("withGlob",thisSA1$X),]
thisSA[nrow(thisSA)+1,] = corvals[grep("Full",thisSA1$X),]
brit_surface = thisSA[,c(1,3,6)]
brit_surface$BJK_ESTIM_AVE = as.numeric(brit_surface$BJK_ESTIM_AVE)
#brit_surface = corvals %>%
#  filter(grepl("(?=.*surface)(?=.*_withGlob)",corvals$X,perl=TRUE)) %>%
#  select(X, BJK_ESTIM_AVE, BJK_ESTIM_PVAL)

#brit_surface$X[8]="Mean_Full_Surface"
#brit_surface=brit_surface[-9,] #remove globScaled row, because it is replicate of Mean_Full_Surface
brit_surface$mergeCol=sapply(brit_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
e3_brit = merge(jason_corvals,brit_surface,by="mergeCol")

#brit_thickness = corvals %>%
#  filter(grepl("(?=.*thicknessDK)(?=.*globalCov)|(?=.*thickness)(?=.*global)",corvals$X,perl=TRUE)) %>%
#  select(X, global_corr_spearman)
#brit_thickness$mergeCol=sapply(brit_thickness$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

# Prep Eur replication results
#eur_surface = corvals_eur %>%
#  filter(grepl("(?=.*surfaceDK)",corvals_eur$X,perl=TRUE)) %>%
#  select(X, BJK_ESTIM_Z, BJK_ESTIM_PVAL)

#eur_surface$X[8]="Mean_Full_Surface"
#eur_surface=eur_surface[-9,] #remove globScaled row, because it is replicate of Mean_Full_Surface
#eur_surface$mergeCol=sapply(eur_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
#e3_eur = merge(jason_corvals,eur_surface,by="mergeCol")

# Correlation plots
ggscatter(e3_brit, x = "BJK_ESTIM_AVE.x", y = "BJK_ESTIM_AVE.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "E3", ylab = "White British subset", title = "SDS estimated average corcoef. comparison (surface area)")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/sds/surface_e3_vs_replication_scatter_corCoef.pdf",width=7,height=7)

ggscatter(e3_eur, x = "BJK_ESTIM_Z.x", y = "BJK_ESTIM_Z.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "E3", ylab = "European subset", title = "Z-score comparison (surface area)")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/plots/sds/surface_e3_vs_eurSub_scatter.pdf",width=7,height=7)

#pdf("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/thickness_e3_vs_replication_scatter2.pdf")
#plot(merged$global_corr_spearman.x,merged$global_corr_spearman.y,xlab="ENIGMA E3",ylab="Replication")
#dev.off()

#mergedTall = merged %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
#ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("ENIGMA MA3", "Replication")) + ggtitle("Comparison of Thickness SDS Correlation Coefficients in MA3 and UKB replication") + ylab("Correlation Coefficient") + xlab("Region")
#ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/thickness_corcoef_comparison.pdf")

#### Read-in left-right SDS corcoef.s

fcorvals_lr = "/data/clusterfs/lag/users/gokala/enigma-evol/sds/SDS_bjk_ancreg_1kblocks_regional_with_global.csv"
corvals_lr = read.csv(fcorvals_lr,stringsAsFactors=F)

#replication_surface = corvals %>%
#  filter(grepl("(?=.*surface)(?=.*globalCov)|(?=.*surface)(?=.*global)",corvals$X,perl=TRUE)) %>%
#  select(X, global_corr_spearman)
#replication_surface$mergeCol=sapply(replication_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

left = corvals_lr %>%
  filter(grepl("(?=.*le_)",corvals_lr$X,perl=TRUE)) %>%
  select(X, BJK_ESTIM_Z)
left$mergeCol=sapply(left$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
left$sample="left"
merged_left=merge(jason_corvals,left,by="mergeCol")

ggscatter(merged_left, x = "BJK_ESTIM_Z.x", y = "BJK_ESTIM_Z.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "E3", ylab = "Left hemisphere", title = "Z-score comparison (surface area)")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/plots/sds/surface_e3_vs_left_scatter.pdf",width=7,height=7)

right = corvals_lr %>%
  filter(grepl("(?=.*re_)",corvals_lr$X,perl=TRUE)) %>%
  select(X, BJK_ESTIM_Z)
right$mergeCol=sapply(right$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
right$sample="right"
merged_right=merge(jason_corvals,right,by="mergeCol")

ggscatter(merged_right, x = "BJK_ESTIM_Z.x", y = "BJK_ESTIM_Z.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "E3", ylab = "Right hemisphere", title = "Z-score comparison (surface area)")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/plots/sds/surface_e3_vs_right_scatter.pdf",width=7,height=7)

#####

#mergedTallLeft = merged_left %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
#comparisonplot=ggplot(mergedTallLeft, aes(mergeCol, CorCoef, fill = Data))+ coord_flip()  + geom_col(position = "dodge") + scale_fill_discrete(name = "Source", labels = c("UKB replication", "Left Hemisphere")) + ggtitle("SDS Cor. Coef.s of UKB Replication Averaged vs. Left Hemisphere GWASes") + ylab("Correlation Coefficient") + xlab("Region")
#ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_replication_vs_left.pdf")

#mergedTallRight = merged_right %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
#comparisonplot=ggplot(mergedTallRight, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("UKB replication", "Right Hemisphere")) + ggtitle("SDS Cor. Coef.s of UKB Replication Averaged vs. Right Hemisphere GWASes") + ylab("Correlation Coefficient") + xlab("Region")
#ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_replication_vs_right.pdf")

merged_1=Reduce(function(...) merge(..., by="mergeCol",all=TRUE), list(jason_corvals,brit_surface,eur_surface))
merged_1=merged_1[!duplicated(merged_1$mergeCol),]
merged_1=merged_1[,c(1,3,6,9)]
colnames(merged_1)=c("mergeCol","E3","Brit","Eur")
merged_1.long = pivot_longer(merged_1, cols=2:4, names_to = "Sample", values_to = "Zscore")
merged_1.long = merged_1.long[c(which(merged_1.long$mergeCol=="Full"),which(merged_1.long$mergeCol!="Full")),] # move Full to the top
#rects = data.frame(ystart = rep(1:34, each=3), yend = rep(1:34, each=3),col=c(0,1))

#merged_1.long$mergeCol <- factor(merged_1.long$mergeCol, levels=unique(merged_1.long$mergeCol))
merged_2=Reduce(function(...) merge(..., by="mergeCol",all=TRUE), list(jason_corvals,left,right))
merged_2=merged_2[!duplicated(merged_2$mergeCol),]
merged_2=merged_2[,c(1,3,6,9)]
colnames(merged_2)=c("mergeCol","E3","Left","Right")
merged_2.long = pivot_longer(merged_2, cols=2:4, names_to = "Sample", values_to = "Zscore")
merged_2.long = merged_2.long[merged_2.long$mergeCol!="Full",]

ggplot(merged_1.long, aes(fct_rev(merged_1.long$mergeCol), merged_1.long$Zscore, fill = fct_rev(merged_1.long$Sample))) +
  geom_col(position=position_dodge()) +
  coord_flip() + ggtitle("Block-jackknife correlation") +
  ylab("Ancestry adjusted Z-scores") + xlab("Region") +
  scale_x_discrete(expand = c(0,0.5)) +
  scale_fill_manual(name = "Source", labels = c("European subset", "E3", "White British subset"), values = c("#E69F00", "#D55E00", "#009E73")) +
  theme_gray() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        panel.border = element_blank())
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_e3_brit_eur.pdf")

ggplot(merged_2.long, aes(fct_rev(mergeCol), Zscore, fill = fct_rev(Sample))) + 
  geom_col(position=position_dodge()) +
  coord_flip() + ggtitle("Block-jackknife correlation") +
  ylab("Ancestry adjusted Z-scores") + xlab("Region") +
  scale_fill_manual(name = "Source", labels = c("Right hemisphere", "Left hemisphere", "E3"), values = c("#F0E442", "#0072B2", "#D55E00")) +
  theme_gray() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        panel.border = element_blank())
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_e3_left_right.pdf")