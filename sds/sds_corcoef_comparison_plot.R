library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

fcorvals = "/data/clusterfs/lag/users/gokala/enigma-evol/sds/SDS_bjk_ancreg_1kblocks.csv"
fjason_corvals = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/SDS_bjk_ancreg_1kblocks_thickness_fromJason.csv"
corvals = read.csv(fcorvals,stringsAsFactors=F)
jason_corvals = read.csv(fjason_corvals,stringsAsFactors=F)

replication_thickness = corvals %>%
  filter(grepl("(?=.*thicknessDK)(?=.*globalCov)|(?=.*Thickness)(?=.*globalValues)",corvals$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)

replace(replication_thickness$X,1,"Mean_Full_Thickness")
replication_thickness$mergeCol=sapply(replication_thickness$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
jason_corvals$mergeCol=sapply(jason_corvals$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
jason_corvals=jason_corvals[,c(1,2,7)]
merged = merge(jason_corvals,replication_thickness,by="mergeCol")

# Correlation plot
ggscatter(merged, x = "global_corr_spearman.x", y = "global_corr_spearman.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "E3", ylab = "Replication", title = "SDS correlation coefficient comparison (thickness)")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/thickness_e3_vs_replication_scatter.pdf")

pdf("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/thickness_e3_vs_replication_scatter2.pdf")
plot(merged$global_corr_spearman.x,merged$global_corr_spearman.y,xlab="ENIGMA E3",ylab="Replication")
dev.off()

mergedTall = merged %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
#ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("ENIGMA MA3", "Replication")) + ggtitle("Comparison of Thickness SDS Correlation Coefficients in MA3 and UKB replication") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/thickness_corcoef_comparison.pdf")

#### Read-in left-right SDS corcoef.s

fcorvals_lr = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/scripts/SDS_bjk_non-ancreg_left_right_withPadj_surface.csv"
corvals_lr = read.csv(fcorvals_lr,stringsAsFactors=F)

replication_surface = corvals %>%
  filter(grepl("(?=.*surface)(?=.*globalCov)|(?=.*surface)(?=.*global)",corvals$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)
replication_surface$mergeCol=sapply(replication_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

left = corvals_lr %>%
  filter(grepl("(?=.*_le)",corvals_lr$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)
left$mergeCol=sapply(left$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[3]})
merged_left=merge(replication_surface,left,by="mergeCol")

right = corvals_lr %>%
  filter(grepl("(?=.*_re)",corvals_lr$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)
right$mergeCol=sapply(right$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[3]})
merged_right=merge(replication_surface,right,by="mergeCol")

mergedTallLeft = merged_left %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
comparisonplot=ggplot(mergedTallLeft, aes(mergeCol, CorCoef, fill = Data))+ coord_flip()  + geom_col(position = "dodge") + scale_fill_discrete(name = "Source", labels = c("UKB replication", "Left Hemisphere")) + ggtitle("SDS Cor. Coef.s of UKB Replication Averaged vs. Left Hemisphere GWASes") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_replication_vs_left.pdf")

mergedTallRight = merged_right %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
comparisonplot=ggplot(mergedTallRight, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("UKB replication", "Right Hemisphere")) + ggtitle("SDS Cor. Coef.s of UKB Replication Averaged vs. Right Hemisphere GWASes") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_corcoef_replication_vs_right.pdf")