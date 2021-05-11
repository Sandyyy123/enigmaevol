# ========================================================
# Making a barplot of the LDSC partitioned heritability
# "proportion of SNP heritability explained" for a single
# annotation. The human gained enhancers data is what
# I'm plotting here, for manuscript Figure 3B.


# ancreg_MAe3ukw3/
# |__ ancreg_Phase3_results/
# |   |__ chimp_PFC_enhancers_hg19/
# |       |__ Mean_Full_SurfArea_ancreg_1KGP3_ancreg.txt.sumstats.gz.results
# |       |__ Mean_Full_Thickness_ancreg_1KGP3_ancreg.txt.sumstats.gz.results
# |   |__ chimp_PFC_promoters_hg1s9/
# |__ results_tables/***
# |__ plots/

# updated for the nonGC version

# ========================================================

library(tidyverse)
library(scales)
library(BSDA)

options(stringsAsFactors=FALSE)

#regionordering <- read.csv("/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/Cerebral_Cortex_revisions/plotting/freesurfer_orderandcolor.csv")
annots = list.files(path = "/data/clusterfs/lag/users/gokala/enigma-evol/ancreg/baseline_subset/replication_v1/munged/results_tables/", full.names = F,recursive = F)

for (i in 1:length(annots)){
print(annots[i])
annotname=str_split(annots[i], pattern = "\\.")[[1]][1]
resultsFile <- Sys.glob(paste0("/data/clusterfs/lag/users/gokala/enigma-evol/ancreg/baseline_subset/replication_v1/munged/results_tables/",annots[i]))
results <- read.table(resultsFile[1], header = TRUE, sep = "\t")
results$Region <- factor(results$Region)
#, levels = regionordering$Region
y_max <- max(results[results$Analysis == "Surface Area", 3])
y_axis_max <- y_max + results[results$Prop._h2 == y_max, 4] + 0.02
label.df <- data.frame(Region = results$Region[results$significant=="Yes"],Prop._h2=results$Prop._h2[results$significant=="Yes"]+0.05)

# Version where the error bars don't push the axis below 0
pSA <- ggplot(data = results[results$Analysis == "Surface Area", ], mapping = aes(Region, Prop._h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "burlywood4") +
  #scale_y_continuous(limits = c(0, y_axis_max), oob = squish) +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Prop._h2 - Prop._h2_std_error, ymax = Prop._h2 + Prop._h2_std_error)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Region",
    y = expression(paste("Proportion ", italic("h"^{2}))),
    title = "Partitioned heritability for cortical surface area, following ancestry regression",
    subtitle = "Proportion of heritability explained"
  )
pSA + geom_text(data = label.df, label = "*")
ggsave(paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/",annotname,"_h2_Prop_SA_ancreg_Phase3_squished.pdf"), width = 6.5, height = 3.25, units = "in")
}

# Version where error bars do go below 0
pSA <- ggplot(data = results[results$Analysis == "Surface Area", ], mapping = aes(Region, Prop._h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "burlywood4") +
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Prop._h2 - Prop._h2_std_error, ymax = Prop._h2 + Prop._h2_std_error)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Region",
    y = expression(paste("Prop. ", italic("h"^{
      2
    }))),
    title = "Partitioned heritability for cortical surface area, following ancestry regression",
    subtitle = "Proportion of heritability explained"
  )
pSA + geom_text(data = label.df, label = "*")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/plots/bothDirections_globalCov_L-R_HSE_7pcw_active_merged_Prop_h2_SA_MA6_ancreg_Phase3.pdf", width = 10, height = 3, units = "in")