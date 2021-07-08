## ----------------------------------------------------
# Partitioned Heritability brainplots for ENIGMA-EVO manuscript
# for HGE 7pcw only, using ancestry regressed, noGC data
## ----------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(diffloop)
library(plotly)

# The functions for plotly plots are here:
#source("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/test_brain_plot.R")
source("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/partherit/plotly_brainplot_functions.R")

# Reading the partitioned heritability summary statistics
HSE_7pcw <- read.delim("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/regional_hemi_spec_glob/HSE_7pcw_active_merged_results_FDR35.txt", header = TRUE)
HSE_7pcw[HSE_7pcw$Region=="global",]$Region="Full"

HSE_7pcw$Region <- factor(HSE_7pcw$Region,levels=regionordering$Region)

factor(HSE_7pcw$Region,levels=regionordering$Region)

HSE_7pcw_SA <- HSE_7pcw %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_TH <- HSE_7pcw %>% dplyr::filter(Analysis == "Thickness")
#rownames(HSE_7pcw_TH)=rownames(HSE_7pcw_SA)

enrich_maxSA <- max(HSE_7pcw_SA$Enrichment)
enrich_maxTH <- max(HSE_7pcw_TH$Enrichment)

## If the FDR-corrected pvalue is > 0.05, set the Enrichment to 0 so it's grey in the brain plot
HSE_7pcw_SA$Enrichment_plot <- if_else(HSE_7pcw_SA$fdr >= 0.05, true = 0, false = HSE_7pcw_SA$Enrichment)
HSE_7pcw_TH$Enrichment_plot <- if_else(HSE_7pcw_TH$fdr >= 0.05, true = 0, false = HSE_7pcw_TH$Enrichment)

brainplot_full(dta = HSE_7pcw_SA,
               out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/replication_whiteBrit_ancreg/mpi_colors/",
               out_suffix = "nonancreg_whiteBrit",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_maxSA,
               low_color = "#ED6B06",
               high_color = "#00786A", 
               nonsig_color = "#BABABC")

brainplot_SA(dta = HSE_7pcw_SA,
             out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/replication_whiteBrit_ancreg/mpi_colors/",
             out_suffix = "ancreg_whiteBrit",
             loop_over = "Annotation",
             region_col = "Region",
             Z_col = "Enrichment_plot",
             analysis_col = "Analysis",
             max_val = enrich_maxSA,
             low_color = "#ED6B06",
             high_color = "#00786A",
             nonsig_color = "#BABABC")


brainplot_TH(dta = HSE_7pcw_TH,
             out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/replication_whiteBrit/",
             out_suffix = "nonancreg_EUR",
             loop_over = "Annotation",
             region_col = "Region",
             Z_col = "Enrichment",
             analysis_col = "Analysis",
             max_val = enrich_maxTH,
             low_color = "#ED6B06",
             high_color = "#00786A")

############################################
# Same plots for L-R partitioned h2 results#
############################################

HSE_7pcw_le = HSE_7pcw[grep("_le", HSE_7pcw$Region),]
HSE_7pcw_re = HSE_7pcw[grep("_re", HSE_7pcw$Region),]
HSE_7pcw_le[] <- lapply(HSE_7pcw_le, function(x) sub("_le", "", x, fixed = TRUE))
HSE_7pcw_re[] <- lapply(HSE_7pcw_re, function(x) sub("_re", "", x, fixed = TRUE))

HSE_7pcw_le$Region <- factor(HSE_7pcw_le$Region,levels=regionordering$Region)
HSE_7pcw_re$Region <- factor(HSE_7pcw_re$Region,levels=regionordering$Region)

HSE_7pcw_SA_le <- HSE_7pcw_le %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_SA_re <- HSE_7pcw_re %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_TH_le <- HSE_7pcw_le %>% dplyr::filter(Analysis == "Thickness")
HSE_7pcw_TH_re <- HSE_7pcw_re %>% dplyr::filter(Analysis == "Thickness")

HSE_7pcw_SA_re$Enrichment = as.double(HSE_7pcw_SA_re$Enrichment)
enrich_max <- max(HSE_7pcw_SA_re$Enrichment) #because we want all the plots on this figure to be using the same color axis

HSE_7pcw_SA_re$Enrichment_plot <- if_else(HSE_7pcw_SA_re$fdr >= 0.05, true = 0, false = HSE_7pcw_SA_re$Enrichment)

brainplot_full(dta = HSE_7pcw_SA_le,
               out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/l-r_hemspecCov/mpi_colors/",
               out_suffix = "left_ancreg_Phase3_noGC",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_max,
               low_color = "midnightblue",
               high_color = "darkred", 
               nonsig_color = "#BABABC")

brainplot_SA(dta = HSE_7pcw_SA_re,
               out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/l-r_hemspecCov/mpi_colors/",
               out_suffix = "right_ancreg_Phase3_noGC",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_max,
               low_color = "#ED6B06",
               high_color = "#00786A",
               nonsig_color = "#BABABC")


brainplot_TH(dta = HSE_7pcw_TH_le,
             out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/plots/partherit_brainplots/l-r_hemspecCov/mpi_colors/",
             out_suffix = "^_left_ancreg_Phase3_noGC",
             loop_over = "Annotation",
             region_col = "Region",
             Z_col = "Enrichment",
             analysis_col = "Analysis",
             max_val = enrich_max,
             low_color = "midnightblue",
             high_color = "darkred")

################

# Reading the partitioned heritability summary statistics
HSE_7pcw <- read.delim("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/regional_with_global/NeanDepleted_results_FDR35.txt", header = TRUE)
HSE_7pcw_le = HSE_7pcw[grep(" le _", HSE_7pcw$Region),]
HSE_7pcw_re = HSE_7pcw[grep(" re _", HSE_7pcw$Region),]
HSE_7pcw_le[] <- lapply(HSE_7pcw_le, function(x) sub(" le _", "", x, fixed = TRUE))
HSE_7pcw_re[] <- lapply(HSE_7pcw_re, function(x) sub(" re _", "", x, fixed = TRUE))
#HSE_7pcw$Region <- factor(HSE_7pcw$Region,levels=regionordering$Region)
HSE_7pcw_le$Region <- factor(HSE_7pcw_le$Region,levels=regionordering$Region)
HSE_7pcw_re$Region <- factor(HSE_7pcw_re$Region,levels=regionordering$Region)

#factor(HSE_7pcw$Region,levels=regionordering$Region)
#HSE_7pcw_SA <- HSE_7pcw %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_SA_le <- HSE_7pcw_le %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_SA_le$Region[34]="Full"
HSE_7pcw_SA_re <- HSE_7pcw_re %>% dplyr::filter(Analysis == "Surface Area")
HSE_7pcw_SA_re$Region[34]="Full"
#HSE_7pcw_TH <- HSE_7pcw %>% dplyr::filter(Analysis == "Thickness")
#HSE_7pcw_TH_le <- HSE_7pcw_le %>% dplyr::filter(Analysis == "Thickness")
#HSE_7pcw_TH_re <- HSE_7pcw_re %>% dplyr::filter(Analysis == "Thickness")
#enrich_max <- max(HSE_7pcw_SA$Enrichment)
HSE_7pcw_SA_re$Enrichment = as.double(HSE_7pcw_SA_re$Enrichment)
enrich_max <- max(HSE_7pcw_SA_re$Enrichment) #because we want all the plots on this figure to be using the same color axis
## If the FDR-corrected pvalue is > 0.05, set the Enrichment to 0 so it's grey in the brain plot
#HSE_7pcw_SA$Enrichment_plot <- if_else(HSE_7pcw_SA$fdr >= 0.05, true = 0, false = HSE_7pcw_SA$Enrichment)
HSE_7pcw_SA_re$Enrichment_plot <- if_else(HSE_7pcw_SA_re$fdr >= 0.05, true = 0, false = HSE_7pcw_SA_re$Enrichment)
#P:/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/

brainplot_full(dta = HSE_7pcw_SA_re,
               out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/partherit/plots/partherit_brainplots/l-r_globalCov/NeanDep",
               out_suffix = "_right_ancreg_Phase3_noGC",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_max,
               low_color = "midnightblue",
               high_color = "darkred", 
               nonsig_color = "#BABABC")

brainplot_SA(dta = HSE_7pcw_SA_re,
             out_prefix = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/partherit/plots/partherit_brainplots/l-r_globalCov/",
             out_suffix = "_right_ancreg_Phase3_noGC",
             loop_over = "Annotation",
             region_col = "Region",
             Z_col = "Enrichment_plot",
             analysis_col = "Analysis",
             max_val = enrich_max,
             low_color = "midnightblue",
             high_color = "darkred",
             nonsig_color = "#BABABC")
