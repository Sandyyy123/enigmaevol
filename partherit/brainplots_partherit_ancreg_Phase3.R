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
library(here)

# The functions for plotly plots are here:
source(here("./plotting", "plotly_brainplot_functions.R"))

# Reading the partitioned heritability summary statistics
HSE_7pcw <- read.delim("./partherit/results_tables/HSE_7pcw_active_merged_results_FDR35.txt", header = TRUE)
HSE_7pcw$Region <- factor(HSE_7pcw$Region,levels=regionordering$Region)
HSE_7pcw_SA <- HSE_7pcw %>% dplyr::filter(Analysis == "Surface Area")
enrich_max <- max(HSE_7pcw_SA$Enrichment) #because we want all the plots on this figure to be using the same color axis


## If the FDR-corrected pvalue is > 0.05, set the Enrichment to 0 so it's grey in the brain plot
HSE_7pcw_SA$Enrichment_plot <- if_else(HSE_7pcw_SA$fdr >= 0.05, true = 0, false = HSE_7pcw_SA$Enrichment)

brainplot_full(dta = HSE_7pcw_SA,
               out_prefix = "P:/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/partherit/plots/partherit_brainplots/",
               out_suffix = "MA6_ancreg_Phase3_noGC",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_max,
               low_color = "midnightblue",
               high_color = "darkred", 
               nonsig_color = "#BABABC")

brainplot_SA(dta = HSE_7pcw_SA,
               out_prefix = "P:/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/partherit/plots/partherit_brainplots/",
               out_suffix = "MA6_ancreg_Phase3_noGC",
               loop_over = "Annotation",
               region_col = "Region",
               Z_col = "Enrichment_plot",
               analysis_col = "Analysis",
               max_val = enrich_max,
               low_color = "midnightblue",
               high_color = "darkred",
               nonsig_color = "#BABABC")


brainplot_TH(dta = HSE_7pcw,
             out_prefix = "P:/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/partherit/plots/partherit_brainplots/",
             out_suffix = "MA^_ancreg_Phase3_noGC",
             loop_over = "Annotation",
             region_col = "Region",
             Z_col = "Enrichment",
             analysis_col = "Analysis",
             max_val = enrich_max,
             low_color = "midnightblue",
             high_color = "darkred")
