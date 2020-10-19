## Making a Cleveland dotplot for ENIGMA_Evolution Figure 2
library(tidyverse)
library(here)


#regionordering = read.csv("plotting/freesurfer_orderandcolor.csv")
intercepts <- read_csv("/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/LDSC_intercepts_w_and_wo_ancreg.csv", col_names = TRUE)
intercepts_surf <- intercepts[intercepts$Surf_Thic == "Surface Area",]
intercepts_thic <- intercepts[intercepts$Surf_Thic == "Thickness",]
#intercepts$Region = factor(intercepts$Region, levels = regionordering$Region)
#intercepts$Region <- fct_rev(intercepts$Region)
View(intercepts_surf)
## Plot for surface area

intercepts_sterr_surf <- intercepts_surf %>% 
  select(Region, Anc_reg, LDSC_int_sterr) %>% 
  spread(Anc_reg, LDSC_int_sterr)
  
intercepts_wide_surf <- intercepts_surf %>% 
  select(Region, Anc_reg, LDSC_intercept) %>% 
  spread(Anc_reg, LDSC_intercept)

plot_data_surf <- left_join(intercepts_wide_surf, intercepts_sterr_surf, by = c("Region"), suffix = c(".intercept", ".sterr"))

Surf_plot_surf <- intercepts_surf %>% 
  filter(Surf_Thic == "Surface Area") %>% 
  ggplot(aes(x = Region, y = LDSC_intercept, color = Anc_reg, group = Anc_reg)) + 
  geom_errorbar(aes(ymax = LDSC_intercept + LDSC_int_sterr,  ymin = LDSC_intercept - LDSC_int_sterr), position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size = 2)+
  scale_color_manual(values=c("#af8dc3", "#7fbf7b"))+
  labs(y = "LDSC Intercept", 
       x = "Region", 
       title = "Surface Area, Change in LDSC intercept",
       subtitle = "Purple = before, green = after") +
  coord_flip() +
  theme_classic()
Surf_plot_surf

ggsave("/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/LDSC_ancreg_Rdata_before_after_w_errorbars_surfaceArea.pdf", width = 7, height = 9, unit = "in")

## Plot for thickness
View(intercepts_sterr_thic)
intercepts_sterr_thic <- intercepts_thic %>% 
  select(Region, Anc_reg, LDSC_int_sterr) %>% 
  spread(Anc_reg, LDSC_int_sterr)

intercepts_wide_thic <- intercepts_thic %>% 
  select(Region, Anc_reg, LDSC_intercept) %>% 
  spread(Anc_reg, LDSC_intercept)

plot_data_thic <- left_join(intercepts_wide_thic, intercepts_sterr_thic, by = c("Region"), suffix = c(".intercept", ".sterr"))

Surf_plot_thic <- intercepts_thic %>% 
  filter(Surf_Thic == "Thickness") %>% 
  ggplot(aes(x = Region, y = LDSC_intercept, color = Anc_reg, group = Anc_reg)) + 
  geom_errorbar(aes(ymax = LDSC_intercept + LDSC_int_sterr,  ymin = LDSC_intercept - LDSC_int_sterr), position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size = 2)+
  scale_color_manual(values=c("#af8dc3", "#7fbf7b"))+
  labs(y = "LDSC Intercept", 
       x = "Region", 
       title = "Thickness, Change in LDSC intercept",
       subtitle = "Purple = before, green = after") +
  coord_flip() +
  theme_classic()
Surf_plot_thic

ggsave("/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/LDSC_ancreg_Rdata_before_after_w_errorbars_thickness.pdf", width = 7, height = 9, unit = "in")




#old_before_intercepts <- read.csv(here("data", "LDSC_ancreg", "LDSC_intercepts_before_after_Phase3_ancreg_old_version.csv"))
################################

intercepts <- intercepts[!intercepts$Region == "Height",]

# Reversing order of the brain regions for the plot, so Full is on top after coord_flip
regionordering = read.csv("freesurfer_orderandcolor.csv")
intercepts$Region = factor(intercepts$Region, levels = regionordering$Region)
intercepts$Region <- fct_rev(intercepts$Region)


reading_traits_plot <- intercepts %>% 
  select(Phenotype, Anc_reg, LDSC_intercept) %>% 
  #filter(Surf_Thic == "Surface Area") %>% 
  spread(Anc_reg, LDSC_intercept) %>% 
  ggplot() + 
  geom_segment(aes(x = Phenotype, xend = Phenotype, y = `TRUE`,  yend = `FALSE`), color = "grey50") +
  geom_point(aes(x = Phenotype, y = `TRUE`), color = "#7fbf7b", size = 3) +
  geom_point(aes(x = Phenotype, y = `FALSE`), color = "#af8dc3", size = 3) +
  labs(y = "LDSC Intercept", 
       x = "Phenotype", 
       title = "Change in LDSC intercept",
       subtitle = "Purple = before, green = after") +
  coord_flip() +
  theme_classic()
Surf_plot
ggsave("P:/workspaces/lg-genlang/Working/Evolution/LDSC_ancreg_Rdata_before_after_dotplot_grouped.pdf", width = 7, height = 9, unit = "in")