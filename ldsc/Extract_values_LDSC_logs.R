library(tidyverse)
library(here)

files <- list.files("/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/replication/ancreg_intercepts/surface", full.names = T, include.dirs = F, pattern = ".log")
files_wo <- list.files("/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/replication/nonancreg_intercepts/surface", full.names = T, include.dirs = F, pattern = ".log") #otherwise it takes directories too

all_files <- c(files, files_wo)

intercepts <- tibble("Region" = as.character(), 
                     "Surf_Thic" = as.character(),
                     "Anc_reg" = as.logical(),
                     "Total_h2" = as.numeric(),
                     "Total_h2_sterr" = as.numeric(),
                     "Lambda_GC" = as.numeric(),
                     "Mean_chisq" = as.numeric(),
                     "LDSC_intercept" = as.numeric(),
                     "LDSC_int_sterr" = as.numeric())

#i=1
for (i in 1:length(all_files)) {
  results <- read_lines(all_files[i], skip = 25, n_max = 4)
  values <- as.numeric(unlist(str_extract_all(string = results, 
                                              pattern = "(?<=[:punct:]{1}\\s?)[:digit:]{1}\\.[:digit:]{2,4}")))
  info1 <- str_split(all_files[i], pattern = "/")
  info2 <- str_split(info1[[1]][12], pattern = "\\.")
  
  Surf_Thic <- case_when(grepl("hick",info1[[1]][11]) ~ "Thickness",
                         grepl("surf",info1[[1]][11]) ~ "Surface Area",
                         TRUE ~ "Not Applicable")
  
  #info2[[1]][1] == "global" & grepl("hick",info2[[1]][3]) ~ "Thickness",
  #info2[[1]][1] == "global" & grepl("urf",info2[[1]][3]) ~ "Surface Area",
  #info2[[1]][2] == "global" & grepl("hick",info2[[1]][3]) ~ "Thickness",
  #info2[[1]][2] == "global" & grepl("urf",info2[[1]][3]) ~ "Surface Area",
  
  Region <- str_split(info2[[1]][1], pattern = "\\.")[[1]][1]
  #info2[[1]][2]
  #paste(info2[[1]][2],str_split(info2[[1]][3], pattern = "\\.")[[1]][1],sep="_")
  Anc_reg <- if_else(grepl("non",info1[[1]][10]), FALSE, TRUE)
  to_add <- tibble("Region" = as.character(Region),
         "Surf_Thic" = as.character(Surf_Thic),
         "Anc_reg" = as.logical(Anc_reg),
         "Total_h2" = as.numeric(values[1]),
         "Total_h2_sterr" = as.numeric(values[2]),
         "Lambda_GC" = as.numeric(values[3]),
         "Mean_chisq" = as.numeric(values[4]),
         "LDSC_intercept" = as.numeric(values[5]),
         "LDSC_int_sterr" = as.numeric(values[6]))
  intercepts <<- rbind(intercepts, to_add)
}

#regionordering = read.csv(here("scripts", "1000Genomes_Phase3_Analysis","plotting","freesurfer_orderandcolor_Height.csv"))

#intercepts$Region = factor(intercepts$Region, levels = regionordering$Region)
#intercepts$Surf_Thic = factor(intercepts$Surf_Thic)
#View(intercepts)
write_csv(intercepts, "/data/clusterfs/lag/users/gokala/enigma-evol/ldsc/replication_LDSC_intercepts_w_and_wo_ancreg.csv")
## For making initial plots

# 
# intercepts_sa <- intercepts %>% 
#   select(Region, Surf_Thic, Anc_reg, LDSC_intercept) %>% 
#   filter(Surf_Thic == "Surface Area") %>% 
#   spread(Anc_reg, LDSC_intercept)
# 
# cleveland <- intercepts_sa %>% 
#   ggplot() + 
#   geom_segment(aes(x = Region, xend = Region, y = `TRUE`,  yend = `FALSE`), color = "grey50") +
#   geom_point(aes(x = Region, y = `TRUE`), color = "green", size = 3) +
#   geom_point(aes(x = Region, y = `FALSE`), color = "purple", size = 3) +
#   labs(y = "LDSC Intercept", 
#        x = "Region", 
#        title = "Change in LDSC intercept, cortical surface area and height",
#        subtitle = "Following ancestry regression (1000 Genomes, Phase 3)\nPurple = before, green = after") +
#   coord_flip() +
#   theme_classic()
# cleveland
# 
# ggsave(filename = "LDSC_ancreg/LDSC_ancreg_before_after_SA.png", width = 7, height = 7, unit = "in")
# 
# 
# ## and for thickness
# intercepts_sa <- intercepts %>% 
#   select(Region, Surf_Thic, Anc_reg, LDSC_intercept) %>% 
#   filter(Surf_Thic == "Surface Area") %>% 
#   spread(Anc_reg, LDSC_intercept)
# # rowwise() %>% mutate( mymean = mean(c(`TRUE`, `FALSE`) )) %>% arrange(desc(mymean))
# 
# cleveland <- intercepts_sa %>% 
#   ggplot() + 
#   geom_segment(aes(x = Region, xend = Region, y = `TRUE`,  yend = `FALSE`), color = "grey50") +
#   geom_point(aes(x = Region, y = `TRUE`), color = "green", size = 3) +
#   geom_point(aes(x = Region, y = `FALSE`), color = "purple", size = 3) +
#   #facet_wrap(Surf_Thic, ncol = 2) +
#   labs(y = "LDSC intercept", 
#        x = "Trait", 
#        title = "Change in LDSC intercept, cortical surface area and height",
#        subtitle = "Following ancestry regression (1000 Genomes, Phase 3)\nPurple = before, green = after") +
#   coord_flip() +
#   theme_classic()
# cleveland
# 
# ggsave(filename = "LDSC_ancreg/LDSC_ancreg_before_after_SA.png", width = 7, height = 7, unit = "in")

#old_files_wo <- list.files(here("data","LDSC_ancreg","wo_1KGP3_ancreg"), full.names = T)

#intercepts <- tibble("Region" = as.character(), 
#                     "Surf_Thic" = as.character(),
#                     "Anc_reg" = as.logical(),
#                     "Total_h2" = as.numeric(),
#                     "Total_h2_sterr" = as.numeric(),
#                     "Lambda_GC" = as.numeric(),
#                     "Mean_chisq" = as.numeric(),
#                     "LDSC_intercept" = as.numeric(),
#                     "LDSC_int_sterr" = as.numeric())

#for (i in 1:length(old_files_wo)) {
#  results <- read_lines(old_files_wo[i], skip = 25, n_max = 4)
#  values <- as.numeric(unlist(str_extract_all(string = results, 
#                                              pattern = "(?<=[:punct:]{1}\\s?)[:digit:]{1}\\.[:digit:]{2,4}")))
#  info1 <- str_split(old_files_wo[i], pattern = "/")
#  # print(info1)
#  # if (grepl("Mean", info1[[1]][11])) {
#  if (length(info1[[1]]) == 10) {
#    info2 <- str_split(info1[[1]][10], pattern = "_")
#    # print(info2)
#    Surf_Thic <- if_else(grepl("hick",info2[[1]][3]), "Thickness", "Surface Area")
#    Region <- info2[[1]][2]
#  } else {
#    info2 <- str_split(info1[[1]][11], pattern = "_")
#    # print(info2)
#    Surf_Thic <- if_else(grepl("hick",info2[[1]][2]), "Thickness", "Surface Area")
#    Region <- info2[[1]][1]
#  }
#  Anc_reg <- if_else(grepl("wo",info1[[1]][9]), FALSE, TRUE)
#  to_add <- tibble("Region" = as.character(Region), 
#                   "Surf_Thic" = as.character(Surf_Thic),
#                   "Anc_reg" = as.logical(Anc_reg),
#                   "Total_h2" = as.numeric(values[1]),
#                   "Total_h2_sterr" = as.numeric(values[2]),
#                   "Lambda_GC" = as.numeric(values[3]),
#                   "Mean_chisq" = as.numeric(values[4]),
#                   "LDSC_intercept" = as.numeric(values[5]),
#                   "LDSC_int_sterr" = as.numeric(values[6]))
#  intercepts <<- rbind(intercepts, to_add)
#}
## Tidying up the results table
#intercepts$Region <- if_else(grepl("GIANT", intercepts$Region), "Height", intercepts$Region)
#intercepts$Region <- if_else(grepl("ancreg", intercepts$Region), "Height", intercepts$Region)
#
#regionordering = read.csv(here("scripts", "1000Genomes_Phase3_Analysis","plotting","freesurfer_orderandcolor_Height.csv"))
#
#intercepts$Region = factor(intercepts$Region, levels = regionordering$Region)
#intercepts$Surf_Thic = factor(intercepts$Surf_Thic)
#
#write_csv(intercepts, here("data", "LDSC_ancreg", "LDSC_intercepts_before_after_Phase3_ancreg_old_version.csv"))