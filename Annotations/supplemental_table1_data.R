#### load libraries
library(tidyverse)
setwd("/data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/MA6/supplemental_table1")


#### load files
# bed <- read_delim("./beds/allAnnots_MAe3ukw3.bed", delim = "\t", col_names = c("CHR", "START", "END", "annot"))
# frqs <- read_table2("./1KGP3.frq") %>% select(-NCHROBS)
extras <- read_table2("1kg_phase3_ns.allpop.unrel2261_eigenvec.P1to20_beta_se_pval.txt") %>% 
  select(CHR, SNP, POS, A1, A2)

#### prep files for use: need to add the bp position via the 'extras' file 
# frqs_bed <- frqs %>% 
#   left_join(extras) %>% 
#   mutate(START = POS, END = POS) %>% 
# select(CHR, START, END, SNP, MAF)
# 
# write_delim(frqs_bed, "1KGP3_snps_w_maf.bed")
# 
# annot_list <- unique(bed$annot)
# 
# #### Note: now go and run bedtools with this command:
# 
# ## Why bedtools? because I can be sure that it will do the overlaps correctly, that's all.
# 
# /data/workspaces/lag/workspaces/lg-neanderthals/raw_data/ENIGMA-EVO/software/bedtools intersect -wo -a 1kg_phase3_ns.allpop.unrel2261_eigenvec.P1to20_beta_se_pval.txt -b ./beds/allAnnots_MAe3ukw3.bed > lkg_phase3_evo_annot_overlaps_from_bedtools.bed
# 
# #### functions
# get_annot_snps <- function(which_annot, frq_file) {
#   mini_bed <- bed %>% 
#     filter(annot == which_annot)
#   
#   the_snps <- frq_file %>% 
#     filter()
# }

#### calculate average MAF for each annot
frq_files <- Sys.glob("1000G.*.frq")

calc_avg_maf <- function(frq) {
  
  print(frq)
  
  frq_file <- read_table2(frq) %>% 
    filter(SNP %in% extras$SNP)
  
  avg_MAF <- mean(frq_file$MAF)
  
  avg_MAF
  
}

calc_nsnps <- function(frq) {
  
  frq_file <- read_table2(frq) %>% 
    filter(SNP %in% extras$SNP)
  
  nSNPs <- nrow(frq_file)
  
  nSNPs
  
}

supp1 <- tibble(annot = frq_files,
                avg_maf = map(frq_files, calc_avg_maf),
                nsnps = map(frq_files, calc_nsnps)
)

# Needed to write an extra column onto this bed file for PLINK
read_table2("./beds/NeanSNPs_3col.bed", col_names = FALSE) %>% 
  mutate(placeholder = "placeholder") %>% 
  write_delim("./beds/NeanSNPs_4col.bed", delim = " ", col_names = FALSE)
  
