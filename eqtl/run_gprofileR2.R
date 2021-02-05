library(gprofiler2)
options(stringsAsFactors = FALSE);
setwd('~/Google Drive/Papers/ENIGMA-Evol/eQTLgenes/');

##Load in global gene list
GlobalSAHGE7PCW = read.csv('GlobalSAHGE7PCW.csv');
GO <- gost(GlobalSAHGE7PCW$ensembl_gene_id,
    organism = "hsapiens",
    exclude_iea = TRUE,
    correction_method = "fdr",
    significant = TRUE,
    sources = c("GO", "KEGG", "REAC"),
    ordered_query = FALSE,
    numeric_ns = "ENTREZGENE_ACC",
    user_threshold = 0.05,
    domain_scope = "annotated"
  )
GO$result
##Make a bar plot


##Load in regional gene list
RegionalSAHGE7PCW = read.csv('RegionalSAHGE7PCW.csv');
GO <- gost(RegionalSAHGE7PCW$ensembl_gene_id,
           organism = "hsapiens",
           exclude_iea = TRUE,
           correction_method = "fdr",
           significant = TRUE,
           sources = c("GO", "KEGG", "REAC"),
           ordered_query = FALSE,
           numeric_ns = "ENTREZGENE_ACC",
           user_threshold = 0.05,
           domain_scope = "annotated"
)




