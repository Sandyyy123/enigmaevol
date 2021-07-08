library(gprofiler2)
options(stringsAsFactors = FALSE)
#setwd('~/Google Drive/Papers/ENIGMA-Evol/eQTLgenes/')

##Load in global gene list
GlobalSAHGE = read.csv("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/GlobalSAfetalHGE.csv")
GO <- gost(GlobalSAHGE$ensembl_gene_id,
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

## Make a bar plot


##Load in regional gene list
RegionalSAHGE = read.csv("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/regionalSAfetalHGE.csv")
GO <- gost(RegionalSAHGE$ensembl_gene_id,
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

## Make a bar plot
