#
# This script will read a GWAS hit region and will generate a dotPlot
# Gokberk Alagoz - 28.10.21

library(seqinr)

nucacd_human <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/human_rs73004625_fordotter.fa")[[1]]
nucacd_chimp <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/chimp_rs73004625_fordotter.fa")[[1]]

dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp)
dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp, wsize = 2, wstep = 2, nmatch = 2)

nucacd_human <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/human_rs73004625_fordotter_2.fa")[[1]]
nucacd_chimp <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/chimp_rs73004625_fordotter_2.fa")[[1]]

dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp)
dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp, wsize = 2, wstep = 2, nmatch = 2)

nucacd_chimp <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/chimp_rs73004625_fordotter_2.fa")[[1]]
nucacd_human <- read.fasta("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/human_rs73004625_fordotter_2.fa")[[1]]

dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp, wsize = 2, wstep = 2, nmatch = 2)
dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp, wsize = 2, wstep = 2, nmatch = 2)
dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp)
dotPlot(seq1 = nucacd_human, seq2 = nucacd_chimp, wsize = 2, wstep = 2, nmatch = 2)