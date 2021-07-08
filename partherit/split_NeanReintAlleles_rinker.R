# 06.05.2021
# Gokberk Alagoz

# Script to refine neanderthal RA SNPs annotation
# from Rinker et al. (2020)

neanSNPs_rinker=read.csv2("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/41559_2020_1261_MOESM3_ESM.csv",stringsAsFactors=FALSE,header = T)
NDA=neanSNPs_rinker[neanSNPs_rinker$CLASS=="NDA",]
RA=neanSNPs_rinker[neanSNPs_rinker$CLASS!="NDA",]
AMH_derived_fixedVar=read.csv2("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/elife-63713-supp1-v1.csv",stringsAsFactors=FALSE,header = T)
nrow(NDA) # allele count matches the number in the paper - noice
nrow(RA) # allele count matches the number in the paper - noice
#write.csv2(NDA,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/41559_2020_1261_MOESM3_ESM_NDA-only.csv",row.names=F)
#write.csv2(RA,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/41559_2020_1261_MOESM3_ESM_RA-only.csv",row.names=F)

# make bed files, BED files should be formatted as "chr  pos-1  pos".
startPos = neanSNPs_rinker[2]-1

neanRA_bed = data.frame(chr = neanSNPs_rinker$chr,
                        startPos = startPos,
                        endPos = neanSNPs_rinker$POS)

write.table(neanRA_bed,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/nean_RA.bed",row.names=F,col.names=F,sep ="\t",quote=F)
#write.table(NDA[,c(1,startPos,2)],file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/NDA.txt",row.names=F,col.names=F,sep ="\t",quote=F)
#write.table(RA[,c(1,2,2)],file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/RA.txt",row.names=F,col.names=F,sep ="\t",quote=F)
#write.table(AMH_derived_fixedVar[,c(4,5,5)],file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/AMH_derived_fixedVariants_hg19-weiss_et_al.sorted.bed",row.names=F,col.names=F,sep ="\t",quote=F) # go manually remove headers from the output
#write.table(AMH_derived_fixedVar[,c(4,6,6)],file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh38/AMH_derived_fixedVariants_hg38-weiss_et_al.sorted.bed",row.names=F,col.names=F,sep ="\t",quote=F) # go manually remove headers from the output

# Reformat .csv files as chr:pos to run
# SNPsnap. SNPsnap will make a matching 
# SNP list (in terms of MAF and LD) to 
# run part. h2.
NDA_snplist=paste0(sub('chr','',NDA$chr),":",NDA$POS)
RA_snplist=paste0(sub('chr','',RA$chr),":",RA$POS)
#write.table(NDA_snplist,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/NDA_snplist_forSNPsnap.txt",row.names=F,col.names=F,quote=F)
#write.table(RA_snplist,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/supplementary_files/RA_snplist_forSNPsnap.txt",row.names=F,col.names=F,quote=F)

weiss_bed=read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/human-spec_noncodingVariants-weiss_et_al.bed",stringsAsFactors=F,header=F)
#write.table(weiss_bed,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/human-spec_noncodingVariants-weiss_et_al.bed",col.names = F,row.names = F,quote = F,sep="\t")

oneKGP3_samples=read.table("/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/integrated_call_samples_v3.20130502.ALL.panel",stringsAsFactors=F)
nrow(oneKGP3_samples)
oneKGP3_samples_EUR=oneKGP3_samples[oneKGP3_samples$V3=="EUR",]
nrow(oneKGP3_samples_EUR)
oneKGP3_samples_EURnonFIN=oneKGP3_samples_EUR[oneKGP3_samples_EUR$V2!="FIN",]
nrow(oneKGP3_samples_EURnonFIN)
#write.table(oneKGP3_samples_EURnonFIN$V1,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/1KGP3_EURnonFIN_samples.txt",col.names = F,row.names = F,quote = F,sep="\t")

# Extract tag Nean. SNPs with the highest
# MAF from each haplotype.

nean_tag_snps=read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/all_tag_snps.EUR.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed")
selected_tag_snps=data.frame()

for (haplotype in (levels(factor(nean_tag_snps$V16)))) {
  max_index=which.max(nean_tag_snps[nean_tag_snps$V16==haplotype,]$V7)
  tmp_tag_snp=data.frame(nean_tag_snps[nean_tag_snps$V16==haplotype,][max_index,])
  selected_tag_snps=rbind(selected_tag_snps,tmp_tag_snp)
}
selected_tag_snps_sorted=selected_tag_snps[order(selected_tag_snps$V1,selected_tag_snps$V2),]
list_for_snpsnap=paste0(substring(selected_tag_snps_sorted$V1,4),":",selected_tag_snps_sorted$V2)
#write.table(list_for_snpsnap,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/nean_tag_snps.txt",col.names = F,row.names = F,quote = F,sep = "\t")

# Extract Nean. RAs and NDAs with the highest
# AF from each haplotype.

selected_NDA=data.frame()
selected_RA=data.frame()

for (haplotype in (levels(factor(NDA$V16_LD_haplotype)))) {
  max_index=which.max(NDA[NDA$V16_LD_haplotype==haplotype,]$INTROG_AF)
  tmp_NDA=data.frame(NDA[NDA$V16_LD_haplotype==haplotype,][max_index,])
  selected_NDA=rbind(selected_NDA,tmp_NDA)
}
for (haplotype in (levels(factor(RA$V16_LD_haplotype)))) {
  max_index=which.max(RA[RA$V16_LD_haplotype==haplotype,]$INTROG_AF)
  tmp_RA=data.frame(RA[RA$V16_LD_haplotype==haplotype,][max_index,])
  selected_RA=rbind(selected_RA,tmp_RA)
}

hist(selected_RA_sorted$POS)

selected_NDA_sorted=selected_NDA[order(selected_NDA$chr,selected_NDA$POS),]
selected_RA_sorted=selected_RA[order(selected_RA$chr,selected_RA$POS),]
snpsnapList_NDA=paste0(substring(selected_NDA_sorted$chr,4),":",selected_NDA_sorted$POS)
snpsnapList_RA=paste0(substring(selected_RA_sorted$chr,4),":",selected_RA_sorted$POS)
length(snpsnapList_NDA)
length(snpsnapList_RA)
#write.table(snpsnapList_NDA,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/nean_snpsnapList_NDA.txt",col.names = F,row.names = F,quote = F,sep = "\t")
#write.table(snpsnapList_RA,"/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/nean_snpsnapList_RA.txt",col.names = F,row.names = F,quote = F,sep = "\t")

# Check allele freq. distribution

hist(NDA$INTROG_AF)
hist(RA$INTROG_AF)
hist(neanSNPs_rinker$INTROG_AF)
NDA_snpsnap=read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/snpsnap/NDA_snplist.txt",stringsAsFactors = F,sep = "\t")
RA_snpsnap=read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/snpsnap/RA_snplist.txt",stringsAsFactors = F,sep = "\t")
nrow(NDA_snpsnap)
nrow(RA_snpsnap)

#pdf("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/new_annotations/histogram_NDA.pdf")
histNDA=hist(NDA_snpsnap$V4)
#dev.off()
#pdf("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/new_annotations/histogram_RA.pdf")
histRA=hist(RA_snpsnap$V4)
#dev.off()