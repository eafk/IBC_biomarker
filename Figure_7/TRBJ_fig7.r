rm(list=ls())
gc(reset = T,full = T)
library(tidyverse)
library(RColorBrewer)
library(msigdbr)

#### heatmap of NMD in PBMC
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)

## get GO:BP genesets
gene_sets = msigdbr(species = "human", collection =  "C5",subcollection = "GO:BP")
gene_sets = split(x = as.character(gene_sets$gene_symbol), f = gene_sets$gs_name)

### extract geneset related to NMD
names(gene_sets)[grep(toupper("nonsense_mediated_decay"),names(gene_sets))]
## get the gene id
names_list<-sort(c(gene_sets[["GOBP_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS_NONSENSE_MEDIATED_DECAY"]],"HNRNPL"))

### add TRBJs and other interested genes
names_list<-c("TRBJ1-1","TRBJ1-2","TRBJ1-6","ANGEL1","ANGEL2","CNP","ANG","ERN1","ENDOU","USB1","N4BP2",names_list)

## read log2 transformed normalized (by DESeq2) counts using genome mapping results
## from DOI: 10.1126/sciadv.adu0031, Sci Adv IBC paper

## TRBJ data were generated after 30 nt filtering step
dat<-read.delim("genome.data")

dat<-dat %>%
  arrange(match(ID, names_list))

pdf("PBMC_NMD_genome.pdf",height=11,width=8)
heatmap.2(as.matrix(dat[,-1]),labRow = dat$ID,dendrogram = "none",
          scale="none",margins = c(10, 20),breaks=seq(-3,3,0.06),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=F,cexCol = 0.5,
          col =heatPalette)
dev.off()

## read log2 transformed normalized (by DESeq2) counts using transcriptome mapping results
## from DOI: 10.1126/sciadv.adu0031, Sci Adv IBC paper

## TRBJ data were generated after 30 nt filtering step, the same as what was used in the above plot

dat<-read.delim("transcriptome.data")

dat<-dat %>%
  arrange(match(ID, names_list))

pdf("PBMC_NMD_transcriptome.pdf",height=11,width=8)
heatmap.2(as.matrix(dat[,-1]),labRow = dat$ID,dendrogram = "none",
          scale="none",margins = c(10, 20),breaks=seq(-3,3,0.06),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=F,cexCol = 0.5,
          col =heatPalette)
dev.off()