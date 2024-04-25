library("dplyr")
library('tidyr')
library('magrittr')
library('stringr')
library("reshape2")
library("DESeq2")
library("plyr")
library("sva")
library(forcats)
library("genefilter")
library(vegan)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(tibble)

###########
# Dalotia #
###########

#read in counts table, header 
dc<- read.delim("E:/Caltech/Beetle Genomes/Dalotia_github/RNAseq_files/Dcor_smartseq_gene_counts.txt",
                head=T, row.names=1)
dc3<-dc[,c(2:4,6:34)]

# read in Dalotia smartseq tpm table and filter for enriched genes
dcMeta <- read.table("E:/Caltech/Beetle Genomes/Dalotia_github/RNAseq_files/smartseq_samples_gland_control.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(dcMeta)

dcMeta$condition<-factor(dcMeta$condition, level=c("control","D1","D2"))
dcMeta$batch <-factor(dcMeta$batch)

# gene length
dc_Len<-data.frame(target_id=rownames(dc),Length=dc$Length)

# orthologues
dcor_alus<-read.table("E:/Caltech/Beetle Genomes/Assemblies/56_Alustrica/smartseq/Dalotia_update.aa__v__Alustrica.aa.tsv",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

# OrthoFinder Results Dcor_Alus
# split to 1-1 relationships
sample_tibble2 <- dcor_alus %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(Alus_Gene = strsplit(as.character(Alustrica.aa), ", ")) %>%
  mutate(Dcor_Gene = strsplit(as.character(Dalotia_update.aa), ", ")) %>%
  unnest(Alus_Gene) %>%
  unnest(Dcor_Gene) %>%
  ungroup()

# variance stabilized counts
ddsMF_all<-DESeqDataSetFromMatrix(countData=dc3,
                                  colData= dcMeta,
                                  design = ~ batch + condition)
vst_all<-varianceStabilizingTransformation(ddsMF_all, blind=FALSE)
vstMat_all<-assay(vst_all)

#################
# Aleochara sp3 #
#################

#read in counts table, header 
al<- read.delim("E:/Caltech/Beetle Genomes/Assemblies/56_Alustrica/smartseq/Alus_smartseq_gene_counts.txt",
                head=T, row.names=1,check.names = FALSE)
al3<-al[,c(-1,-5,-11)]

# read in Aleochara smartseq tpm table and filter for enriched genes
alMeta <- read.table("E:/Caltech/Beetle Genomes/Assemblies/56_Alustrica/sample_metadata.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(alMeta)

alMeta$condition<-factor(alMeta$condition, level=c("control","D1","D2"))
alMeta$batch <- 5

# gene length
al_Len<-data.frame(target_id=rownames(al),Length=al$Length)

# variance stabilized counts
ddsAL_all<-DESeqDataSetFromMatrix(countData=al3,
                                  colData= alMeta,
                                  design = ~ condition)

vstAL_all<-varianceStabilizingTransformation(ddsAL_all, blind=FALSE)
vstMatAL_all<-assay(vstAL_all)

###################
# joined metadata #
###################

#create new metadata for the joint data
all_samples<-rbind.fill(alMeta,dcMeta)
all_samples$species<-c(rep("Alus",14),rep("Dcor",32))
rownames(all_samples)<-all_samples$sample

all_samples$condition<-factor(all_samples$condition, level=c("D1","D2","control"))

#######################
# identify shared OGs #
#######################
# merge tables
OGs_dc_alus<-sample_tibble2 %>%
  dplyr::select(Orthogroup, group,Dcor_Gene,Alus_Gene) %>%
  left_join(FC_al %>% dplyr::select(Protein_ID,change,beta_D1vC,beta_D2vC,padj_D1vC,padj_D2vC,Annotation), by=c("Alus_Gene"="Protein_ID")) %>%
  set_colnames(c("OG","group","Dcor_Gene","Alus_Gene","Alus_change","Alus_beta_D1","Alus_beta_d2","Alus_padj_D1vC","Alus_padj_D2vC","Annotation")) %>%
  left_join(FC_dc %>% dplyr::select(Protein_ID,change,beta_D1vC,beta_D2vC,padj_D1vC,padj_D2vC), by=c("Dcor_Gene"="Protein_ID")) %>%
  set_colnames(c("OG","group","Dcor_Gene","Alus_Gene","Alus_change","Alus_beta_D1","Alus_beta_d2","Alus_padj_D1vC","Alus_padj_D2vC","Annotation",
                 "Dcor_change","Dcor_beta_D1","Dcor_beta_d2","Dcor_padj_D1vC","Dcor_padj_D2vC")) %>%
  mutate(Annotation=ifelse(is.na(Annotation),anno_dc$Annotation[match(.$Dcor_Gene, anno_dc$Protein_ID)], Annotation)) %>%
  mutate(Annotation=ifelse(is.na(Annotation),anno$Annotation[match(.$Alus_Gene, anno$query)], Annotation))

head(OGs_dc_alus)

# find single ortholog pair between Aleochara and Dalotia, n=9314
sub_OGs_TPM<-OGs_dc_alus %>%
  left_join(al_Len %>% dplyr::select(target_id, Length) %>% distinct(), by=c("Alus_Gene"="target_id")) %>% 
  dplyr::rename(Alus_len = Length) %>%
  left_join(dc_Len %>% dplyr::select(target_id, Length) %>% distinct(), by=c("Dcor_Gene"="target_id")) %>%
  dplyr::rename(Dcor_len = Length) %>%
  left_join(data.frame(vstMatAL_all,check.names=F) %>% mutate(alus_id=rownames(.)), by=c("Alus_Gene"="alus_id")) %>%
  left_join(data.frame(vstMat_all,check.names=F) %>% mutate(dcor_id=rownames(.)), by=c("Dcor_Gene"="dcor_id")) %>%
  dplyr::arrange(ifelse(Dcor_padj_D1vC<Dcor_padj_D2vC,Dcor_padj_D1vC,Dcor_padj_D2vC),ifelse(Alus_padj_D1vC<Alus_padj_D2vC,Alus_padj_D1vC,Alus_padj_D2vC)) %>%
  group_by(group)%>%
  dplyr::slice(1) %>%
  ungroup()

#write.table(sub_OGs_TPM,"./08_DEGs_analysis/sub_OGs_TPM_update.txt", sep="\t", row.names=FALSE)

#######################
# import previous table
og <-read.table("./08_DEGs_analysis/sub_OGs_TPM_update.txt", sep="\t", header=T)

# extract orthologs with shared enriched cell-type expression for Aleochara and Dalotia
OGs_dc_alus<-sample_tibble2 %>%
  dplyr::select(Orthogroup, group,Dcor_Gene,Alus_Gene) %>%
  left_join(dc3 %>% rownames_to_column(var="Dcor_Gene"), by=c("Dcor_Gene")) %>%
  left_join(al3 %>% rownames_to_column(var="Alus_Gene"), by=c("Alus_Gene")) 


fi <- filter(OGs_dc_alus, Dcor_Gene %in% og$Dcor_Gene & Alus_Gene %in% og$Alus_Gene)

count_matrix <- as.matrix(fi[,5:50])

all_samples <- all_samples[match(colnames(count_matrix),all_samples$sample),]

# adjust count matrix for batch/species effect.
adjusted <- ComBat_seq(count_matrix , batch=all_samples$species, group=NULL)

#######
# PCA #
#######

# extract all shared orthologs (OGs) with shared expression
ddsMF_all<-DESeqDataSetFromMatrix(countData=adjusted,
                                  colData= all_samples,
                                  design = ~ condition)
vst_all<-varianceStabilizingTransformation(ddsMF_all, blind=FALSE)
vstMat_all<-assay(vst_all)

plotPCA(vst_all, intgroup = c("species", "condition"), ntop=9000)

# plot PCA based on 
pc<-prcomp(t(vstMat_all), center=TRUE, scale=FALSE)
summary(pc)

percentVar <- pc$sdev^2/sum(pc$sdev^2)

d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3= pc$x[, 3],PC4= pc$x[, 4],species = all_samples$species, 
                condition=all_samples$condition, names = all_samples$sample)

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "condition",  shape="species")) + 
  geom_point(size = 8,alpha = 0.8) + scale_shape_manual(values=c(17,19,18)) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=14))+theme(legend.title=element_blank()) 

ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "condition",  shape="species")) + 
  geom_point(size = 4,alpha = 0.8) + scale_shape_manual(values=c(17,19,18)) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  ylab(paste0("PC3: ", round(percentVar[3] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=14))+theme(legend.title=element_blank())

