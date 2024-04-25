# Load packages
library("dplyr")
library('tidyr')
library('magrittr')
library('stringr')
library("reshape2")
library("DESeq2")
library("variancePartition")
library(forcats)
library(cowplot)
library(pheatmap)
library(ggplot2)
library(viridis)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library("genefilter")
library(vegan)
library(RColorBrewer)
library("plyr")
library("sva")
library("tibble")

# Load in Liometoxenus annotation files
anno_liom <- read.table("./08_DEGs_analysis/Liometoxenus/Liom_emapper_annotations.txt",
                        header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

anno_dc <- read.table("./08_DEGs_analysis/Dalotia/SMARTseq_Gland_analysis/Dcor_v3_annotation_only.txt",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

dcor_liom<-read.table("./08_DEGs_analysis/Liometoxenus/Dalotia_update.aa__v__Liometoxenus.aa.tsv",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

# convert Dcor to Liom table for one gene per row
sample_tibble <- dcor_liom %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(Liom_Gene = strsplit(as.character(Liometoxenus.aa), ", ")) %>%
  mutate(Dcor_Gene = strsplit(as.character(Dalotia_update.aa), ", ")) %>%
  unnest(Liom_Gene) %>%
  unnest(Dcor_Gene) %>%
  dplyr::ungroup()

head(sample_tibble)

# join tables
anno_liom<-anno_liom[,c(1,2,9,11)] %>%
  full_join(sample_tibble %>% dplyr::select(Liom_Gene,Dcor_Gene,Orthogroup), by=c("query"="Liom_Gene")) %>%
  left_join(anno_dc, by=c("Dcor_Gene"="Protein_ID"))

#read in counts table, header 
li<- read.delim("./08_DEGs_analysis/Liometoxenus/Liom_smartseq_gene_counts.txt",
                head=T, row.names=1,check.names = FALSE) 

# remove two outliers
li2<-li[,c(-4,-14)]
head(li2)

# coverage filter, 10 reads for at least 5 samples
keep <- rowSums(li2 >= 10) >= 5
li2 <- li2[keep,]

# load in Liometoxenus metadata
lioMeta <- read.table("./08_DEGs_analysis/Liometoxenus/sample_metadata_Liometoxenus.txt", 
                      header = TRUE, stringsAsFactors=FALSE, sep="\t")

# dataframe with all samples except two outliers
lioMeta <-lioMeta  %>% filter(sample !="22510") %>%filter(sample !="22512")
lioMeta 

lioMeta$condition<-factor(lioMeta$condition, level=c("control","D1","D2"))
lioMeta$batch <-factor(lioMeta$batch)

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsLN<-DESeqDataSetFromMatrix(countData=li2,
                              colData= lioMeta,
                              design = ~ batch + condition)

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsLN<-DESeq(ddsLN, betaPrior=FALSE)

# size factor added
colData(ddsLN)

#information on model tested by DESeq, should be standard
attr(ddsLN, "modelMatrixType")
attr(ddsLN, "modelMatrix")

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsLN)

# D1 cells v control cells
LN_resD1vC<-results(ddsLN, name="condition_D1_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(LN_resD1vC,alpha=0.05)
LN_resD1vC_tab<-as.data.frame(LN_resD1vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_liom %>% select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
LN_resD1vC_tabOrder<-LN_resD1vC_tab[order(LN_resD1vC_tab$padj),]
head(LN_resD1vC_tabOrder)
LN_resSigD1vC<-subset(LN_resD1vC_tab, padj <0.05) 
#write.csv(LN_resD1vC_tab, file= "Liom_D1vC_tab_ALL.csv")

# D2 cells v control cells
LN_resD2vC<-results(ddsLN, name="condition_D2_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(LN_resD2vC,alpha=0.05)
LN_resD2vC_tab<-as.data.frame(LN_resD2vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_liom %>% select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
LN_resD2vC_tabOrder<-LN_resD2vC_tab[order(LN_resD2vC_tab$padj),]
head(LN_resD2vC_tabOrder)
LN_resSigD2vC<-subset(LN_resD2vC_tab, padj <0.05) 
#write.csv(LN_resD2vC_tab, file= "Liom_D2vC_tab_ALL.csv")

# D2 cells v D1 cells
LN_resD1vD2<-results(ddsLN, contrast=list("condition_D2_vs_control","condition_D1_vs_control"),alpha=0.05,cooksCutoff=FALSE)
summary(LN_resD1vD2, alpha=0.05)
LN_resD1vD2_tab<-as.data.frame(LN_resD1vD2)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_liom %>% select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
LN_resD1vD2_tabOrder<-LN_resD1vD2_tab[order(LN_resD1vD2_tab$padj),]
head(LN_resD1vD2_tabOrder)
LN_resSigD1vD2<-subset(LN_resD1vD2_tab, padj <0.05) 
#write.csv(LN_resD1vD2_tab, file= "Liom_D1vD2_tab_ALL.csv")

###################################################
# variancePartition analysis
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(ddsLN)>1) >= 0.5 * ncol(ddsLN)

#compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( ddsLN )[isexpr,] + 1)

# Define formula (batch = machine= reads)
form <- ~  (1|batch) + (1|condition)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, lioMeta)

# plot results
vp <- sortCols(varPart)
a<-plotVarPart(vp, col=c("grey90","black","white"))

# get median % variation by random effect
summary(vp)

# variability among replicates by cell type, batch, etc.
# extract normalized read counts
cd <- counts(ddsLN,normalized=TRUE)

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(lioMeta$condition[match(.$ind,lioMeta$sample)]),
         batch=as.factor(lioMeta$batch[match(.$ind,lioMeta$sample)]))

# median, mean and variance calculated by cell type
# calculate by individual and then cell type
med.dat<- dat %>% 
  group_by(ind) %>% 
  dplyr::mutate(varcount =var(values, na.rm = TRUE),median = median(log2(values)+1, na.rm = TRUE),
                avg=mean(values, na.rm = TRUE)) %>%
  mutate(varcount=log2(varcount)+1, avg=log2(avg)+1) %>%
  group_by(condition) %>% 
  dplyr::select(ind,varcount,median,avg) %>%
  distinct()

ggplot(med.dat, aes(x=condition,y=varcount))+
  geom_boxplot()+geom_point() +ylim(0,35)

# plot normalized counts by cell type
b<-ggplot(dat,aes(x = fct_reorder(ind,log2(values)+1,.fun='median'), y = log2(values)+1, fill=batch)) + 
  geom_boxplot(outlier.size = 0.1,lwd=0.1)+ xlab("Samples") + ylab(expression("Log"["2"]~"Normalized Counts"))+
  scale_fill_grey(start = 1,end = 0,)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")+
  facet_grid(~ condition, scales="free")

# principal component analysis on VST transformed counts

# variance stabilized transformed counts
vst<-varianceStabilizingTransformation(ddsLN, blind=FALSE)
vstMat<-assay(vst)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

# plotPCA from DESeq2 uses top 500 genes by default
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))

# PC variance explained
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsLN)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d_500 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, names = colnames(ddsLN))

c<- ggplot(data = d_500, aes_string(x = "PC1", y = "PC2", color = "condition", shape="batch")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# repeat but will all genes
pcaALL <- prcomp(t(assay(vst)))

# PC variance explained
percentVarALL <- pcaALL$sdev^2/sum(pcaALL$sdev^2)

intgroup.df <- as.data.frame(colData(ddsLN)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
dALL <- data.frame(PC1 = pcaALL$x[, 1], PC2 = pcaALL$x[, 2], group = group, 
                   intgroup.df, names = colnames(ddsLN))

d<- ggplot(data = dALL, aes_string(x = "PC1", y = "PC2", color = "condition", shape="batch")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVarALL[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVarALL[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# put all the pieces together
plot_grid(a,b,c,d,nrow = 1, ncol = 4, rel_widths = c(0.3,0.4,0.25,0.25))

###########
## Volcano plots

vhigh<-read.table("./08_DEGs_analysis/Dcor_Alus_Liom_highlightVolcano_genes.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t")

# add column to indicate significant result
LN_res <- dplyr::mutate(LN_resD1vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
LN_res2<-LN_resD1vC_tab[LN_resD1vC_tab$Protein_ID %in% vhigh$Liom_gene_id,] %>% 
  left_join(vhigh %>% select(Liom_gene_id,Gene_group_liom, Cell_type_ltx), by=c("Protein_ID"="Liom_gene_id")) %>%
  distinct(Protein_ID,.keep_all = TRUE)

# plot all things
p <- ggplot(LN_res, aes(log2FoldChange, -log10(padj))) +theme_classic()
p <- p + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
p <- p + geom_point(data=LN_res2,aes(log2FoldChange, -log10(padj),color=Cell_type_ltx), alpha = 0.5, size=1.5) 
p <- p + geom_text_repel(data=LN_res2,aes(log2FoldChange, -log10(padj), label=Gene_group_liom,color=Cell_type_ltx), fontface=4,size = 4) 
p <- p + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#00D500",S="#B326E2",G="#FEF02E",terpene="#1811FD"))
p <- p + xlab('BQ cells/Control cells (Log2)') 
p <- p + ylab('-log10(adjusted p-value)') 
p <- p + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p

# add column to indicate significant result
LN_res3 <- dplyr::mutate(LN_resD2vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
LN_res4<-LN_resD2vC_tab[LN_resD2vC_tab$Protein_ID %in% vhigh$Liom_gene_id,] %>% 
  left_join(vhigh %>% select(Liom_gene_id,Gene_group_liom, Cell_type_ltx), by=c("Protein_ID"="Liom_gene_id")) %>%
  distinct(Protein_ID,.keep_all = TRUE)

# plot all things D2 v control
q <- ggplot(LN_res3, aes(log2FoldChange, -log10(padj))) +theme_classic()
q <- q + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
q <- q + geom_point(data=LN_res4,aes(log2FoldChange, -log10(padj),color=Cell_type_ltx), alpha = 0.5, size=1.5) 
q <- q + geom_text_repel(data=LN_res4,aes(log2FoldChange, -log10(padj), label=Gene_group_liom,color=Cell_type_ltx), fontface=4,size = 4) 
q <- q + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#00D500",S="#B326E2",G="#FEF02E",terpene="#1811FD"))
q <- q + xlab('Solvent cells/Control cells (Log2)')+xlim(-11,11) 
q <- q + ylab('-log10(adjusted p-value)') 
q <- q + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

q

#grid.arrange(p,q,nrow=1)
#ggarrange(p, q + rremove("ylab"), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

# D1 v D2
# add column to indicate significant result
LN_res5 <- dplyr::mutate(LN_resD1vD2_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
LN_res6<-LN_resD1vD2_tab[LN_resD1vD2_tab$Protein_ID %in% vhigh$Liom_gene_id,] %>% 
  left_join(vhigh %>% select(Liom_gene_id,Gene_group_liom, Cell_type_ltx), by=c("Protein_ID"="Liom_gene_id")) %>%
  distinct(Protein_ID,.keep_all = TRUE) %>% filter(Cell_type_ltx !="G")

r <- ggplot(LN_res5, aes(log2FoldChange, -log10(padj))) +theme_classic()
r <- r + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
r <- r + geom_point(data=LN_res6 ,aes(log2FoldChange, -log10(padj),color=Cell_type_ltx), alpha = 0.5, size=1.5) 
r <- r + geom_point(data=LN_res6 %>% filter(!is.na(Gene_group_liom)),aes(log2FoldChange, -log10(padj),color=Cell_type_ltx), alpha = 1, size=3) 
r <- r + geom_text_repel(data=LN_res6,aes(log2FoldChange, -log10(padj), label=Gene_group_liom,color=Cell_type_ltx), fontface=4,size = 4) 
r <- r + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#00D500",S="#B326E2",terpene="#1811FD"))
r <- r + xlab('Solvent cells/BQ cells (Log2)')+xlim(-15,15)
r <- r + ylab('-log10(adjusted p-value)') 
r <- r + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
r

###################################################
# transform counts
vstLN<-varianceStabilizingTransformation(ddsLN, blind=FALSE)
vstMatLN<-assay(vstLN)

#principle component analysis
# plotPCA uses top 500 genes by default
rv <- rowVars(assay(vstLN))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vstLN)[select, ]))

# all genes
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsLN)[, c("condition"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = colnames(ddsLN))

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label="names")) + 
  geom_point(size = 4,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#00D500",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  geom_text_repel()+
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=14))+theme(legend.title=element_blank()) 

###############################################
# calculate fold change of the TPMs by cell type
FC_LN<-as.data.frame(vstMatLN) %>% 
  mutate(avg_ctrl=round((`22503`+`22505`+`22507`+`22513`+`25`+`26`)/6,3),
         avg_D1=round((`27`+`29`+`30`+`31`+`32`+`22508`)/6,3),
         avg_D2=round((`22504`+`22506`+`22509`+`22511`+`22514`)/5,3),
         ctrl_d1_ratio=round(avg_D1-avg_ctrl,3), ctrl_d2_ratio=round(avg_D2-avg_ctrl,3),
         d1_d2_ratio=round(avg_D1-avg_D2,3),
         Protein_ID=rownames(.)) %>%
  full_join(li %>% mutate(Protein_ID=rownames(.)) %>% select(Protein_ID), by="Protein_ID")%>%
  full_join(anno_liom %>% select(query,Dcor_Gene, Annotation) ,by=c("Protein_ID"="query")) %>%
  left_join(LN_resD1vC_tab %>% select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D1vC = log2FoldChange, padj_D1vC =padj) %>%
  left_join(LN_resD2vC_tab %>% select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vC = log2FoldChange, padj_D2vC =padj) %>%
  left_join(LN_resD1vD2_tab %>% select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vD1 = log2FoldChange, padj_D2vD1 =padj) %>%
  distinct(Protein_ID,.keep_all = TRUE) %>%
  mutate(change=ifelse(beta_D1vC > 2 & d1_d2_ratio > 2 & padj_D1vC < 0.05 , "D1-enriched",
                       ifelse(beta_D2vC > 2 & d1_d2_ratio < -2 & padj_D2vC < 0.05, "D2-enriched",
                              ifelse(beta_D1vC > 2 & beta_D2vC > 2 & d1_d2_ratio > -2 & d1_d2_ratio < 2 & padj_D1vC < 0.05 & padj_D2vC < 0.05,"gland-enriched",
                                     ifelse(beta_D1vC < -2 & beta_D2vC < -2 & padj_D1vC < 0.05 & padj_D2vC < 0.05, "ctrl-enriched",
                                            ifelse(padj_D1vC < 0.05 | padj_D2vC < 0.05 | padj_D2vD1 < 0.05, "other","not_sig")))))) %>%
  mutate(change=ifelse(is.na(padj_D1vC)|is.na(padj_D2vC)|is.na(padj_D2vD1), "not_tested", change))

#write.csv(FC_LN,"FC_enrichment_Liometoxenus.csv")

# tabulate number of genes per category above
summary_FC_LN <- FC_LN %>% select(change, Protein_ID) %>%
  group_by(change) %>% dplyr::summarize(n = n())
summary_FC_LN

####################
# plot heat map of significant genes
# create dataframe for heat map
sub_DEGs<- FC_LN %>% filter(change == "D1-enriched") %>%
  arrange(desc(beta_D1vC)) %>%
  select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
         -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1,-Dcor_Gene) %>%
  distinct(.[,1:17],.keep_all = TRUE)

sub_DEGs<- FC_LN %>% filter(change == "D2-enriched") %>%
  arrange(desc(beta_D2vC)) %>%
  select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
         -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1,-Dcor_Gene) %>%
  distinct(.[,1:17],.keep_all = TRUE)

rownames(sub_DEGs)<-sub_DEGs$Protein_ID

dim(sub_DEGs)

# set column annotation
col.anno<- data.frame(row.names=colnames(sub_DEGs[1:75,c(-18,-19)]),
                      condition=as.factor(lioMeta$condition[match(colnames(sub_DEGs[,c(-18,-19)]),lioMeta$sample)]),
                      segment=as.factor(lioMeta$segment[match(colnames(sub_DEGs[,c(-18,-19)]),lioMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#00D500",D2="#B326E2"),
               segment=c(seg6="grey",seg7="black"))


# diverging palette is good for scaled TPMs
mypalette <- magma(100)

mat<-as.matrix(sub_DEGs[1:75,c(-18,-19)])
rownames(mat) <- paste(rownames(mat),sub_DEGs[1:75,19])
rownames(mat) <- paste(sub_DEGs[1:75,19])


p <- pheatmap(mat,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")

##################################################

###################
# joined metadata #
###################

#create new metadata for the joint data
all_samples<-rbind.fill(alMeta,dcMeta,lioMeta)
all_samples$species<-c(rep("Alus",14),rep("Dcor",32), rep("Ltx",17))
rownames(all_samples)<-all_samples$sample

# extract orthologs with shared enriched cell-type expression for Aleochara and Dalotia
og <-read.table("./08_DEGs_analysis/sub_OGs_TPM_update.txt", sep="\t", header=T)

##############
# shared OGs #
##############

# extract orthologs with shared enriched cell-type expression for Aleochara and Dalotia
OGs_dc_alus<-sample_tibble %>%
  select(Orthogroup, group,Dcor_Gene,Liom_Gene) %>%
  left_join(sample_tibble2 %>% select(Dcor_Gene, Alus_Gene), by=c("Dcor_Gene")) %>%
  left_join(dc3 %>% rownames_to_column(var="Dcor_Gene"), by=c("Dcor_Gene")) %>%
  left_join(li2 %>% rownames_to_column(var="Liom_Gene"), by=c("Liom_Gene")) %>%
  left_join(al3 %>% rownames_to_column(var="Alus_Gene"), by=c("Alus_Gene")) 

fi <- filter(OGs_dc_alus, Dcor_Gene %in% og$Dcor_Gene & Alus_Gene %in% og$Alus_Gene & !is.na(`25` ))

count_matrix <- as.matrix(fi[,6:68])

all_samples <- all_samples[match(colnames(count_matrix),all_samples$sample),]

batch <- c(rep(1, 32), rep(2, 14))

# remove batch effect

adjusted <- ComBat_seq(count_matrix , batch=all_samples$species, group=NULL)

#####################
# three species PCA #
#####################

all_samples$condition<-factor(all_samples$condition, level=c("control","D1","D2"))

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
  geom_point(size = 4,alpha = 0.8) + scale_shape_manual(values=c(17,19,18)) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=14))+theme(legend.title=element_blank()) 

ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "condition",  shape="species")) + 
  geom_point(size = 4,alpha = 0.8) + scale_shape_manual(values=c(17,19,18)) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC3: ", round(percentVar[3] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=14))+theme(legend.title=element_blank())


#####################

# extract orthologs with shared enriched cell-type expression for Aleochara and Dalotia
HOGs_dc_al_liom<-sample_tibble %>%
  select(Orthogroup, group,Dcor_Gene,Liom_Gene) %>%
  left_join(FC_LN %>% select(Protein_ID,change,beta_D1vC,beta_D2vC,padj_D1vC,padj_D2vC,Annotation), by=c("Liom_Gene"="Protein_ID")) %>%
  set_colnames(c("OG","group","Dcor_Gene","Liom_Gene","Liom_change","Liom_beta_D1","Liom_beta_d2","Liom_padj_D1vC","Liom_padj_D2vC","Annotation")) %>%
  left_join(FC_dc %>% select(Protein_ID,change,beta_D1vC,beta_D2vC,padj_D1vC,padj_D2vC), by=c("Dcor_Gene"="Protein_ID")) %>%
  set_colnames(c("OG","group","Dcor_Gene","Liom_Gene","Liom_change","Liom_beta_D1","Liom_beta_d2","Liom_padj_D1vC","Liom_padj_D2vC","Annotation",
                 "Dcor_change","Dcor_beta_D1","Dcor_beta_d2","Dcor_padj_D1vC","Dcor_padj_D2vC")) %>%
  left_join(FC_al %>% select(Protein_ID,change,beta_D1vC,beta_D2vC,padj_D1vC,padj_D2vC, Dcor_Gene), by=c("Dcor_Gene")) %>%
  set_colnames(c("OG","group","Dcor_Gene","Liom_Gene","Liom_change","Liom_beta_D1","Liom_beta_d2","Liom_padj_D1vC","Liom_padj_D2vC","Annotation",
                 "Dcor_change","Dcor_beta_D1","Dcor_beta_d2","Dcor_padj_D1vC","Dcor_padj_D2vC",
                 "Alus_Gene","Alus_change","Alus_beta_D1","Alus_beta_d2","Alus_padj_D1vC","Alus_padj_D2vC")) %>%
  mutate(Annotation=ifelse(is.na(Annotation),anno_dc$Annotation[match(.$Dcor_Gene, anno_dc$Protein_ID)], Annotation)) %>%
  mutate(Annotation=ifelse(is.na(Annotation),anno$Annotation[match(.$Liom_Gene, anno$query)], Annotation))

head(HOGs_dc_al_liom)

#write.table(HOGs_dc_liom, "./08_DEGs_analysis/Liometoxenus/Dcor_Liom_orthologs_enriched_expression.txt", sep="\t", row.names=FALSE)

shared_D1<-HOGs_dc_al_liom %>% filter(Dcor_change == "D1-enriched" & Liom_change == "D1-enriched" ) 
shared_D2<-HOGs_dc_al_liom %>% filter(Dcor_change == "D2-enriched" & Liom_change == "D2-enriched" ) 

#write.table(shared_D1, "./08_DEGs_analysis/Liometoxenus/Dcor_Liom_orthologs_D1-enriched_expression.txt", sep="\t", row.names=FALSE)
#write.table(shared_D2, "./08_DEGs_analysis/Liometoxenus/Dcor_Liom_orthologs_D2-enriched_expression.txt", sep="\t", row.names=FALSE)

