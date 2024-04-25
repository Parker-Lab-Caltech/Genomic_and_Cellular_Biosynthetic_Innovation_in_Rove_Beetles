library("plyr")
library("dplyr")
library('tidyr')
library('magrittr')
library('stringr')
library("reshape2")
library("DESeq2")
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpubr)

#######################
## Dalotia Bulk RNAseq gland v. control segments
# load in annotation for Dalotia
anno_dc <- read.table("./08_DEGs_analysis/Dalotia/SMARTseq_Gland_analysis/Dcor_v3_annotation_only.txt",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

#read in counts table, header 
dc<- read.delim("./08_DEGs_analysis/Dalotia/Bulk_Gland_analysis/Dcor_bulk_gland-control_gene_counts.txt",
                head=T, row.names=1) 

dc2<-dc[,c(2:13)]
head(dc2)

# coverage filter, 10 reads for at least 3 samples
keep <- rowSums(dc2 >= 10) >= 3
dc2 <- dc2[keep,]

# read in Dalotia smartseq tpm table and filter for enriched genes
dcMeta <- read.table("./08_DEGs_analysis/Dalotia/Bulk_Gland_analysis/bulk_samples_gland_control.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(dcMeta)

dcMeta$condition<-factor(dcMeta$condition, level=c("control","seg7"))
dcMeta$tecRep <-factor(dcMeta$tecRep)


#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsB<-DESeqDataSetFromMatrix(countData=dc2,
                             colData= dcMeta,
                             design = ~ condition)

#collapse replicates
ddsColl <- collapseReplicates(ddsB, ddsB$sample, ddsB$tecRep)

# examine the colData and column names of the collapsed data
colData(ddsColl)
colnames(ddsColl)

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsColl<-DESeq(ddsColl, betaPrior=FALSE)

# size factor added
colData(ddsColl)

#information on model tested by DESeq, should be standard
attr(ddsColl, "modelMatrixType")
attr(ddsColl, "modelMatrix")

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsColl)

# gland cells v control cells
res7vC<-results(ddsColl, name="condition_seg7_vs_control",alpha=0.05,cooksCutoff=TRUE)
summary(res7vC,alpha=0.05)
res7vC_tab<-as.data.frame(res7vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_dc,by="Protein_ID")
res7vC_tabOrder<-res7vC_tab[order(res7vC_tab$padj),]
head(res7vC_tabOrder)
resSig7vC<-subset(res7vC_tab, padj <0.05) 
#write.csv(res7vC_tab, file= "./08_DEGs_analysis/Dalotia/Bulk_Gland_analysis/Dalotia_BULK_seg7vC_tab_ALL.csv")

###################################################
# variancePartition analysis
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(ddsB)>1) >= 0.5 * ncol(ddsB)

#compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( ddsB )[isexpr,] + 1)

# Define formula
form <- ~  (1|tecRep) + (1|condition)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, dcMeta)

# plot results
vp <- sortCols(varPart)
a<-plotVarPart(vp, col=c("grey90","black","white"))

# get median % variation by random effect
summary(vp)

# variability among replicates by cell type, batch, etc.
# extract normalized read counts
ddsB<-DESeq(ddsB, betaPrior=FALSE)
cd <- counts(ddsB,normalized=TRUE)

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(dcMeta$condition[match(.$ind,dcMeta$id2)]),
         batch=as.factor(dcMeta$tecRep[match(.$ind,dcMeta$id2)]))

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

vst<-varianceStabilizingTransformation(ddsB, blind=FALSE)
vstMat<-assay(vst)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

# plotPCA from DESeq2 uses top 500 genes by default
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))

# PC variance explained
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsB)[, c("condition","tecRep"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d_500 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, names = colnames(ddsB))

c<- ggplot(data = d_500, aes_string(x = "PC1", y = "PC2", color = "condition", shape="tecRep", label="names")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",seg7="#0ECA7C"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  geom_text_repel(size = 3,box.padding = unit(0.5, "lines"))+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# repeat but will all genes
pcaALL <- prcomp(t(assay(vst)))

# PC variance explained
percentVarALL <- pcaALL$sdev^2/sum(pcaALL$sdev^2)

intgroup.df <- as.data.frame(colData(ddsB)[, c("condition","tecRep"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
dALL <- data.frame(PC1 = pcaALL$x[, 1], PC2 = pcaALL$x[, 2], group = group, 
                   intgroup.df, names = colnames(ddsB))

d<- ggplot(data = dALL, aes_string(x = "PC1", y = "PC2", color = "condition", shape="tecRep",label = "names")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",seg7="#0ECA7C"))+
  xlab(paste0("PC1: ", round(percentVarALL[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVarALL[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  geom_text_repel(size = 3,box.padding = unit(0.5, "lines"))+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# put all the pieces together
plot_grid(a,b,c,d,nrow = 1, ncol = 4, rel_widths = c(0.3,0.4,0.25,0.25))

#######################################
## Volcano plots
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggpubr)

vhigh<-read.table("./08_DEGs_analysis/Dcor_Alus_Liom_highlightVolcano_genes.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t")

# add column to indicate significant result
res <- dplyr::mutate(res7vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
res2<-res7vC_tab[res7vC_tab$Protein_ID %in% vhigh$Dcor_v3_gene_id,] %>% 
  left_join(vhigh %>% select(Dcor_v3_gene_id,Gene_group_dcor, Cell_type), by=c("Protein_ID"="Dcor_v3_gene_id"))

# plot all things
p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +theme_classic()
p <- p + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
p <- p + geom_point(data=res2%>% filter(padj < 0.05),aes(log2FoldChange, -log10(padj),color=Cell_type), alpha = 0.5, size=1.5) 
p <- p + geom_point(data=res2 %>% filter(!is.na(Gene_group_dcor),padj < 0.05),aes(log2FoldChange, -log10(padj),color=Cell_type), alpha = 1, size=3)
p <- p + geom_text_repel(data=res2%>% filter(padj < 0.05),aes(log2FoldChange, -log10(padj), label=Gene_group_dcor,color=Cell_type), fontface=4,size = 4) 
p <- p + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#00D500",S="#B326E2"))
p <- p + xlab('Gland Segment/Control (Log2)') 
p <- p + ylab('-log10(adjusted p-value)') 
p <- p + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p

## variance stabilized transformed counts on collapsed replicates
vst<-varianceStabilizingTransformation(ddsColl, blind=FALSE)
vstMat<-assay(vst)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

#principle component analysis
library("genefilter")
library(vegan)

# plotPCA uses top 500 genes by default
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))

# all genes
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsColl)[, c("condition"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = colnames(ddsColl))

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group",label="names")) + 
  geom_point(size = 4,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",seg7="black"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  theme_bw()+
  geom_text_repel()+
  stat_ellipse(type="t", linetype=2)
theme(text = element_text(size=14))+theme(legend.title=element_blank())

###############################################
# calculate fold change of the TPMs by segment
FC_dc<-as.data.frame(vstMat) %>% 
  mutate(avg_ctrl=round((control_1+control_2+control_3)/3,3),
         avg_seg7=round((gland_2+gland_2+gland_3)/3,3),
         ctrl_seg7_ratio=round(avg_seg7-avg_ctrl,3),
         Protein_ID=rownames(.)) %>%
  full_join(anno_dc,by="Protein_ID") %>%
  left_join(res7vC_tab %>% select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% 
  dplyr::rename(beta_seg7vC = log2FoldChange, padj_seg7vC =padj) %>%
  mutate(change=ifelse(beta_seg7vC > 2 & padj_seg7vC < 0.05 , "seg7-enriched",
                       ifelse(beta_seg7vC < -2 & padj_seg7vC < 0.05 , "ctrl-enriched",
                              ifelse(padj_seg7vC < 0.05 , "other","not_sig")))) %>%
  mutate(change=ifelse(is.na(padj_seg7vC), "not_tested", change))

#write.csv(FC_dc,"./08_DEGs_analysis/Dalotia/Bulk_Gland_analysis/FC_enrichment_Dalotia_bulk_RNAseq.csv")

# tabulate number of genes per category above
summary_FC_dc <- FC_dc %>% select(change, Protein_ID) %>%
  group_by(change) %>% dplyr::summarize(n = n())
summary_FC_dc