library("dplyr")
library('tidyr')
library('magrittr')
library('stringr')
library("reshape2")
library("DESeq2")
library("plyr")
library("variancePartition")
library(forcats)
library("genefilter")
library(vegan)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(pairwiseAdonis)

####################
# Dalotia SMARTseq #
####################

# load in annotation for Dalotia
anno_dc <- read.table("./08_DEGs_analysis/Dalotia/SMARTseq_Gland_analysis/Dcor_v3_annotation_only.txt",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

#read in counts table, header 
dc<- read.delim("./08_DEGs_analysis/Dalotia/SMARTseq_Gland_analysis/Dcor_smartseq_gene_counts.txt",
                head=T, row.names=1) 

dc2<-dc[,c(2:4,6:34)]
head(dc2)

# coverage filter, 10 reads for at least 10 samples
keep <- rowSums(dc2 >= 10) >= 10
dc2 <- dc2[keep,]

# read in Dalotia smartseq tpm table and filter for enriched genes
dcMeta <- read.table("./08_DEGs_analysis/Dalotia/SMARTseq_Gland_analysis/smartseq_samples_gland_control.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(dcMeta)

dcMeta$condition<-factor(dcMeta$condition, level=c("control","D1","D2"))
dcMeta$batch <-factor(dcMeta$batch)

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsMF<-DESeqDataSetFromMatrix(countData=dc2,
                              colData= dcMeta,
                              design = ~ batch + condition)

#dataframe information 
colData(ddsMF)

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsMF<-DESeq(ddsMF, betaPrior=FALSE)

# size factor added
colData(ddsMF)

#information on model tested by DESeq, should be standard
attr(ddsMF, "modelMatrixType")
attr(ddsMF, "modelMatrix")

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsMF)

# D1 cells v control cells
resD1vC<-results(ddsMF, name="condition_D1_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(resD1vC,alpha=0.05)
resD1vC_tab<-as.data.frame(resD1vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_dc,by="Protein_ID")
resD1vC_tabOrder<-resD1vC_tab[order(resD1vC_tab$padj),]
head(resD1vC_tabOrder)
resSigD1vC<-subset(resD1vC_tab, padj <0.05) 
#write.csv(resD1vC_tab, file= "Dalotia_D1vC_tab_ALL.csv")

# D2 cells v control cells
resD2vC<-results(ddsMF, name="condition_D2_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(resD2vC,alpha=0.05)
resD2vC_tab<-as.data.frame(resD2vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_dc,by="Protein_ID")
resD2vC_tabOrder<-resD2vC_tab[order(resD2vC_tab$padj),]
head(resD2vC_tabOrder)
resSigD2vC<-subset(resD2vC_tab, padj <0.05) 
#write.csv(resD2vC_tab, file= "Dalotia_D2vC_tab_ALL.csv")

# D2 cells v D1 cells
resD1vD2<-results(ddsMF, contrast=list("condition_D2_vs_control","condition_D1_vs_control"),alpha=0.05,cooksCutoff=FALSE)
summary(resD1vD2, alpha=0.05)
resD1vD2_tab<-as.data.frame(resD1vD2)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno_dc,by="Protein_ID")
resD1vD2_tabOrder<-resD1vD2_tab[order(resD1vD2_tab$padj),]
head(resD1vD2_tabOrder)
resSigD1vD2<-subset(resD1vD2, padj <0.05) 
#write.csv(resD1vD2_tab, file= "Dalotia_D1vD2_tab_ALL.csv")

###################################################
# variancePartition analysis
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(ddsMF)>1) >= 0.5 * ncol(ddsMF)

#compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( ddsMF )[isexpr,] + 1)

# Define formula
form <- ~  (1|batch) + (1|condition)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, dcMeta)

# plot results
vp <- sortCols(varPart)
a<-plotVarPart(vp, col=c("grey90","black","white"))

# get median % variation by random effect
summary(vp)

# variability among replicates by cell type, batch, etc.
# extract normalized read counts
cd <- counts(ddsMF,normalized=TRUE)

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(dcMeta$condition[match(.$ind,dcMeta$sample)]),
         batch=as.factor(dcMeta$batch[match(.$ind,dcMeta$sample)]))

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

vst<-varianceStabilizingTransformation(ddsMF, blind=FALSE)
vstMat<-assay(vst)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

# plotPCA from DESeq2 uses top 500 genes by default
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))

# PC variance explained
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsMF)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d_500 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, names = colnames(ddsMF))

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

intgroup.df <- as.data.frame(colData(ddsMF)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
dALL <- data.frame(PC1 = pcaALL$x[, 1], PC2 = pcaALL$x[, 2], group = group, 
                   intgroup.df, names = colnames(ddsMF))

d<- ggplot(data = dALL, aes_string(x = "PC1", y = "PC2", color = "condition", shape="batch")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVarALL[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVarALL[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# correlation of VST transformed counts
dmat <- as.matrix(cor(vstMat,method="spearman"))

col.anno<- data.frame(row.names=colnames(dmat),
                      condition=as.factor(dcMeta$condition[match(colnames(dmat),dcMeta$sample)]),
                      batch=as.factor(dcMeta$batch[match(colnames(dmat),dcMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               batch=c(`1`="black",`2`="orange",`3`="red",`4`="blue"))

mypalette <- rev(magma(100))
p <- pheatmap(dmat,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")

# put all the pieces together
plot_grid(a,b,c,d,nrow = 1, ncol = 4, rel_widths = c(0.3,0.4,0.25,0.25))

###########
# Volcano plots
vhigh<-read.table("./08_DEGs_analysis/Dcor_Alus_Liom_highlightVolcano_genes.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t")

# add column to indicate significant result
res <- dplyr::mutate(resD1vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
res2<-resD1vC_tab[resD1vC_tab$Protein_ID %in% vhigh$Dcor_v3_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Dcor_v3_gene_id,Gene_group_dcor, Cell_type_dcor), by=c("Protein_ID"="Dcor_v3_gene_id")) %>%
  filter(Cell_type_dcor=="BQ") 

# plot all things
p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +theme_classic()
p <- p + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
p <- p + geom_point(data=res2 ,aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 0.5, size=1.5) 
p <- p + geom_point(data=res2 %>% filter(!is.na(Gene_group_dcor)),aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 1, size=3) 
p <- p + geom_text_repel(data=res2,aes(log2FoldChange, -log10(padj), label=Gene_group_dcor,color=Cell_type_dcor), fontface=4,size = 4, max.overlaps=20)
p <- p + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
p <- p + xlab('BQ cells/Control cells (Log2)') + xlim(-30,20)+ylim(0,40)
p <- p + ylab('-log10(adjusted p-value)') 
p <- p + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p


# add column to indicate significant result
res3 <- dplyr::mutate(resD2vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
res4<-resD2vC_tab[resD2vC_tab$Protein_ID %in% vhigh$Dcor_v3_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Dcor_v3_gene_id,Gene_group_dcor, Cell_type_dcor), by=c("Protein_ID"="Dcor_v3_gene_id")) %>%
  filter(Cell_type_dcor=="S")

# plot all things D2 v control
q <- ggplot(res3, aes(log2FoldChange, -log10(padj))) +theme_classic()
q <- q + geom_point(aes(color=significant), alpha = 0.3, shape=19)
q <- q + geom_point(data=res4 ,aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 0.5, size=1.5) 
q <- q + geom_point(data=res4 %>% filter(!is.na(Gene_group_dcor)),aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 1, size=3) 
q <- q + geom_text_repel(data=res4,aes(log2FoldChange, -log10(padj), label=Gene_group_dcor,color=Cell_type_dcor), fontface=4,size = 4, max.overlaps=20) 
q <- q + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
q <- q + xlab('Solvent cells/Control cells (Log2)')+xlim(-11,11) +ylim(0,120)
q <- q + ylab('-log10(adjusted p-value)') 
q <- q + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
q

#grid.arrange(p,q,nrow=1)
#ggarrange(p, q + remove("ylab"), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

# D1 v D2
# add column to indicate significant result
res5 <- dplyr::mutate(resD1vD2_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
res6<-resD1vD2_tab[resD1vD2_tab$Protein_ID %in% vhigh$Dcor_v3_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Dcor_v3_gene_id,Gene_group_dcor, Cell_type_dcor), by=c("Protein_ID"="Dcor_v3_gene_id")) %>%
  distinct(Protein_ID,.keep_all = TRUE) 

r <- ggplot(res5, aes(log2FoldChange, -log10(padj))) +theme_classic()
r <- r + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
r <- r + geom_point(data=res6 ,aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 0.5, size=1.5) 
r <- r + geom_point(data=res6 %>% filter(!is.na(Gene_group_dcor)),aes(log2FoldChange, -log10(padj),color=Cell_type_dcor), alpha = 1, size=3) 
r <- r + geom_text_repel(data=res6,aes(log2FoldChange, -log10(padj), label=Gene_group_dcor,color=Cell_type_dcor), fontface=4,size = 4, max.overlaps=20) 
r <- r + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
r <- r + xlab('Solvent cells/BQ cells (Log2)')+xlim(-26,21) +ylim(0,60)
r <- r + ylab('-log10(adjusted p-value)') 
r <- r + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
r

# solvent cell gene expression profile (GEP) genes
geps<-read.table("./08_DEGs_analysis/GEP9_17_Dalotia.txt", 
                 header = TRUE, stringsAsFactors=FALSE, sep="\t")

gepsFC<-resD2vC_tab[resD2vC_tab$Protein_ID %in% geps$Dcor_v3_gene_id,] %>% 
  mutate(significant = padj < 0.05) %>% 
  left_join(geps %>% dplyr::select(Dcor_v3_gene_id,GEP_id), by=c("Protein_ID"="Dcor_v3_gene_id"))

newdata17 <- gepsFC[order(-gepsFC$log2FoldChange),] %>% filter(GEP_id =="GEP17") %>%  filter(significant =="TRUE") %>%
  mutate(rankID = with_order(order_by =rev(log2FoldChange) , fun = row_number, x = rev(log2FoldChange)))

newdata9 <- gepsFC[order(-gepsFC$log2FoldChange),] %>% filter(GEP_id =="GEP9") %>%  filter(significant =="TRUE") %>%
  mutate(rankID = with_order(order_by =rev(log2FoldChange) , fun = row_number, x = rev(log2FoldChange)))

ggplot(newdata17,aes(x=rankID, y=log2FoldChange))+
  geom_point(alpha=0.5)+
  geom_point(data=newdata9,aes(x=rankID, y=log2FoldChange), color="blue", alpha=0.5)+
  theme_bw()

#write.csv(gepsFC,"Dcor_GEP9_GEP17_D2vC_expression.txt")

###############################################
# calculate fold change of the TPMs by cell type
FC_dc<-as.data.frame(vstMat) %>% 
  mutate(avg_ctrl=round((lib_03+lib_26+lib_27+lib_28+lib_29+lib_30+lib_46+lib_47+lib_48+lib_49+lib_410)/11,3),
         avg_D1=round((lib_02+lib_06+lib_08+lib_11+lib_12+lib_13+lib_15+lib_16+lib_S1+lib_S2)/10,3),
         avg_D2=round((lib_01+lib_411+lib_43+lib_44+lib_45+lib_25+lib_24+lib_23+lib_22+lib_20+lib_07)/11,3),
         ctrl_d1_ratio=round(avg_D1-avg_ctrl,3), ctrl_d2_ratio=round(avg_D2-avg_ctrl,3),
         d1_d2_ratio=round(avg_D1-avg_D2,3),
         Protein_ID=rownames(.)) %>%
  full_join(anno_dc,by="Protein_ID") %>%
  left_join(resD1vC_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D1vC = log2FoldChange, padj_D1vC =padj) %>%
  left_join(resD2vC_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vC = log2FoldChange, padj_D2vC =padj) %>%
  left_join(resD1vD2_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vD1 = log2FoldChange, padj_D2vD1 =padj) %>%
  mutate(change=ifelse(beta_D1vC > 2 & d1_d2_ratio > 2 & padj_D1vC < 0.05 , "D1-enriched",
                       ifelse(beta_D2vC > 2 & d1_d2_ratio < -2 & padj_D2vC < 0.05, "D2-enriched",
                              ifelse(beta_D1vC > 2 & beta_D2vC > 2 & d1_d2_ratio > -2 & d1_d2_ratio < 2 & padj_D1vC < 0.05 & padj_D2vC < 0.05,"gland-enriched",
                                     ifelse(beta_D1vC < -2 & beta_D2vC < -2 & padj_D1vC < 0.05 & padj_D2vC < 0.05, "ctrl-enriched",
                                            ifelse(padj_D1vC < 0.05 | padj_D2vC < 0.05 | padj_D2vD1 < 0.05, "other","not_sig")))))) %>%
  mutate(change=ifelse(is.na(padj_D1vC)|is.na(padj_D2vC)|is.na(padj_D2vD1), "not_tested", change))

#write.csv(FC_dc,"FC_enrichment_Dalotia.csv")

# tabulate number of genes per category above
summary_FC_dc <- FC_dc %>% dplyr::select(change, Protein_ID) %>%
  group_by(change) %>% dplyr::summarize(n = n())
summary_FC_dc

####################
# plot heat map of significant genes
# create data frame for heat map
# D1 DEGs
sub_DEGs<- FC_dc %>% filter(change == "D1-enriched") %>%
  arrange(desc(beta_D1vC)) %>%
  dplyr::select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
                -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1) %>%
  distinct(.[,1:32],.keep_all = TRUE)

rownames(sub_DEGs)<-sub_DEGs$Protein_ID

dim(sub_DEGs)

# set column annotation
col.anno<- data.frame(row.names=colnames(sub_DEGs[1:75,c(-33,-34)]),
                      condition=as.factor(dcMeta$condition[match(colnames(sub_DEGs[,c(-33,-34)]),dcMeta$sample)]),
                      segment=as.factor(dcMeta$segment[match(colnames(sub_DEGs[,c(-33,-34)]),dcMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               segment=c(seg6="grey",seg7="black"))


# diverging palette is good for scaled TPMs
mypalette <- magma(100)

mat<-as.matrix(sub_DEGs[1:75,c(-33,-34)])
rownames(mat) <- paste(rownames(mat),sub_DEGs[1:75,34])
rownames(mat) <- paste(sub_DEGs[1:75,34])

p <- pheatmap(mat,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")

# D2 DEGs
sub_DEGsD2<- FC_dc %>% filter(change == "D2-enriched") %>%
  arrange(desc(beta_D2vC)) %>%
  dplyr::select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
                -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1) %>%
  distinct(.[,1:32],.keep_all = TRUE)

rownames(sub_DEGsD2)<-sub_DEGsD2$Protein_ID

dim(sub_DEGsD2)

# set column annotation
col.anno<- data.frame(row.names=colnames(sub_DEGsD2[1:75,c(-33,-34)]),
                      condition=as.factor(dcMeta$condition[match(colnames(sub_DEGsD2[,c(-33,-34)]),dcMeta$sample)]),
                      segment=as.factor(dcMeta$segment[match(colnames(sub_DEGsD2[,c(-33,-34)]),dcMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               segment=c(seg6="grey",seg7="black"))


# diverging palette is good for scaled CPMs
mypalette <- magma(100)

matD2<-as.matrix(sub_DEGsD2[1:75,c(-33,-34)])
rownames(matD2) <- paste(rownames(mat),sub_DEGsD2[1:75,34])
rownames(matD2) <- paste(sub_DEGsD2[1:75,34])

p <- pheatmap(matD2,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")