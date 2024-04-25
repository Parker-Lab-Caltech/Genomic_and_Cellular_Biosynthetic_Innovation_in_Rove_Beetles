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

###########
# Holobus #
###########
# load in metadata for Holobus 
holMeta <- read.table("./08_DEGs_analysis/Holobus/sample_metadata.txt", 
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
# remove head sample
holMeta_2 <- holMeta %>% filter(!sample=="22785", !sample=="23386") 

# factor seg6 and seg7
holMeta_2$condition<-factor(holMeta_2$segment, levels=c("control","seg7"))
#holMeta$segment<-factor(holMeta$segment, levels=c("control","seg7"))
holMeta_2$batch<-factor(holMeta_2$batch)

# dcor to holobus annotation
anno <- read.table("./08_DEGs_analysis/Holobus/Hol.emapper.annotations.txt",
                   header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

dcor_hol<-read.table("./08_DEGs_analysis/Holobus/Dalotia_update.aa__v__Holobus.aa.tsv",
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

sample_tibble2 <- dcor_hol %>%
  group_by(Dalotia_update.aa) %>%
  dplyr::slice(1) %>%
  group_by(Holobus.aa) %>%
  dplyr::slice(1) %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(Hol_Gene = strsplit(as.character(Holobus.aa), ", ")) %>%
  mutate(Dcor_Gene = strsplit(as.character(Dalotia_update.aa), ", ")) %>%
  unnest(Hol_Gene) %>%
  unnest(Dcor_Gene) %>%
  ungroup()

anno<-anno[,c(1,2,8,9)] %>%
  left_join(sample_tibble2 %>% dplyr::select(Hol_Gene,Dcor_Gene), by=c("query"="Hol_Gene")) %>%
  left_join(anno_dc, by=c("Dcor_Gene"="Protein_ID")) %>%
  mutate(Annotation=ifelse(is.na(Annotation),Preferred_name,Annotation)) %>%
  group_by(query) %>% 
  dplyr::slice(1) %>%
  ungroup()

ho<- read.delim("./08_DEGs_analysis/Holobus/Hol_smartseq_gene_counts.txt",
                head=T, row.names=1,check.names = FALSE) 

ho2<-ho[,c(-1,-7,-11)]
head(ho2)

# coverage filter, 10 reads for at least 4 samples
keep <- rowSums(ho2 >= 10) >= 4
ho2 <- ho2[keep,]

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsHO<-DESeqDataSetFromMatrix(countData=ho2,
                              colData= holMeta_2,
                              design = ~ batch + condition)

#dataframe information 
colData(ddsHO)

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsHO<-DESeq(ddsHO, betaPrior=FALSE)

# size factor added
colData(ddsHO)

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsHO)

# gland cells v control cells
HO_res7vC<-results(ddsHO, name="segment_seg7_vs_control",alpha=0.05,cooksCutoff=TRUE)
summary(HO_res7vC,alpha=0.05)
HO_res7vC_tab<-as.data.frame(HO_res7vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno %>% select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
HO_res7vC_tabOrder<-HO_res7vC_tab[order(HO_res7vC_tab$padj),]
head(HO_res7vC_tabOrder)
HO_resSig7vC<-subset(HO_res7vC_tab, padj <0.05) 
#write.csv(HO_res7vC_tab, file= "Holobus_BULK_seg7vC_tab_ALL.csv")

###################################################
# variancePartition analysis
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(ddsHO)>1) >= 0.5 * ncol(ddsHO)

#compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( ddsHO )[isexpr,] + 1)

# Define formula
form <- ~  (1|batch) + (1|reads) + (1|condition)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, holMeta_2)

# plot results
vp <- sortCols(varPart, FUN="mean")
a<-plotVarPart(vp, col=c("grey90","grey50","black","white"))

# get median % variation by random effect
summary(vp)

# variability among replicates by cell type, batch, etc.
# extract normalized read counts
cd <- counts(ddsHO,normalized=TRUE)

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(holMeta_2$condition[match(.$ind,holMeta_2$sample)]),
         batch=as.factor(holMeta_2$batch[match(.$ind,holMeta_2$sample)]),
         reads=as.factor(holMeta_2$reads[match(.$ind,holMeta_2$sample)]))

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

vst<-varianceStabilizingTransformation(ddsHO, blind=FALSE)
vstMat<-assay(vst)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

# plotPCA from DESeq2 uses top 500 genes by default
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))

# PC variance explained
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsHO)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d_500 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, names = colnames(ddsHO))

c<- ggplot(data = d_500, aes_string(x = "PC1", y = "PC2", color = "condition", shape="batch")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",seg7="#0ECA7C"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# repeat but will all genes
pcaALL <- prcomp(t(assay(vst)))

# PC variance explained
percentVarALL <- pcaALL$sdev^2/sum(pcaALL$sdev^2)

intgroup.df <- as.data.frame(colData(ddsHO)[, c("condition","batch"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
dALL <- data.frame(PC1 = pcaALL$x[, 1], PC2 = pcaALL$x[, 2], group = group, 
                   intgroup.df, names = colnames(ddsHO))

d<- ggplot(data = dALL, aes_string(x = "PC1", y = "PC2", color = "condition", shape="batch")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",seg7="#0ECA7C"))+
  xlab(paste0("PC1: ", round(percentVarALL[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVarALL[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# put all the pieces together
plot_grid(a,b,c,d,nrow = 1, ncol = 4, rel_widths = c(0.3,0.4,0.25,0.25))


###############################################
# calculate fold change of the TPMs by segment
FC_ho<-as.data.frame(vstMatHO) %>% 
  mutate(avg_ctrl=round((`22784`+`23406`+`H13`+`H14`)/4,3),
         avg_seg7=round((`22783`+`23404`+`H11`+`H12`)/4,3),
         ctrl_seg7_ratio=round(avg_seg7-avg_ctrl,3),
         Protein_ID=rownames(.)) %>%
  full_join(anno %>% select(query,Dcor_Gene, Annotation,Description),by=c("Protein_ID"="query")) %>%
  left_join(HO_res7vC_tab %>% select(Protein_ID,log2FoldChange,padj), by=c("Protein_ID")) %>% 
  dplyr::rename(beta_seg7vC = log2FoldChange, padj_seg7vC =padj) %>%
  mutate(change=ifelse(beta_seg7vC > 2 & padj_seg7vC < 0.05 , "seg7-enriched",
                       ifelse(beta_seg7vC < -2 & padj_seg7vC < 0.05 , "ctrl-enriched",
                              ifelse(padj_seg7vC < 0.05 , "other","not_sig")))) %>%
  mutate(change=ifelse(is.na(padj_seg7vC), "not_tested", change))

#write.csv(FC_ho,"./08_DEGs_analysis/Holobus/FC_enrichment_Holobus_bulk_RNAseq.csv")

# tabulate number of genes per category above
summary_FC_ho <- FC_ho %>% select(change, Protein_ID) %>%
  group_by(change) %>% dplyr::summarize(n = n())
summary_FC_ho