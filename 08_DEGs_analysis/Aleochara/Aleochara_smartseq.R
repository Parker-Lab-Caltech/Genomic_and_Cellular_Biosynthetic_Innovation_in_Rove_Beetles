#####################
# Aleochara SMARTseq #
######################

# load in annotation for Aleochara
anno <- read.table("./08_DEGs_analysis/Aleochara/Alus.emapper.annotations.txt",
                   header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

dcor_alus<-read.table("./08_DEGs_analysis/Aleochara/Dalotia_update.aa__v__Alustrica.aa.tsv",
                      header = TRUE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")

sample_tibble2 <- dcor_alus %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(Alus_Gene = strsplit(as.character(Alustrica.aa), ", ")) %>%
  mutate(Dcor_Gene = strsplit(as.character(Dalotia_update.aa), ", ")) %>%
  unnest(Alus_Gene) %>%
  unnest(Dcor_Gene) %>%
  ungroup()

anno<-anno[,c(1,2,9,11)] %>%
  left_join(sample_tibble2 %>% dplyr::select(Alus_Gene,Dcor_Gene), by=c("query"="Alus_Gene")) %>%
  left_join(anno_dc, by=c("Dcor_Gene"="Protein_ID")) %>%
  mutate(Annotation=ifelse(is.na(Annotation),Tcas_annotation,Annotation)) %>%
  mutate(Annotation=ifelse(is.na(Annotation),Preferred_name,Annotation)) %>%
  group_by(query) %>% 
  dplyr::slice(1) %>%
  ungroup()

#read in counts table, header 
al<- read.delim("./08_DEGs_analysis/Aleochara/Alus_smartseq_gene_counts.txt",
                head=T, row.names=1,check.names = FALSE) 

al2<-al[,c(-1,-5,-11)]
head(al2)

# coverage filter, 10 reads for at least 4 samples
keep <- rowSums(al2 >= 10) >= 4
al2 <- al2[keep,]

# read in Aleochara smartseq tpm table and filter for enriched genes
alMeta <- read.table("./08_DEGs_analysis/Aleochara/sample_metadata.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(alMeta)

alMeta$condition<-factor(alMeta$condition, level=c("control","D1","D2"))

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsAL<-DESeqDataSetFromMatrix(countData=al2,
                              colData= alMeta,
                              design = ~ condition)

#dataframe information 
colData(ddsAL)

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsAL<-DESeq(ddsAL, betaPrior=FALSE)

# size factor added
colData(ddsAL)

#information on model tested by DESeq, should be standard
attr(ddsAL, "modelMatrixType")
attr(ddsAL, "modelMatrix")

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsAL)

# D1 cells v control cells
AL_resD1vC<-results(ddsAL, name="condition_D1_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(AL_resD1vC,alpha=0.05)
AL_resD1vC_tab<-as.data.frame(AL_resD1vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno %>% dplyr::select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
AL_resD1vC_tabOrder<-AL_resD1vC_tab[order(AL_resD1vC_tab$padj),]
head(AL_resD1vC_tabOrder)
AL_resSigD1vC<-subset(AL_resD1vC_tab, padj <0.05) 
#write.csv(AL_resD1vC_tab, file= "Alus_D1vC_tab_ALL.csv")

# D2 cells v control cells
AL_resD2vC<-results(ddsAL, name="condition_D2_vs_control",alpha=0.05,cooksCutoff=FALSE)
summary(AL_resD2vC,alpha=0.05)
AL_resD2vC_tab<-as.data.frame(AL_resD2vC)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno %>% dplyr::select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
AL_resD2vC_tabOrder<-AL_resD2vC_tab[order(AL_resD2vC_tab$padj),]
head(AL_resD2vC_tabOrder)
AL_resSigD2vC<-subset(AL_resD2vC_tab, padj <0.05) 
#write.csv(AL_resD2vC_tab, file= "Alus_D2vC_tab_ALL.csv")

# D2 cells v D1 cells
AL_resD1vD2<-results(ddsAL, contrast=list("condition_D2_vs_control","condition_D1_vs_control"),alpha=0.05,cooksCutoff=FALSE)
summary(AL_resD1vD2, alpha=0.05)
AL_resD1vD2_tab<-as.data.frame(AL_resD1vD2)%>% mutate(Protein_ID=rownames(.)) %>% left_join(anno %>% dplyr::select(query, Annotation,Dcor_Gene),by=c("Protein_ID"="query"))
AL_resD1vD2_tabOrder<-AL_resD1vD2_tab[order(AL_resD1vD2_tab$padj),]
head(AL_resD1vD2_tabOrder)
AL_resSigD1vD2<-subset(AL_resD1vD2, padj <0.05) 
#write.csv(AL_resD1vD2_tab, file= "Alus_D1vD2_tab_ALL.csv")

###################################################
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(ddsAL)>1) >= 0.5 * ncol(ddsAL)

#compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( ddsAL )[isexpr,] + 1)

# Define formula
form <- ~  (1|ext_method) + (1|condition)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, alMeta)

# plot % variation
vp <- sortCols( varPart )
a<-plotVarPart(vp, col=c("grey90","black","white"))

# get median % variation
summary(vp)

# variability among replicates by cell type, batch, etc.
#normalized counts
cd <- counts(ddsAL,normalized=TRUE)

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(alMeta$condition[match(.$ind,alMeta$sample)]),
         ext_method=as.factor(alMeta$ext_method[match(.$ind,alMeta$sample)]))

med.dat<- dat %>% 
  group_by(ind) %>% 
  dplyr::mutate(varcount =var(values, na.rm = TRUE),median = median(log2(values)+1, na.rm = TRUE),
                avg=mean(values, na.rm = TRUE)) %>%
  mutate(varcount=log2(varcount)+1, avg=log2(avg)+1) %>%
  group_by(condition) %>% 
  dplyr::select(ind,varcount,median,avg) %>%
  distinct()

ggplot(med.dat, aes(x=condition,y=varcount))+
  geom_boxplot()+geom_point() 

b<-ggplot(dat,aes(x = fct_reorder(ind,log2(values)+1,.fun='median'), y = log2(values)+1, fill=ext_method)) + 
  geom_boxplot(lwd=0.1, outlier.size=0.1)+ xlab("Sample Libraries") + ylab(expression("Log"["2"]~"Normalized Counts"))+
  scale_fill_grey(start = 1,end = 0,)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")+
  facet_grid(~ condition, scales="free")

#principal component analysis on VST transformed counts
# vst counts
vstAL<-varianceStabilizingTransformation(ddsAL, blind=FALSE)
vstMatAL<-assay(vstAL)
#write.csv(as.data.frame(vstMat), file= "vsd_star.csv")

# plotPCA from DESeq2 uses top 500 genes by default
rv <- rowVars(assay(vstAL))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vstAL)[select, ]))

# PC variance explained
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsAL)[, c("condition","ext_method"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d_500 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, names = colnames(ddsAL))

c<- ggplot(data = d_500, aes_string(x = "PC1", y = "PC2", color = "condition", shape="ext_method")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

# repeat but will all genes
pcaALL <- prcomp(t(assay(vstAL)))

# PC variance explained
percentVarALL <- pcaALL$sdev^2/sum(pcaALL$sdev^2)

intgroup.df <- as.data.frame(colData(ddsAL)[, c("condition","ext_method"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
dALL <- data.frame(PC1 = pcaALL$x[, 1], PC2 = pcaALL$x[, 2], group = group, 
                   intgroup.df, names = colnames(ddsAL))

d<- ggplot(data = dALL, aes_string(x = "PC1", y = "PC2", color = "condition", shape="ext_method")) + 
  geom_point(size = 3,alpha = 0.8) + 
  scale_color_manual(values=c(control="grey",D1="#0ECA7C",D2="#B326E2"))+
  xlab(paste0("PC1: ", round(percentVarALL[1] *100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVarALL[2] *100), "% variance")) + 
  stat_ellipse(type="t", linetype=2) + theme_bw()+
  theme(text = element_text(size=10))+theme(legend.title=element_blank(),
                                            legend.position = "bottom") 

plot_grid(a,b,c,d,nrow = 1, ncol = 4, rel_widths = c(0.3,0.4,0.25,0.25))

# correlation of VST transformed counts
dmat <- as.matrix(cor(vstMatAL,method="spearman"))

col.anno<- data.frame(row.names=colnames(dmat),
                      condition=as.factor(alMeta$condition[match(colnames(dmat),alMeta$sample)]),
                      ext_method=as.factor(alMeta$ext_method[match(colnames(dmat),alMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               ext_method=c(NEBNext="black",TRIZOL="orange"))

mypalette <- rev(magma(100))
p <- pheatmap(dmat,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")
###########
## Volcano plots
vhigh<-read.table("./08_DEGs_analysis/Dcor_Alus_Liom_highlightVolcano_genes.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t")

# add column to indicate significant result
AL_res <- dplyr::mutate(AL_resD1vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
AL_res2<-AL_resD1vC_tab[AL_resD1vC_tab$Protein_ID %in% vhigh$Alus_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Alus_gene_id,Gene_group_alus, Cell_type_alus), by=c("Protein_ID"="Alus_gene_id")) %>%
  filter(Cell_type_alus=="BQ")

# plot all things
p <- ggplot(AL_res, aes(log2FoldChange, -log10(padj))) +theme_classic()
p <- p + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
p <- p + geom_point(data=AL_res2 ,aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 0.5, size=1.5) 
p <- p + geom_point(data=AL_res2 %>% filter(!is.na(Gene_group_alus)),aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 1, size=3) 
p <- p + geom_text_repel(data=AL_res2,aes(log2FoldChange, -log10(padj), label=Gene_group_alus,color=Cell_type_alus), fontface=4,size = 4, max.overlaps=20)
p <- p + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
p <- p + xlab('BQ cells/Control cells (Log2)') + xlim(-30,20)+ylim(0,150)
p <- p + ylab('-log10(adjusted p-value)') 
p <- p + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p

# add column to indicate significant result
AL_res3 <- dplyr::mutate(AL_resD2vC_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
AL_res4<-AL_resD2vC_tab[AL_resD2vC_tab$Protein_ID %in% vhigh$Alus_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Alus_gene_id,Gene_group_alus, Cell_type_alus), by=c("Protein_ID"="Alus_gene_id")) %>%
  filter(Cell_type_alus=="S")

# plot all things D2 v control
q <- ggplot(AL_res3, aes(log2FoldChange, -log10(padj))) +theme_classic()
q <- q + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
q <- q + geom_point(data=AL_res4 ,aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 0.5, size=1.5) 
q <- q + geom_point(data=AL_res4 %>% filter(!is.na(Gene_group_alus)),aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 1, size=3) 
q <- q + geom_text_repel(data=AL_res4,aes(log2FoldChange, -log10(padj), label=Gene_group_alus,color=Cell_type_alus), fontface=4,size = 4, max.overlaps=20)
q <- q + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
q <- q + xlab('Solvent cells/Control cells (Log2)')+xlim(-25,15) 
q <- q + ylab('-log10(adjusted p-value)') 
q <- q + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

q

#grid.arrange(p,q,nrow=1)
#ggarrange(p, q + rremove("ylab"), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

# D1 v D2
# add column to indicate significant result
AL_res5 <- dplyr::mutate(AL_resD1vD2_tab, significant = padj < 0.05)

# extract just the biosynthesis genes of interests
AL_res6<-AL_resD1vD2_tab[AL_resD1vD2_tab$Protein_ID %in% vhigh$Alus_gene_id,] %>% 
  left_join(vhigh %>% dplyr::select(Alus_gene_id,Gene_group_alus, Cell_type_alus), by=c("Protein_ID"="Alus_gene_id")) %>%
  distinct()

r <- ggplot(AL_res5, aes(log2FoldChange, -log10(padj))) +theme_classic()
r <- r + geom_point(aes(color=significant), alpha = 0.3, shape=19) 
r <- r + geom_point(data=AL_res6 ,aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 0.5, size=1.5) 
r <- r + geom_point(data=AL_res6 %>% filter(!is.na(Gene_group_alus)),aes(log2FoldChange, -log10(padj),color=Cell_type_alus), alpha = 1, size=3) 
r <- r + geom_text_repel(data=AL_res6,aes(log2FoldChange, -log10(padj), label=Gene_group_alus,color=Cell_type_alus), fontface=4,size = 4,max.overlaps=15) 
r <- r + scale_color_manual(values = c(`FALSE`='grey95',`TRUE`='grey85',BQ="#0ECA7C",S="#B326E2",G="#FEF02E"))
r <- r + xlab('Solvent cells/BQ cells (Log2)')+xlim(-15,15)
r <- r + ylab('-log10(adjusted p-value)') 
r <- r + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
r

# solvent cell gene expression profile (GEP) genes, RUN DALOTIA first!!
geps<-read.table("./08_DEGs_analysis/GEP9_17_Dalotia.txt", 
                 header = TRUE, stringsAsFactors=FALSE, sep="\t")

gepsFC_al<-AL_resD2vC_tab[AL_resD2vC_tab$Dcor_Gene %in% geps$Dcor_v3_gene_id,] %>% 
  mutate(significant = padj < 0.05) %>% 
  left_join(geps %>% dplyr::select(Dcor_v3_gene_id,GEP_id), by=c("Dcor_Gene"="Dcor_v3_gene_id"))


newdata17_al <- gepsFC_al[order(-gepsFC_al$log2FoldChange),] %>% filter(GEP_id =="GEP17") %>% filter(significant =="TRUE") %>%
  left_join(newdata17 %>% dplyr::select(rankID, Protein_ID,log2FoldChange), by=c("Dcor_Gene"="Protein_ID"))  %>%  filter(!is.na(log2FoldChange.y)) %>%
  mutate(diffDcor=log2FoldChange.x-log2FoldChange.y)

newdata9_al <- gepsFC_al[order(-gepsFC_al$log2FoldChange),] %>% filter(GEP_id =="GEP9") %>%  filter(significant =="TRUE") %>%
  left_join(newdata9 %>% dplyr::select(rankID, Protein_ID,log2FoldChange), by=c("Dcor_Gene"="Protein_ID")) %>%  filter(!is.na(log2FoldChange.y)) %>%
  mutate(diffDcor=log2FoldChange.x-log2FoldChange.y)

#dc_al_combo<-AL_resD2vC_tab %>% filter(padj  <=0.05) %>%
#  left_join(resD2vC_tab %>% filter(padj  <=0.05) %>% dplyr::select(Protein_ID,log2FoldChange), by=c("Dcor_Gene"="Protein_ID")) %>%
#  filter(!is.na(log2FoldChange.y))

dc_al_combo<-AL_resD2vC_tab %>%
  left_join(resD2vC_tab %>% dplyr::select(Protein_ID,log2FoldChange), by=c("Dcor_Gene"="Protein_ID")) %>%
  filter(!is.na(log2FoldChange.y)) %>%
  left_join(gepsFC_al %>% dplyr::select(Dcor_Gene,GEP_id), by=c("Dcor_Gene")) %>%
  distinct()


a<-ggplot(newdata17,aes(x=rankID, y=log2FoldChange))+
  geom_point(data=newdata17_al ,aes(x=rankID, y=log2FoldChange.x), color="grey50", alpha=0.8,size=2)+
  geom_segment(data=newdata17_al, aes(x=rankID, xend=rankID, y=log2FoldChange.y, yend=log2FoldChange.x),color="grey50", alpha=0.8)+
  geom_point(alpha=0.3,size=1,color="black")+ 
  theme_bw()+xlab("") +ylab("log2 Fold Change")+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #geom_point(data=newdata17_al %>% filter(log2FoldChange.y >2 &log2FoldChange.x >2),aes(x=rankID, y=log2FoldChange.x), shape=8,size=4)+
  #geom_text_repel(data=newdata17_al %>% filter(log2FoldChange.y >2 &log2FoldChange.x >2),aes(x=rankID, y=log2FoldChange.x, label=Annotation))+
  scale_y_continuous(limits=c(-10,10),expand=c(0,0))

b<-ggplot(newdata9,aes(x=rankID, y=log2FoldChange))+
  geom_point(data=newdata9_al ,aes(x=rankID, y=log2FoldChange.x), color="red", alpha=0.8,size=2)+
  geom_segment(data=newdata9_al, aes(x=rankID, xend=rankID, y=log2FoldChange.y, yend=log2FoldChange.x),color="red", alpha=0.8)+
  geom_point(alpha=0.3,size=1,color="black")+ 
  theme_bw()+xlab("") +ylab("log2 Fold Change") +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #geom_point(data=newdata9_al %>% filter(log2FoldChange.y >2 &log2FoldChange.x >2),aes(x=rankID, y=log2FoldChange.x), shape=8,size=4)+
  #geom_text_repel(data=newdata9_al %>% filter(log2FoldChange.y >2 &log2FoldChange.x >2),aes(x=rankID, y=log2FoldChange.x, label=Annotation))+
  scale_y_continuous(limits=c(-10,10),expand=c(0,0))

c<-ggplot(newdata9_al,aes(x=GEP_id, y=diffDcor))+ geom_hline(yintercept = 0) +geom_violin(fill="red",alpha=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.5, color="red")+
  theme_bw() + ylim(-10,10)+xlab("")+scale_y_continuous(position = "right", limits=c(-10,10),expand=c(0,0))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line.y.right = element_line(size = 0.5))
d<-ggplot(newdata17_al,aes(x=GEP_id, y=diffDcor))+ geom_hline(yintercept = 0) + geom_violin(fill="grey50",alpha=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=0.5, color="grey50")+
  theme_bw()+xlab("")+ scale_y_continuous(position = "right", limits=c(-10,10),expand=c(0,0))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line.y.right = element_line(size = 0.5))

grid.arrange(b,c,a,d, ncol=2, nrow=2, widths=c(4,1))

# correlation test
cor.test(newdata9_al$log2FoldChange.x,newdata9_al$log2FoldChange.y, method = "spearman")
cor.test(newdata17_al$log2FoldChange.x,newdata17_al$log2FoldChange.y, method = "spearman")

# correlation of fold change values
ggplot(dc_al_combo, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=factor(GEP_id)))+
  geom_point(color="grey50", alpha=0.1,size=2)+
  geom_point(data=newdata17_al, aes(x=log2FoldChange.x, y=log2FoldChange.y), color="black", alpha=0.3,size=3)+
  geom_point(data=newdata9_al, aes(x=log2FoldChange.x, y=log2FoldChange.y),color="red", alpha=0.3,size=3)+
  scale_color_manual(values=c(GEP17="grey50", GEP9="red"),na.translate = F)+
  theme_bw() +xlab("Aleochara log2 Fold Change") +ylab("Dalotia log2 Fold Change")+
  geom_smooth(method = "lm", se = FALSE)

#write.csv(gepsFC,"./08_DEGs_analysis/Aleochara/Alus_GEP9_GEP17_D2vC_expression.txt")

###############################################
# calculate fold change of the TPMs by cell type
FC_al<-as.data.frame(vstMatAL) %>% 
  mutate(avg_ctrl=round((.$`23391`+.$`23393`+.$`23394`+.$`23395`)/4,3),
         avg_D1=round((.$`23388`+.$`23385`+.$`23390`+.$`23389`+.$`23387`)/5,3),
         avg_D2=round((.$`23400`+.$`23396`+.$`23398`+.$`23399`+.$`23397`)/5,3),
         ctrl_d1_ratio=round(avg_D1-avg_ctrl,3), ctrl_d2_ratio=round(avg_D2-avg_ctrl,3),
         d1_d2_ratio=round(avg_D1-avg_D2,3),
         Protein_ID=rownames(.)) %>%
  full_join(al %>% mutate(Protein_ID=rownames(.)) %>% dplyr::select(Protein_ID), by="Protein_ID")%>%
  full_join(anno %>% select(query,Dcor_Gene, Annotation),by=c("Protein_ID"="query")) %>%
  left_join(AL_resD1vC_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D1vC = log2FoldChange, padj_D1vC =padj) %>%
  left_join(AL_resD2vC_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vC = log2FoldChange, padj_D2vC =padj) %>%
  left_join(AL_resD1vD2_tab %>% dplyr::select(Protein_ID,log2FoldChange,padj), by="Protein_ID")%>% dplyr::rename(beta_D2vD1 = log2FoldChange, padj_D2vD1 =padj) %>%
  mutate(change=ifelse(beta_D1vC > 2 & d1_d2_ratio > 2 & padj_D1vC < 0.05 , "D1-enriched",
                       ifelse(beta_D2vC > 2 & d1_d2_ratio < -2 & padj_D2vC < 0.05, "D2-enriched",
                              ifelse(beta_D1vC > 2 & beta_D2vC > 2 & d1_d2_ratio > -2 & d1_d2_ratio < 2 & padj_D1vC < 0.05 & padj_D2vC < 0.05,"gland-enriched",
                                     ifelse(beta_D1vC < -2 & beta_D2vC < -2 & padj_D1vC < 0.05 & padj_D2vC < 0.05, "ctrl-enriched",
                                            ifelse(padj_D1vC < 0.05 | padj_D2vC < 0.05 | padj_D2vD1 < 0.05, "other","not_sig")))))) %>%
  mutate(change=ifelse(is.na(padj_D1vC)|is.na(padj_D2vC)|is.na(padj_D2vD1), "not_tested", change))

#write.csv(FC_al,"./08_DEGs_analysis/Aleochara/FC_enrichment_Aleochara.csv")

# tabulate number of genes per category above
summary_FC_al <- FC_al %>% dplyr::select(change, Protein_ID) %>%
  group_by(change) %>% dplyr::summarize(n = n())
summary_FC_al

####################
# plot heat map of significant genes

# create dataframe for heat map
sub_DEGs<- FC_al %>% filter(change == "D1-enriched") %>%
  arrange(desc(beta_D1vC)) %>%
  dplyr::select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
                -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1,-Dcor_Gene) %>%
  distinct(.[,1:14],.keep_all = TRUE)

sub_DEGs<- FC_al %>% filter(change == "D2-enriched") %>%
  arrange(desc(beta_D2vC)) %>%
  dplyr::select(-change,-padj_D1vC,-padj_D2vC,-padj_D2vD1,-d1_d2_ratio, -ctrl_d1_ratio,-ctrl_d2_ratio,-avg_ctrl,
                -avg_D1, -avg_D2, -beta_D1vC,-beta_D2vC,-beta_D2vD1,-Dcor_Gene) %>%
  distinct(.[,1:14],.keep_all = TRUE)

rownames(sub_DEGs)<-sub_DEGs$Protein_ID

dim(sub_DEGs)

# set column annotation
col.anno<- data.frame(row.names=colnames(sub_DEGs[1:75,c(-15,-16)]),
                      condition=as.factor(alMeta$condition[match(colnames(sub_DEGs[,c(-15,-16)]),alMeta$sample)]),
                      segment=as.factor(alMeta$segment[match(colnames(sub_DEGs[,c(-15,-16)]),alMeta$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               segment=c(seg6="grey",seg7="black"))


# diverging palette is good for scaled TPMs
mypalette <- magma(100)

mat<-as.matrix(sub_DEGs[1:75,c(-15,-16)])
rownames(mat) <- paste(rownames(mat),sub_DEGs[1:75,16])
rownames(mat) <- paste(sub_DEGs[1:75,16])

p <- pheatmap(mat,
              annotation_col = col.anno, 
              annotation_colors = col.cell,
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")

##############
# shared OGs #
##############

# extract orthologs with shared enriched cell-type expression for Aleochara and Dalotia
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

#write.table(OGs_dc_alus, "./08_DEGs_analysis/Aleochara/Dcor_Alus_orthologs_enriched_expression.txt", sep="\t", row.names=FALSE)

shared_D1<-OGs_dc_alus %>% filter(Dcor_change == "D1-enriched" & Alus_change == "D1-enriched" ) 
shared_D2<-OGs_dc_alus %>% filter(Dcor_change == "D2-enriched" & Alus_change == "D2-enriched" ) 

#write.table(shared_D1, "./08_DEGs_analysis/Aleochara/Dcor_Alus_orthologs_D1-enriched_expression.txt", sep="\t", row.names=FALSE)
#write.table(shared_D2, "./08_DEGs_analysis/Aleochara/Dcor_Alus_orthologs_D2-enriched_expression.txt", sep="\t", row.names=FALSE)
