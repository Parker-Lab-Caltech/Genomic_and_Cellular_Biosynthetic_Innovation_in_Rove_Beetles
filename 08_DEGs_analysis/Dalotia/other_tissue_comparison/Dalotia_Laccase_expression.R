library(pheatmap)
library(viridis)
library("dplyr")
library('tidyr')
library('magrittr')
library('stringr')
library("reshape2")
library("DESeq2")
library("plyr")

#########################
## Tissue-specific CPM ##
#########################
# read in Dalotia smartseq cpm table and filter for enriched genes
dcTis <- read.table("Dcor_otherTissues_gene_counts.txt", 
                    header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="", row.names=1) %>% dplyr::select(-Length)

# coverage filter, 10 reads for at least 4 samples
head(dcTis)
keep <- rowSums(dcTis >= 10) >= 4
dcTis <- dcTis[keep,]

# read in Dalotia smartseq tpm table and filter for enriched genes
dc2 <- read.table("rnaseq_allLibraries.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
dc2<-dc2[c(33:51),c(1:2,5:8)]
head(dc2)
dim(dc2)
dc2$condition<-factor(dc2$condition, level=c("m_body","fem_body","m_antenna","fem_antenna","brain","larvae","pupae","control","seg7"))
dc2$tecRep <-factor(dc2$tecRep)

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsMF<-DESeqDataSetFromMatrix(countData=dcTis,
                              colData= dc2,
                              design = ~ condition)

#collapse replicates
ddsColl <- collapseReplicates(ddsMF, ddsMF$sample, ddsMF$tecRep)

#dataframe information 
colData(ddsColl)

vst<-varianceStabilizingTransformation(ddsColl, blind=FALSE)
vstMat<-as.data.frame(assay(vst))

#average the replicates for the gland, control and brain samples
vstMat$gland_avg<-(vstMat$gland_1_rep1+ vstMat$gland_1_rep2+ vstMat$gland_2_rep1+ vstMat$gland_2_rep2+ vstMat$gland_3_rep1+ vstMat$gland_3_rep2)/6
vstMat$control_avg<-(vstMat$control_1_rep1+ vstMat$control_1_rep2+ vstMat$control_2_rep1+ vstMat$control_2_rep2+ vstMat$control_3_rep1+ vstMat$control_3_rep2)/6

head(vstMat)

# heatmap of laccase genes
genelist<-c("Dcor_evm.model.hic_scaffold_4_pilon.560",
            "Dcor_evm.model.hic_scaffold_1_pilon.330",
            "Dcor_evm.model.hic_scaffold_7_pilon.1222",
            "Dcor_evm.model.hic_scaffold_5_pilon.1298",
            "Dcor_evm.model.hic_scaffold_5_pilon.1299",
            "Dcor_evm.model.hic_scaffold_5_pilon.1300",
            "Dcor_evm.model.hic_scaffold_5_pilon.1302",
            "Dcor_evm.model.hic_scaffold_3_pilon.1487",
            "Dcor_evm.model.hic_scaffold_3_pilon.1486",
            "Dcor_evm.model.hic_scaffold_3_pilon.532"
            
)

lab<-c("Lac1",
       "HAL1",
       "Dmd",
       "HAL3",
       "HAL4",
       "HAL5",
       "HAL6",
       "Lac2A",
       "Lac2B",
       "MCO"
       
)

# filter for gene list
lcpm_gl<- vstMat[genelist,c(1,8,9,16:21)]

# reorder table
lcpm_gl2<- lcpm_gl %>% relocate("larvae", .before="pupae")
lcpm_gl2$gene_id<-rownames(lcpm_gl2)

# set palette
mypalette <- viridis(100)

# plot tissues
mat<-as.matrix(lcpm_gl2[,c(2,4,1,8,9)])
rownames(mat) <- paste(lab)
brks <- seq(7,15,length.out=100) 
p <- pheatmap(t(mat),
              #annotation_col = col.anno, 
              #annotation_colors = col.cell,
              color = mypalette,
              border_color = NA,
              breaks=brks, 
              cluster_cols = FALSE,
              cluster_rows = FALSE)

# plot life stages
mat<-as.matrix(lcpm_gl2[,c(6,7,3,5)])
rownames(mat) <- paste(lab)
p <- pheatmap(t(mat),
              #annotation_col = col.anno, 
              #annotation_colors = col.cell,
              color = mypalette,
              border_color = NA,
              breaks=brks,
              cluster_cols = FALSE,
              cluster_rows = TRUE)

########################
# variability among replicates by tissue type, replicates, etc.
library(forcats)
#normalized counts
dds<-estimateSizeFactors(ddsColl)
cd <- counts(dds,normalized=TRUE)
boxplot(log2(as.matrix(cd)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="DESeq2")

dat <- stack(as.data.frame(cd)) %>% 
  mutate(condition=as.factor(dc2$condition[match(.$ind,dc2$sample)]),
         company=as.factor(dc2$company[match(.$ind,dc2$sample)]))

med.dat<- dat %>% group_by(condition) %>%  dplyr::summarise(median = median(log2(values)+1, na.rm = TRUE))

a<-ggplot(dat %>% dplyr::filter(condition != "control",condition != "seg7" ),
          aes(x = fct_reorder(ind,log2(values)+1,.fun='median'), y = log2(values)+1, fill=condition)) + 
  geom_boxplot()+ xlab("Sample Libraries") + ylab(expression("Log"["2"]~"Normalized Counts"))+
  scale_fill_grey(start = 1,end = 0,)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

# correlation of VST transformed counts, non-gland samples
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(cowplot)

dmat <- as.matrix(cor(vstMat[,c(1,8:9,16:19)],method="spearman"))

col.anno<- data.frame(row.names=colnames(dmat),
                      condition=as.factor(dc2$condition[match(colnames(dmat),dc2$sample)]),
                      company=as.factor(dc2$company[match(colnames(dmat),dc2$sample)]))

col.cell<-list(condition=c(control="grey",D1="#0ECA7C",D2="#B326E2"),
               company=c(`Omega`="black",`Macrogen`="orange"))

mypalette <- rev(magma(100))
b <- pheatmap(dmat,
              main = "Spearman Correlation",
              color = mypalette,
              cluster_cols = TRUE,
              cluster_rows = TRUE,border_color = "grey40")

plot_grid(a,b[[4]],nrow = 1, ncol = 2, rel_widths = c(1, 1, 1))
