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

#################################
# Sex-enriched expression by chromosome

#read in counts table, header 
dc_ot<- read.delim("./08_DEGs_analysis/Dalotia/other_tissue_comparison/Dcor_otherTissues_gene_counts.txt",
                   head=T, row.names=1) 

dc2<-dc_ot[,c(2:8)]
head(dc2)

# read in Dalotia smartseq tpm table and filter for enriched genes
dcMeta <- read.table("./08_DEGs_analysis/Dalotia/other_tissue_comparison/rnaseq_allLibraries.txt", 
                     header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

# look at the generated table
head(dcMeta)

otMeta<-dcMeta[33:39,]

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsOT<-DESeqDataSetFromMatrix(countData=dc2,
                              colData= otMeta,
                              design = ~ condition)

ddsOT <- estimateSizeFactors(ddsOT)
sizeFactors(ddsOT)
normalized_counts <- as.data.frame(assay(normTransform(ddsOT, f = log2, pc = 1)))

cutoff <- 1
# set threshold for MW and FW samples for sex-specific expression
drop <- which(apply(normalized_counts[,1:2], 1, max) < cutoff)
d <- normalized_counts[-drop,] 
dim(d) # number of genes left

scaf_log<-tibble::rownames_to_column(d,'geneId') %>%
  mutate(geneId=gsub("evm.model.", "",geneId)) %>%
  separate(geneId, c("chr","geneNum"), sep="\\.")

head(scaf_log)

list<-c("Dcor_hic_scaffold_1_pilon", "Dcor_hic_scaffold_2_pilon","Dcor_hic_scaffold_3_pilon","Dcor_hic_scaffold_4_pilon",
        "Dcor_hic_scaffold_5_pilon","Dcor_hic_scaffold_6_pilon","Dcor_hic_scaffold_7_pilon","Dcor_hic_scaffold_8_pilon",
        "Dcor_hic_scaffold_9_pilon","Dcor_hic_scaffold_10_pilon")

sum_counts<-scaf_log %>%
  filter(chr %in% list) %>%
  group_by(chr) %>%
  summarize(avg_fw = mean(FW), avg_mw = mean(MW), FC_fm_body= avg_fw-avg_mw, avgCPM=(avg_fw+avg_mw)/2)

# fold-change average log2 counts summarized over each chromosome
p <- ggplot(sum_counts, aes(avgCPM, FC_fm_body, label = chr)) +
  geom_text_repel(data= subset(sum_counts, FC_fm_body > 0.2)) +
  geom_text_repel(data= subset(sum_counts, FC_fm_body < -0.6)) +
  geom_point(fill=ifelse(sum_counts$FC_fm_body > 0.2, "#DA5D42", 
                         ifelse(sum_counts$FC_fm_body < -0.6, "#7fadd1", "grey30")), size=5, shape=21) + 
  theme_bw()+
  ylab("Female/Male Fold-change Avg. Log2 Counts")+ xlab("Average Log2 Counts")
p


sexBiasExp<-scaf_log %>% 
  filter(chr %in% list)%>%
  group_by(chr) %>%
  mutate(numGeneChr=n(), enriched= ifelse(MW > FW + 10, "10-fold_male",
                                          ifelse(MW > FW + 5, "5-fold_male",
                                                 ifelse(MW > FW + 2, "2-fold_male",
                                                        ifelse(FW > MW + 10, "10-fold_female",
                                                               ifelse(FW > MW + 5, "5-fold_female",
                                                                      ifelse(FW > MW + 2, "2-fold_female","unbiased"))))))) %>%
  group_by(chr,enriched) %>%
  summarize(biasSum=n(), per_bias=biasSum/numGeneChr, numGeneChr) %>%
  distinct()

col_bias<-c("#2D5276","#7fadd1","#B2ddeb","grey90","#FCA85E","#DA5D42","#B71126")

sexBiasExp$enriched<-factor(sexBiasExp$enriched, levels=c("10-fold_female","5-fold_female","2-fold_female","unbiased","2-fold_male","5-fold_male", "10-fold_male"))

sexBiasExp$chr<-factor(sexBiasExp$chr, levels=c("Dcor_hic_scaffold_1_pilon", "Dcor_hic_scaffold_2_pilon","Dcor_hic_scaffold_3_pilon","Dcor_hic_scaffold_4_pilon",
                                                "Dcor_hic_scaffold_5_pilon","Dcor_hic_scaffold_6_pilon","Dcor_hic_scaffold_7_pilon","Dcor_hic_scaffold_8_pilon",
                                                "Dcor_hic_scaffold_9_pilon","Dcor_hic_scaffold_10_pilon"))

ggplot(sexBiasExp, aes(fill=enriched, y=biasSum, x=chr)) + 
  geom_bar(position="fill", stat="identity") + theme_bw()+
  scale_fill_manual(values=col_bias) +
  scale_y_continuous(expand=c(0,0))+ xlab("")+ylab("Num. of Genes")+
  theme(axis.text.x = element_text(angle = 30, hjust=1))

#sig over/under representation
tab<-as.data.frame(spread(sexBiasExp[,1:3], key = enriched, value = biasSum))

rownames(tab)<-tab$chr

tab[is.na(tab)] <- 0

tab<-as.array(as.matrix(tab[,-1]))
addmargins(tab)

myChiSq <- chisq.test(tab)
myChiSq

round(myChiSq$observed - myChiSq$expected, 2)

myChiSq$observed
myChiSq$expected
myChiSq$stdres

library(chisq.posthoc.test)
chisq.posthoc.test(tab, method = "fdr")

sexBiasExpOverall<-scaf_log %>% 
  filter(chr %in% list)%>%
  group_by(chr) %>%
  mutate(numGeneChr=n(), enriched= ifelse(MW > FW + 2, "male",
                                        ifelse(FW > MW + 2, "female",
                                              ifelse(FW > MW + 2, "2-fold_female","unbiased")))) %>%
           group_by(chr,enriched) %>%
           summarize(biasSum=n(), per_bias=biasSum/numGeneChr, numGeneChr) %>%
           distinct()

#sig over/under representation
tab<-as.data.frame(spread(sexBiasExp[,1:3], key = enriched, value = biasSum))

rownames(tab)<-tab$chr

tab[is.na(tab)] <- 0

tab<-as.array(as.matrix(tab[,-1]))
addmargins(tab)

myChiSq <- chisq.test(tab)
myChiSq

round(myChiSq$observed - myChiSq$expected, 2)

myChiSq$observed
myChiSq$expected
myChiSq$stdres

library(chisq.posthoc.test)
chisq.posthoc.test(tab, method = "fdr")

#####
# promer dot plots
library(GenomicRanges)

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

mumgp = readDelta("protein_PcognY_DcorChr.delta")
mumgp = readDelta("protein_TcasX_DcorChr.delta")
mumgp = readDelta("protein_PcognX_DcorChr.delta")
mumgp = readDelta("protein_OolensX_DcorChr.delta")

filterMum <- function(df, minl=200, flanks=1000){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

chrom<-c("hic_scaffold_1_pilon","hic_scaffold_2_pilon","hic_scaffold_3_pilon","hic_scaffold_4_pilon","hic_scaffold_5_pilon",
         "hic_scaffold_6_pilon","hic_scaffold_7_pilon","hic_scaffold_8_pilon","hic_scaffold_9_pilon","hic_scaffold_10_pilon")

mumgp.filt <- filterMum(mumgp, minl=200) %>%
  filter(qid %in% chrom)

mumgp.filt$qid<- factor(mumgp.filt$qid, levels=chrom)

ggplot(mumgp.filt, aes(x=rs/1000000, xend=re/1000000, y=qs/1000000, yend=qe/1000000, colour=strand)) + geom_segment() +
  geom_point(alpha=0.5) + facet_grid(rows=mumgp.filt$qid, scales = "free_y")+
  theme_bw() + theme(strip.text.y=element_text(angle=0, size=8),
                     legend.position="none", 
                     strip.background=element_blank()) +
  xlab('P. cognatus X chromosome') + ylab('Dalotia chromosomes') + scale_colour_manual(values=c("#7fadd1","#DA5D42"))

ggplot(mumgp.filt, aes(x=rs/1000000, xend=re/1000000, y=qs/1000000, yend=qe/1000000, colour=strand)) + geom_segment() +
  geom_point(alpha=0.5) + facet_grid(rows=mumgp.filt$qid, scales = "free_y")+
  theme_bw() + theme(strip.text.y=element_text(angle=0, size=8),
                     legend.position="none", 
                     strip.background=element_blank()) +
  xlab('O. olens X chromosome') + ylab('Dalotia chromosomes') + scale_colour_manual(values=c("#7fadd1","#DA5D42"))

ggplot(mumgp.filt, aes(x=rs/1000000, xend=re/1000000, y=qs/1000000, yend=qe/1000000, colour=strand)) + geom_segment() +
  geom_point(alpha=0.5) + facet_grid(rows=mumgp.filt$qid, scales = "free_y")+
  theme_bw() + theme(strip.text.y=element_text(angle=0, size=8),
                     legend.position="none", 
                     strip.background=element_blank()) +
  xlab('T. castaneum X chromosome') + ylab('Dalotia chromosomes') + scale_colour_manual(values=c("#7fadd1","#DA5D42"))