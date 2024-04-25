# Venn diagram of gland expression

library(ggplot2)
library(dplyr)

Input <- read.table("./08_DEGs_analysis/DEO_comparison/UpSet_plot/venn_diagram_input_update_noFC_filtering_2.txt",
                    check.names=FALSE, header=T, 
                    na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(Input)

# upper left bar plot
total<-colSums(Input[,10:13])
total<-as.data.frame(total) %>% mutate(species=rep(c("Dalotia","Aleochara"),times=2), cellType=rep(c("D1","D2"),each = 2))%>%
  filter(cellType=="D1" | cellType=="D2")

total$cellType <-factor(total$cellType, levels=c("D2","D1"))

upperleft <- total %>% 
  ggplot(aes(x = interaction(cellType,species), y= total)) +
  geom_bar(stat = "identity", aes(fill = interaction(cellType,species)), position = position_dodge(), width=0.8) +
  geom_text(aes(label = as.character(total)), size = 6, angle = 90, hjust = -0.1, y = 1, fontface = "bold") +
  scale_fill_manual(values = c("#B326E2","#00D500",  "#B326E2","#00D500")) +           
  scale_x_discrete(labels = NULL,  expand = c(0, .2)) +
  scale_y_continuous(labels = NULL, expand=c(0,0)) +
  labs(x = NULL,
       y = "DEOs") +
  theme_minimal() +
  geom_hline(yintercept = -Inf, size = 1.5) +
  theme(legend.position = "none") +
  theme(text = element_text(size= 14, face="bold")) +
  theme(axis.text.x=element_text(colour = "black", angle = 45, hjust = 1)) +
  theme(axis.text.y=element_text(colour = "black")) +
  theme(panel.grid = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))+
  coord_cartesian(clip = 'off')
  
upperleft

# lower right matrix
mat<-Input %>%
  mutate(D1shared=ifelse(Dalotia_D1==1 & Aleochara_D1==1& Dalotia_D2==0 & Aleochara_D2==0, 1,0),
         D2shared=ifelse(Dalotia_D2==1 & Aleochara_D2==1 & Dalotia_D1==0 & Aleochara_D1==0, 1,0),
         Aleo_only=ifelse(Aleochara_D1==1 & Aleochara_D2==1 & Dalotia_D2==0 & Dalotia_D1==0, 1,0),
         Dcor_only=ifelse(Dalotia_D1==1 & Dalotia_D2==1 & Aleochara_D2==0 & Aleochara_D1==0, 1,0),
         mixedDA=ifelse(Dalotia_D1==1 & Aleochara_D2==1 & Dalotia_D2==0 & Aleochara_D1==0, 1,0), 
         mixedAD=ifelse(Dalotia_D2==1 & Aleochara_D1==1 & Dalotia_D1==0 & Aleochara_D2==0, 1,0),
         AleoD1_only=ifelse(Aleochara_D1==1 & Dalotia_D1==0 & Dalotia_D2==0 & Aleochara_D2==0, 1,0),
         AleoD2_only=ifelse(Aleochara_D2==1 & Dalotia_D1==0 & Dalotia_D2==0 & Aleochara_D1==0, 1,0),
         DcorD1_only=ifelse(Dalotia_D1==1 & Aleochara_D1==0 & Aleochara_D2==0 & Dalotia_D2==0, 1,0),
         DcorD2_only=ifelse(Dalotia_D2==1 & Aleochara_D1==0 & Aleochara_D2==0 & Dalotia_D1==0, 1,0),
         D2_3a=ifelse(Dalotia_D2==1 & Aleochara_D2==1 & Dalotia_D1==1  & Aleochara_D1==0 , 1,0),
         D2_3b=ifelse(Dalotia_D2==1 & Aleochara_D2==1 & Dalotia_D1==0  & Aleochara_D1==1 , 1,0),
         D1_3a=ifelse(Dalotia_D1==1 & Aleochara_D1==1 & Dalotia_D2==1  & Aleochara_D2==0 , 1,0),
         D1_3b=ifelse(Dalotia_D1==1 & Aleochara_D1==1 & Dalotia_D2==0  & Aleochara_D2==1 , 1,0),
         allcombo=ifelse(Dalotia_D1==1 & Aleochara_D1==1 & Dalotia_D2==1  & Aleochara_D2==1 , 1,0)) %>%
  select(AleoD1_only,AleoD2_only,DcorD1_only,DcorD2_only,D1shared,D2shared,Aleo_only,
         mixedDA,Dcor_only,mixedAD, D2_3a,D2_3b, D1_3a, D1_3b,allcombo)

mat2<-melt(colSums(mat)) %>% mutate(Set=rownames(.))

mat2$Set<-factor(mat2$Set, levels=rev(c("AleoD2_only","AleoD1_only","DcorD2_only","DcorD1_only","D2shared","D1shared",
                                    "Aleo_only","mixedAD","Dcor_only", "mixedDA","D2_3a","D2_3b","D1_3a","D1_3b","allcombo")))

lowerright <- ggplot(mat2,aes(x = Set, y = value, fill=Set)) +
  geom_bar(stat = "identity",  color = NA, alpha = 0.8) +
  geom_text(aes(label = value, y = value+1), size = 5, hjust = 0, vjust = 0.5, fontface = "bold") +
  scale_y_continuous(labels=NULL, expand=c(0,0)) + #I left white space here for better alignment w/ extended plots
  scale_x_discrete(labels = NULL,expand = c(0, 0)) +
  labs(y = "No. of DEOs",
       x = NULL) +
  scale_fill_manual(values = c("grey80","grey80","grey80","grey80", "grey80","grey80","grey80","grey80" ,
                                       "grey80","#00D500","#B326E2","grey80","grey80" ,"grey80","grey80")) +
  theme_minimal() +
  geom_hline(yintercept = -Inf, size = 1.5) +
  #geom_vline(xintercept = -Inf, size = 1.5) +
  theme(text = element_text(size= 14, face="bold"),legend.position = "none") +
  theme(axis.text.x=element_text(colour = "black", angle = 45, hjust = 1)) +
  theme(axis.text.y=element_text(colour = "black")) +
  theme(panel.grid = element_blank(),plot.margin = unit(c(0,2,0,0), "cm")) + 
  coord_flip(clip = 'off')
lowerright

# lower left matrix
mat3<-read.table("./08_DEGs_analysis/DEO_comparison/UpSet_plot/matrix_pa.txt",header=T, na.strings="na", sep="\t")

mat3$Set<-factor(mat3$Set, levels=rev(c("AleoD2","AleoD1","DcorD2","DcorD1","D2_only","D1_only",
                                        "Aleo_only","mixedAD","Dcor_only", "mixedDA","D2_3a","D2_3b","D1_3a","D1_3b","allcombo")))

mat3<-melt(mat3)

mat3$variable<-factor(mat3$variable, levels=c("Aleo_d2","Aleo_d1","Dcor_d2","Dcor_d1"))

lowerleft<-ggplot(mat3,aes(x = variable,  y=Set, fill=Set, alpha=value))+
  geom_tile(color = "white", size = 1, alpha=c(mat3$value)) +
  scale_y_discrete(labels = NULL,expand = c(0, 0)) +
  scale_x_discrete(labels = NULL,expand = c(0, 0)) + #I left white space here for better alignment w/ extended plots 
  labs(x = " ", #white space for better alignment w/ right side plots 
       y = "overlap") +
  scale_fill_manual(values = c("grey80","grey80","grey80","grey80", "grey80","grey80","grey80","grey80" ,
                                       "grey80", "#00D500","#B326E2","grey80","grey80","grey80","grey80" )) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(text = element_text(size= 14, face="bold")) +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  theme(panel.grid = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))
lowerleft

# putting it together
library(cowplot)
plot_grid(upperleft, NULL,lowerleft, lowerright, 
          nrow = 2, 
          ncol = 2,
          rel_heights = c(2, 4), #the more rows in the lower part, the longer it should be
          rel_widths = c(1, 1))

pdf(paste0("upset_celltype_041823.pdf"),paper = "letter", width = 5, height = 5)
plot_grid(upperleft, NULL,lowerleft, lowerright, 
          nrow = 2, 
          ncol = 2,
          rel_heights = c(2, 4), #the more rows in the lower part, the longer it should be
          rel_widths = c(1, 1))
dev.off()


## install GO database
install.packages("./01_Dalotia_Genome_Assembly_Annotation/org.Dcoriaria5.eg.db", repos=NULL,type="source")

library(clusterProfiler)
library(enrichplot)

tab<-read.table("./08_DEGs_analysis/DEO_comparison/UpSet_plot/venn_diagram_input_update_noFC_filtering_2.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")

# D1- BQ cells
D1<-tab %>% filter(Dalotia_D1=="1"&Aleochara_D1=="1") #238 genes


# Cellular Component
egoCC <- enrichGO(gene         = D1$Dcor_Gene,
                 OrgDb         = "org.Dcoriaria5.eg.db",
                 keyType       = 'GID',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(egoCC)
cc2 <- simplify(egoCC, cutoff=0.7, by="p.adjust", select_fun=min)
cc2
barplot(cc2)

egoCCdat<-as.data.frame(cc2)

cc3 <- pairwise_termsim(cc2)
emapplot(cc3,cex_category=1.5)

# Biological Process
egoBP <- enrichGO(gene         = D1$Dcor_Gene,
                  OrgDb         = "org.Dcoriaria5.eg.db",
                  keyType       = 'GID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(egoBP)
bp2 <- simplify(egoBP, cutoff=0.7, by="p.adjust", select_fun=min)
bp2
barplot(bp2)

egoBPdat<-as.data.frame(bp2)

bp3 <- pairwise_termsim(bp2)
emapplot(bp3,cex_category=1.5)

# Molecular Function
egoMF <- enrichGO(gene         = D1$Dcor_Gene,
                  OrgDb         = "org.Dcoriaria5.eg.db",
                  keyType       = 'GID',
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(egoMF)
mf2 <- simplify(egoMF, cutoff=0.7, by="p.adjust", select_fun=min)
mf2
barplot(mf2)

egoMFdat<-as.data.frame(mf2)

mf3 <- pairwise_termsim(mf2)
emapplot(mf3,cex_category=1.5)

alldat<-rbind(egoCCdat,egoBPdat,egoMFdat)

#write.table(alldat,"D1_GO_enrichment.txt",sep="\t")

# D2- Solvent cells
D2<-tab %>% filter(Dalotia_D2=="1"& Aleochara_D2=="1") #364 genes

# DEGs of Dalotia with >2-fold upregulation in D2
dc_D2<-FC_dc %>% filter(beta_D2vC >=2 & padj_D2vC < 0.05)
# DEGs of Aleochara with >2-fold upregulation in D2
al_D2<-FC_al %>% filter(beta_D2vC >=2 & padj_D2vC < 0.05) %>% filter(!is.na(Dcor_Gene))

# Cellular Component
egoCC <- enrichGO(gene         = D2$Dcor_Gene,
                  OrgDb         = "org.Dcoriaria5.eg.db",
                  keyType       = 'GID',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(egoCC)
cc2 <- simplify(egoCC, cutoff=0.7, by="p.adjust", select_fun=min)
cc2
barplot(cc2)
barplot(egoCC)

egoCCdat<-as.data.frame(cc2)

cc3 <- pairwise_termsim(cc2)
emapplot(cc3,cex_category=1.5)

# Biological Process
egoBP <- enrichGO(gene         = D2$Dcor_Gene,
                  OrgDb         = "org.Dcoriaria5.eg.db",
                  keyType       = 'GID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(egoBP)
bp2 <- simplify(egoBP, cutoff=0.7, by="p.adjust", select_fun=min)
bp2
barplot(bp2)
barplot(egoBP)

egoBPdat<-as.data.frame(bp2)

bp3 <- pairwise_termsim(bp2)
emapplot(bp3,cex_category=1.5)

# Molecular Function
egoMF <- enrichGO(gene         = D2$Dcor_Gene,
                  OrgDb         = "org.Dcoriaria5.eg.db",
                  keyType       = 'GID',
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(egoMF)
mf2 <- simplify(egoMF, cutoff=0.7, by="p.adjust", select_fun=min)
mf2
barplot(mf2)
barplot(egoMF)

egoMFdat<-as.data.frame(mf2)

mf3 <- pairwise_termsim(mf2)
emapplot(mf3,cex_category=1.5)

alldat<-rbind(egoCCdat,egoBPdat,egoMFdat)

#write.table(alldat,"D2_GO_enrichment.txt",sep="\t")



##############
anno<-read.table("./01_Dalotia_Genome_Assembly_Annotation/Dcor_GeneAnnotation_combined_v3_2022.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")

head(anno)

ktab_D1<-D1 %>%
  left_join(anno %>% select(Protein_ID,KEGG_orthology,GO_terms_eggNOG,GO_terms_other),by=c("Dcor_Gene"="Protein_ID"))

# KEGG enrichment (need gene list with KO ids)
kk <- enrichKEGG(gene         = ktab_D1$KEGG_orthology,
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk)

kkdat_D1<-as.data.frame(kk)

write.table(kkdat_D1,"D1_KEGG_enrichment.txt",sep="\t")

ktab_D2<-D2 %>%
  left_join(anno %>% select(Protein_ID,KEGG_orthology,GO_terms_eggNOG,GO_terms_other),by=c("Dcor_Gene"="Protein_ID"))

kk2 <- enrichKEGG(gene         = ktab_D2$KEGG_orthology,
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk2)

kkdat_D2<-as.data.frame(kk2)

write.table(kkdat_D2,"D2_KEGG_enrichment.txt",sep="\t")