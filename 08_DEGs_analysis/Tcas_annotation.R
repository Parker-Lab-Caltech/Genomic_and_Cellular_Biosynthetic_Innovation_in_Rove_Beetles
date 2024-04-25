########################
# Tribolium Annotation #
########################

library("dplyr")

#Tcas ensembl to ncbi annotation
library("biomaRt")
ensembl <- useEnsemblGenomes(biomart = "metazoa_mart")
searchDatasets(ensembl, pattern = "Tribolium")
ensembl_tcas<- useDataset("tcastaneum_eg_gene", mart=ensembl)
listAttributes(ensembl_tcas)

# exported annotation from ENSEMBL online
tcas<-read.table("E:/Caltech/Beetle Genomes/Tcas_gene_description.txt",
                 header = FALSE, stringsAsFactors=FALSE, sep="\t", as.is = TRUE, quote = "")
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "refseq_peptide"),
  values=tcas$V1,
  mart=ensembl_tcas)

# Get NCBI IDs and annotation
library(AnnotationHub)
ah <-AnnotationHub()

ahs <- query(ah, c('NCBI', 'Tribolium castaneum'))
ahss<-ahs[["AH102319"]]

# Join NCBI and ENSEMBL annotations
tcas2<-tcas %>%
  left_join(genes, by=c("V1"="ensembl_gene_id"))

accNum_tcas<-select(ahss, tcas2$refseq_peptide,c("ACCNUM", "GENENAME"), "ACCNUM") %>% drop_na() %>% distinct()

tcas3<-tcas2 %>%
  left_join(accNum_tcas, by=c("refseq_peptide"="ACCNUM"))

#write.table(tcas3,"E:/Caltech/Beetle Genomes/Tcas_REFseq_id_gene_description.txt", sep="\t",row.names=FALSE)