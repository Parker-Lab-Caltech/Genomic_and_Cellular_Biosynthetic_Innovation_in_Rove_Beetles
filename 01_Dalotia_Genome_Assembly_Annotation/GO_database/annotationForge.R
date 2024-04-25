########################################
## Not necessary to re-run
library(dplyr)
library(tidyr)

#Convert GO table into correct format
dGO<-read.table("go_terms_4.txt",
                check.names=FALSE, header=F, 
                na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(dGO)

colnames(dGO)<-c("GID","GO_terms")

sample_tibble <- dGO %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  dplyr::mutate(GO = strsplit(GO_terms, ",")) %>%
  unnest(GO) %>%
  dplyr::mutate(EVIDENCE="blast")%>%
  dplyr::select(-group, -GO_terms) %>%
  dplyr::filter(GO != 'NA')%>%
  dplyr::distinct() 

head(sample_tibble)

dGO2<-as.data.frame(sample_tibble[,2:4])
head(dGO2)
dim(dGO2)

# Write out the formatted table to not repeat everytime
#write.table(dGO2,"Dcor_v4_GO.txt", sep="\t", row.names = F, quote=FALSE)

########################################################

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("AnnotationForge")
#BiocManager::install("clusterProfiler")
#BiocManager::install("AnnotationHub")

library(AnnotationForge)

#import the gene symbol table
dSym<-read.table("Dcor_v4_SYM.txt",
             check.names=FALSE, header=T, 
             na.strings=c(""), stringsAsFactors=FALSE, sep="\t",quote="")
head(dSym)
dim(dSym)

dGO2<-read.table("Dcor_v4_GO.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(dGO2)
dim(dGO2)

## Then call the function
makeOrgPackage(gene_info=dSym, go=dGO2,
               version="0.5",
               maintainer="Sheila <so@someplace.org>",
               author="Sheila <so@someplace.org>",
               outputDir = ".",
               tax_id="866043",
               genus="Dalotia",
               species="coriaria5",
               goTable="go",
               verbose=TRUE)


## then you can call install.packages based on the return value
install.packages("./org.Dcoriaria5.eg.db", repos=NULL,type="source")


