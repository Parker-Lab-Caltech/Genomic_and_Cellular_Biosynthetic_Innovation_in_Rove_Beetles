library(forcats)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(cowplot)

gi<-read.table("./02_Other_Rove_Beetles_Assembly_Annotation/genomeInfo.txt", header=T,stringsAsFactors=TRUE, sep = "\t",quote = "")

gi
dim(gi)

# genome size plot
gi_gs<-gi[c(12:34,36:38),c(1,13:17)]
gi_gsm<-melt(gi_gs)

gi_gsm$species<-factor(gi_gsm$species, levels=rev(c("Tachinus","Sepedophilus","Coproporus",
                                          "Gymnusa","Deinopsis","Adinopsis","Cypha","Holobus","Oligota","Aleochara 1",
                                          "Aleochara 2","Aleochara 3","Leptusa","Oxypoda","Liometoxenus","Myllaena","Falagria","Lissagria",
                                          "Drusilla", "Earota", "Geostiba", "Atheta", "Dalotia 1", 
                                          "Ecitophya", "Ecitomorpha", "Ecitodaemon")))

col_kmer<-rev(c("#2D5276","#7fadd1","#9EC790","#FCBF49","#B71126"))

mean_gs<-gi_gsm %>% group_by(species) %>%
  summarize(mean = mean(value/1000000, na.rm = TRUE))

ggplot() + 
  geom_point(as.data.frame(gi_gsm), mapping=aes(fill=variable, y=value/1000000, x=species),size=3,shape=21, alpha=0.5) + 
  coord_flip()+theme_classic()+xlab("")+ylab("Genome Size (Mb)") +
  scale_fill_manual(values=col_kmer) +
  scale_y_continuous(expand = c(0.01, 0.1), labels = scales::comma,limits=c(0,3000), breaks = seq(0, 3000, by = 500))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(face = "italic")) +
  geom_point(data=mean_gs, mapping=aes(x=species, y=mean),shape=4,color="black",size=3,stroke=1.5)

#BUSCO
mel<-gi %>%
  select(species,BUSCO_S, BUSCO_D, BUSCO_F, BUSCO_M) %>%
  melt(value.name="value")

#par(mar=c(4, 4, 4, 4), mfrow=c(1,1),oma = c(1, 12, 1, 0.5),cex.lab=1.5)

mel<-filter(mel, species!= "Dalotia 2")

col_class<-rev(c("#3370A3","#B2ddeb","#FACA66","#B71126"))

mel$species<-factor(mel$species, levels=rev(c("Agrilus","Photinus","Tribolium","Aethina","Dendroctonus","Leptinotarsa","Anoplophora",
                                                "Onthophagus",
                                                "Nicrophorus","Ocypus","Philonthus","Tachinus","Sepedophilus","Coproporus",
                                                "Gymnusa","Deinopsis","Adinopsis","Cypha","Holobus","Oligota","Aleochara 1",
                                                "Aleochara 2","Aleochara 3","Leptusa","Oxypoda","Liometoxenus","Myllaena","Falagria","Lissagria",
                                                "Drusilla", "Earota", "Geostiba", "Atheta", "Dalotia 1", 
                                                "Ecitophya", "Ecitomorpha", "Ecitodaemon")))

mel$variable <- factor(mel$variable, levels=c("BUSCO_M","BUSCO_F","BUSCO_D","BUSCO_S"))

ggplot(as.data.frame(mel), aes(fill=variable, y=value, x=species)) + 
        geom_bar(position="stack", stat="identity") + 
        coord_flip()+theme_classic()+xlab("")+ylab("")+
        scale_fill_manual(values=col_class,labels = rev(c("Complete & single copy", "Complete & duplicated",
                                                      "Fragmented", "Missing")))+
        scale_y_continuous(limits = c(0,101), expand = expansion(mult = c(0, .1)))+
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.text.y = element_text(face = "italic"))
