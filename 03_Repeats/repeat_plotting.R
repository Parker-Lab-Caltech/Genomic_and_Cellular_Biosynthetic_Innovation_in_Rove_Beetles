library(forcats)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(cowplot)

gi<-read.table("./02_Other_Rove_Beetles_Assembly_Annotation/genomeInfo.txt", header=T,stringsAsFactors=TRUE, sep = "\t",quote = "")

gi
dim(gi)

#repeat content proportion (dnaPipTE results)
gi_rep<-gi %>%
  filter(species!="Ocypus", species !="Philonthus", species !="Tribolium")


col_beetles<-c(Agrilus="grey",Photinus="grey",Aethina="grey",Dendroctonus="grey",Leptinotarsa="grey",Anoplophora="grey",
               Onthophagus="grey",
               Nicrophorus="#D75F49",Tachinus="#D75F49",Sepedophilus="#D75F49",Coproporus="#D75F49",
               Gymnusa="black",Deinopsis="black",Adinopsis="black",Cypha="#067AA8",Holobus="#067AA8",Oligota="#067AA8",`Aleochara 1`="#067AA8",
               `Aleochara 2`="#067AA8",`Aleochara 3`="#067AA8",Leptusa="#067AA8",Oxypoda="#067AA8",Liometoxenus="#067AA8",Myllaena="#067AA8", Falagria="#067AA8",Lissagria="#067AA8",
               Drusilla="#067AA8", Earota="#067AA8", Geostiba="#067AA8", Atheta="#067AA8", `Dalotia 1`="#067AA8", `Dalotia 2`="#067AA8", 
               Ecitophya="#067AA8", Ecitomorpha="#067AA8", Ecitodaemon="#067AA8")

gi_rep$species<-factor(gi_rep$species, levels=c("Agrilus","Photinus","Aethina","Dendroctonus","Leptinotarsa","Anoplophora",
                                                "Onthophagus",
                                                "Nicrophorus","Tachinus","Sepedophilus","Coproporus",
                                                "Gymnusa","Deinopsis","Adinopsis","Cypha","Holobus","Oligota","Aleochara 1",
                                                "Aleochara 2","Aleochara 3","Leptusa","Oxypoda","Liometoxenus","Myllaena","Falagria","Lissagria",
                                                "Drusilla", "Earota", "Geostiba", "Atheta", "Dalotia 1", "Dalotia 2", 
                                                "Ecitophya", "Ecitomorpha", "Ecitodaemon"))

a<-ggplot(gi_rep, aes(fill=species, y=dnapipete_repeat, x=species)) + 
  geom_bar(stat="identity", col="black")+ coord_flip() +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by = 0.25),expand = c(0, 0), labels=c(0,0.25,0.5,0.75,1)) +
  scale_x_discrete(limits = rev(levels(gi_rep$species)),  expand = c(0, 0)) +
  scale_fill_manual(values=col_beetles)+theme_classic()+
  ylab("Proportion Repetitve")+
  theme(legend.position="none",
        axis.title.y=element_blank())

# repeat classes
gi2<-read.table("./03_Repeats/dnapipeTE_repeat_by_class.txt", header=T,stringsAsFactors=TRUE, sep = "\t",quote = "")

gi2

gi2$Class<-factor(gi2$Class, levels = rev(c("LTR",
                                            "LINE",
                                            "SINE",
                                            "DNA",
                                            "Helitron",
                                            "rRNA",
                                            "Low_Complexity/Satellites",
                                            "Simple_repeat",
                                            "Dc-Sat1",
                                            "Unknown",
                                            "non-repeat")))

gi2$Organism<-factor(gi2$Organism, levels=c("Agrilus","Photinus","Aethina","Dendroctonus","Leptinotarsa","Anoplophora",
                                            "Onthophagus",
                                            "Nicrophorus","Tachinus","Sepedophilus","Coproporus",
                                            "Gymnusa","Deinopsis","Adinopsis","Cypha","Holobus","Oligota","Aleochara 3",
                                            "Aleochara 2","Aleochara 1","Leptusa","Oxypoda","Liometoxenus","Myllaena","Falagria","Lissagria",
                                            "Drusilla", "Earota", "Geostiba", "Atheta", "Dalotia 1", "Dalotia 2", 
                                            "Ecitophya", "Ecitomorpha", "Ecitodaemon"))

col_class<-rev(c("#2D5276","#7fadd1","#B2ddeb","#d9edd4","#8FCCC3","#FFFdbc","#FED38D","#FCA85E","#DA5D42","#B71126", "grey90"))

b<-ggplot(gi2, aes(fill=Class, y=Per_dnaPipeTE_genome, x=Organism)) + 
  geom_bar(position="fill", stat="identity")+ coord_flip() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),expand = c(0, 0), labels=c(0,0.25,0.5,0.75,1))+
  scale_x_discrete(limits = rev(levels(gi2$Organism)),  expand = c(0, 0))+
  scale_fill_manual(values=col_class)+theme_classic()+
  ylab("Proportion of each Repeat Class ")+
  theme(axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.4, 'cm'),
        legend.position = "right")
b


# RepeatProfiler plots
gi3<-read.csv("./03_Repeats/Rep_Profiler_Dcor_sat1_4x_allSamples.csv", 
              header=T,stringsAsFactors=TRUE,quote = "",check.names = FALSE)

gi3_melt<-melt(gi3, id.vars = "Position")

gi3_melt$variable<-factor(gi3_melt$variable, levels=c("Agrilus","Photinus","Aethina","Dendroctonus","Leptinotarsa","Anoplophora",
                                                      "Onthophagus",
                                                      "Nicrophorus","Tachinus","Sepedophilus","Coproporus",
                                                      "Gymnusa","Deinopsis","Adinopsis","Cypha","Holobus","Oligota","Aleochara 3",
                                                      "Aleochara 2","Aleochara 1","Leptusa","Oxypoda","Liometoxenus","Myllaena","Falagria","Lissagria",
                                                      "Drusilla", "Earota", "Geostiba", "Atheta", "Dalotia 1", "Dalotia 2", 
                                                      "Ecitophya", "Ecitomorpha", "Ecitodaemon"))
mid<-max(gi3_melt$value)/2

y_min<-0

min_gi3<-list(gi3_melt %>% group_by(variable) %>% summarize(min(value)))
max_gi3<-list(gi3_melt %>% group_by(variable) %>% summarize(max(value)))

library(ggh4x)
c<-ggplot(gi3_melt, aes(y=value, x=Position)) + 
  geom_segment(aes(xend=Position, yend=0, color=value)) +
  geom_line()+ theme_classic()+
  scale_color_gradient2(low="#2D5276", 
                        mid="#FFFdbc",  
                        high="#B71126",
                        midpoint = mid) +
  facet_grid2(rows=vars(variable), scales="free_y",independent = "y", axes="y") +  
  facetted_pos_scales(x = min_gi3$`min(value)`, y =  max_gi3$`max(value)`)+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y=element_blank(),
        legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.4, 'cm'),
        legend.position = "none",
        panel.spacing=unit(0, "pt")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand =c(0,0), position = "right")+
  geom_vline(xintercept = c(147,294,441), linetype="dotted", size=1.2)
c

ggdraw() +
  draw_plot(b, 0, 0, 0.5, 1) +
  draw_plot(c, 0.5, 0, 0.5, 1)


# repeat landscape sat1
sat<-read.table("./03_Repeats/Dcor_sat_Kimura_estimates.divsum", 
                header=T,stringsAsFactors=TRUE,quote = "",sep=" ",check.names = FALSE)

sat<-sat[,c(-3,-4)]

sat$sat1<-sat$sat1/755000000
sat$species<-"Dalotia"

ggplot(sat, aes(x=Div, y=sat1, fill=species)) + 
  geom_bar(stat="identity", col="black")+
  scale_fill_manual(values=c("#DA5D42"))+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits=c(-0.5,20),expand = c(0, 0)) +
  theme_classic()+
  ylab("Proportion of Genomic Reads")+
  xlab("Kimura Substitution Level (%)")+
  theme(legend.position="none")


# repeat regions
region<-data.frame(region=c("exon","intron","intergenic"),
                   bp_coverage=c(6469, 20487,65801725))

region$region<- factor(region$region, levels=c("exon","intron","intergenic"))

ggplot(region,aes(x=region,y=log10(bp_coverage)))+
  geom_bar(stat="identity", fill="#DA5D42", col="black",width=0.7)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  ylab("bp coverage (log10)")+ xlab("")


ggplot(region,aes(x=region,y=bp_coverage))+
  geom_bar(stat="identity", fill="#DA5D42", col="black",width=0.7)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  ylab("bp coverage (log10)")+ xlab("")

# dcsat location along the chromosomes
library(karyoploteR)
library(stringr)
library(zoo)

chr<-read.table("./03_Repeats/Dcorv3.idx",
                check.names=FALSE, header=F, 
                na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")

head(chr)
colnames(chr)<-c("scaff","start","end", "chr")
chr2<-as.data.frame(cbind(chr$chr, chr$start, chr$end))
custom.genome <- toGRanges(chr2)
head(custom.genome)

bins <- tileGenome(seqinfo(custom.genome), tilewidth=50000,cut.last.tile.in.chrom=TRUE)

repgff<-read.table("./03_Repeats/Dcor_assembly_v3_220711.repeats.gff",
                   check.names=FALSE, header=F, 
                   na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(repgff)
repgff<-repgff[,c(-2,-3,-6,-7,-8)]
colnames(repgff)<-c("scaffold","start","end","TE", "Class")

repgff$scaffold<-str_replace(repgff$scaffold, pattern="hic_scaffold_", replacement="chr")
repgff$scaffold<-str_replace(repgff$scaffold, pattern="_pilon", replacement="")
rep2<-toGRanges(repgff %>% filter(scaffold %in% chr$chr) %>% mutate(bp=end-start))
head(rep2)


sat1<-repgff %>%
  filter(TE=="rnd-5_family-549")
s1<-toGRanges(sat1)
head(s1)

#repSat1<-read.table("./03_Repeats/Dcor_assembly_v3_220711_DcSat1.gff",
#                    check.names=FALSE, header=F, 
#                    na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
#head(repSat1)
#repSat1<-repSat1[,c(-2,-3,-6,-7,-8)] %>%
#  filter(V9=="Motif:Dcor_sat1_4x")
#colnames(repSat1)<-c("scaffold","start","end","TE")

#repSat1$scaffold<-str_replace(repSat1$scaffold, pattern="hic_scaffold_", replacement="chr")
#repSat1$scaffold<-str_replace(repSat1$scaffold, pattern="_pilon", replacement="")
#rep3<-toGRanges(repSat1)
#head(rep3)

LTR<-toGRanges(repgff %>% filter(Class =="LTR")%>% mutate(bp=end-start))
head(LTR)

LINE<-toGRanges(repgff %>% filter(Class =="LINE")%>% mutate(bp=end-start))
head(LINE)

SINE<-toGRanges(repgff %>% filter(Class =="SINE")%>% mutate(bp=end-start))
head(SINE)

DNA<-toGRanges(repgff %>% filter(Class =="DNA")%>% mutate(bp=end-start))
head(DNA)

Hel<-toGRanges(repgff %>% filter(Class =="RC")%>% mutate(bp=end-start))
head(Hel)

pp <- getDefaultPlotParams(plot.type=6)
pp$data1outmargin <- 5
pp$data1inmargin <- 2
pp$leftmargin<- 0.15
pp$topmargin <-50
pp$bottommargin <-50

newC<-col2rgb("#4B7CB7", alpha = 0.5)
gc.colors <- c("#4B7CB795", "#B7112670",  "#FCA85E", "#FFFDBC", "#B2DDEB")

kp <- plotKaryotype(genome=custom.genome,cex=0.8,plot.type=6, plot.params=pp,ideogram.plotter = NULL)
kpDataBackground(kp, color = "#FFFFFFFA")
kpAddBaseNumbers(kp, tick.dist = 1000000, tick.len = 1, tick.col="black", cex=0.5)

#kpDataBackground(kp, data.panel = 1, col="#AACBFF")
#kpAddLabels(kp, labels="% GC", col= "grey30",data.panel = 1, cex=0.6,label.margin = 0.05,srt=270,pos=1, side=2)
#kpAxis(kp, data.panel=1, cex=0.5,side = 2,col= "grey30", lwd=1.5, text.col="grey30", ymin=0.2, ymax=0.5)
#kpAddLabels(kp, labels="% repeat", data.panel =2, cex=0.6,label.margin = 0.03,srt=90,pos=1, side=1)
kpAxis(kp, data.panel=2, cex=0.5,side = 1)

kpAbline(kp, h=0.5, lty=2, col=gc.colors[1],data.panel="ideogram")

#repeats
kp<-kpPlotDensity(kp, data=rep2, col="grey80",window.size=50000,data.panel="ideogram", border="grey80")
wind<-kp$latest.plot$computed.values$windows
wind$reps<-kp$latest.plot$computed.values$density

kp<-kpPlotDensity(kp, data=DNA, col="#FCA85E80", r0=0, r1=1,window.size=50000,border="white")
wind$DNA<-(kp$latest.plot$computed.values$density)/wind$reps

kp<-kpPlotDensity(kp, data=Hel, col="#8FCCC380", r0=0.2, r1=1.2,window.size=50000,border="white")
wind$Hel<-(kp$latest.plot$computed.values$density)/wind$reps

kp<-kpPlotDensity(kp, data=LTR, col="#2D527680", r0=0.4, r1=1.4,window.size=50000,border="white")
wind$LTR<-(kp$latest.plot$computed.values$density)/wind$reps

kp<-kpPlotDensity(kp, data=LINE, col="#B7112680", r0=0.6, r1=1.6,window.size=50000,border="white")
wind$LINE<-(kp$latest.plot$computed.values$density)/wind$reps

kp<-kpPlotDensity(kp, data=SINE, col="#d9edd480", r0=0.8, r1=1.8,window.size=50000,border="white")
wind$SINE<-(kp$latest.plot$computed.values$density)/wind$reps

kp<-kpPlotDensity(kp, data=s1, col="#DA5D42", r0=1, r1=2,window.size=50000)
wind$DcSat1<-(kp$latest.plot$computed.values$density)/wind$reps

kp <- plotKaryotype(genome=custom.genome,cex=1.2,plot.type=6, plot.params=pp,ideogram.plotter = NULL)
kpDataBackground(kp, color = "grey99", r0=0, r1=1)
kpAddBaseNumbers(kp, tick.dist = 2000000, tick.len = 1, tick.col="black", cex=0.5)
#kpAxis(kp, data.panel=2, cex=0.5,side = 1)
kpPlotRegions(kp, data=LTR, col="#4B7CB755", border="#4B7CB755", r0=0.2, r1=0.4,num.layers=1,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=DNA, col="#FCA85E55", border="#FCA85E55", r0=0.4, r1=0.6,num.layers=1,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=Hel, col="#B2DDEB55", border="#B2DDEB55", r0=0.6, r1=0.8,num.layers=1,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=LINE, col="#B33B7755", border="#B33B7755",r0=0.8, r1=1,num.layers=1,avoid.overlapping=FALSE)
#kpPlotRegions(kp, data=SINE,  col="#4B7CB780", border="#4B7CB780",r0=0.6, r1=0.8)
#kpPlotRegions(kp, data=s1, col="#DA5D42", r0=0.8, r1=1)
kp<-kpPlotDensity(kp, data=s1, col="#DA5D42", border="grey30",r0=0, r1=1,window.size=50000)
#kpPlotRibbon(kp, data=wind, y0=wind$SINE, y1=wind$DcSat1, col="#DA5D42")

par(new=T)
legend(x = "topright", fill = c("#4B7CB790", "#FCA85E90","#B2DDEB90","#B33B7790","#DA5D42"), 
       legend = c("LTR", "DNA","Helitron","LINE","DcSat1"),xpd=TRUE)


#kpLines(kp, data=wind, y=wind$DNA, col="#FCA85E",lwd=1.5)
#kpLines(kp, data=wind, y=wind$Hel, col="#8FCCC3",lwd=1.5)
#kpLines(kp, data=wind, y=wind$LTR, col="#2D5276",lwd=1.5)
#kpLines(kp, data=wind,y=wind$LINE, col="#B71126",lwd=1.5)
#kpLines(kp, data=wind, y=wind$SINE, col="#4B7CB7",lwd=1.5)
#kpLines(kp, data=wind, y=wind$DcSat1, col="#DA5D42",lwd=2)


######
# repeat frequency heatmap
library(gplots)
library(RColorBrewer)
hm<-read.table("repeat_counts_2.txt",
                check.names=FALSE, header=T, 
                na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")

rownames(hm) <- hm$Family

# heatmap of the core vs random gene sets
heatmap.2(as.matrix(hm[,2:13]), scale="row", Colv = NA,margins = c(6,15),density.info = "none",trace = "none",
          dendrogram = "row",colsep=c(1,2,3,4,5,6,7,8,9,10,11),rowsep=1:nrow(hm),sepcolor = "grey",
          col = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(44),keysize=0.75, key.par = list(cex=0.5),
          breaks=seq(-2.2,2.2,0.1))

# heatmap of the genomewide estimates
heatmap.2(as.matrix(hm[,14:18]), scale="row", Colv = NA,margins = c(6,15),density.info = "none",trace = "none",
          dendrogram = "row",colsep=c(1,2,3,4,5),rowsep=1:nrow(hm),sepcolor = "grey",
          col =  colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(44),keysize=0.75, key.par = list(cex=0.5),
          breaks=seq(-2.2,2.2,0.1))


library(chisq.posthoc.test)
chisq.posthoc.test(hm[,c(17,18)],method = "bonferroni")

v<-chisq.test(t(hm[,c(5,6)]))
v

head(hm)
