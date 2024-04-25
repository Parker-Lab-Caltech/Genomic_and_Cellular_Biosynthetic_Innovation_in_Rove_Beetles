########################
## chemical heat map ##
#######################
library(reshape2)
library(ggplot2)
library(patchwork)
library(grid)


chem<- read.table("PA_chemical_compounds_staphs_AB_TN_3.txt",
                  sep="\t", header=T)
head(chem)

chem$chemical_group<-factor(chem$chemical_group, 
                            levels=c("BQs","alkane_alkenes","alcohol","aldehydes",
                                     "terpenes","ketones",
                                     "aliphatic esters",
                                     "fatty acid",
                                     "aromatics","unidentified compounds"))

chem$chemical_group

chem2<-melt(chem,id.vars = c("chemical","chem_order","chemical_group"))
head(chem2)

chem2$variable<-factor(chem2$variable, levels=rev(unique(chem2$variable)))


chem2 <- chem2[order(as.numeric(chem2$chem_order),chem2$chemical ), ]
chem2$chemical <- factor(chem2$chemical , levels = unique(chem2$chemical))

ggplot(chem2, aes(x=variable, y=chemical, fill=as.factor(value)))+
  geom_tile(color=NA) +theme_bw() + coord_flip(expand = FALSE) +
  scale_fill_manual(values=c("#EAEAEA","#96E1FF","#067AA8"))+
  scale_y_discrete(position = "right",expand = c(0,0))+ 
  scale_x_discrete(expand = c(0,0))+theme_grey()+
  theme(plot.background = element_rect(fill = "transparent", colour = NA),legend.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0, size=8),
        plot.margin = unit(c(1,1,1.8,1), "cm"),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face="italic"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.8, "mm"),) +
  facet_grid(~chemical_group,scales = "free", space="free")
