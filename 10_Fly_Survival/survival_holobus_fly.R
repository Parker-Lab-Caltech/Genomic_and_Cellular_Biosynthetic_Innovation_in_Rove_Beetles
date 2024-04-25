fly<-read.table("Fly_larvae_survival.txt", sep="\t", header=T)

head(fly)

library(ggplot2)
library(gridExtra)

fly$treatment<-factor(fly$treatment, levels=rev(unique(fly$treatment)))

p1<-ggplot(fly, aes(x=treatment,y=survival_1h*10, group=treatment))+
  geom_boxplot(fill="grey",outlier.shape = NA) +theme_classic() + 
  geom_jitter(width = 0.2)+scale_y_continuous(expand=c(0,0))+
  scale_fill_grey()+
  ylab("Fly larval survival post-immersion (%)")

p2<-ggplot(fly, aes(x=treatment,y=eclosed*10, group=treatment))+
  geom_boxplot(fill="grey",outlier.shape = NA) +theme_classic() + 
  geom_jitter(width = 0.2)+scale_y_continuous(expand=c(0,0))+
  scale_fill_grey()+
  ylab("Adult fly survival post-immersion (%)")

grid.arrange(p1,p2, nrow = 1, ncol = 2)

library(multcompView)

mod<-aov(fly$survival_1h~fly$treatment)
summary(mod)
TH<-TukeyHSD(mod)

multcompLetters(TH$`fly$treatment`[,4])

mod<-aov(fly$eclosed~fly$treatment)
summary(mod)
TH<-TukeyHSD(mod)

multcompLetters(TH$`fly$treatment`[,4])

