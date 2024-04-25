r1<-read.table("HAL_relax.txt", sep="\t", header=T)


r2<-r1 %>%
  group_by(Group) %>%
  mutate(p.adj= p.adjust(PVALUE, method = "fdr"))

library(dplyr)
library(ggplot2)
library(ggridges)

r3<-r2 %>% mutate(direction= ifelse(p.adj < 0.1, direction, "NA"))

r3$Group<-as.factor(r3$Group)
r3$direction<-as.factor(r3$direction)

ggplot(r3 %>% filter(direction != "NA") %>% filter(p.adj <= 0.1) , aes(x=log2_K, y=Group, fill=Group,group=Group)) + 
  #geom_vline(aes(xintercept=-1.56), color="blue",linetype="dashed", size=1) +
  #geom_vline(aes(xintercept=-2.36), color="red",linetype="dashed", size=1) +
  geom_density_ridges(alpha=0.5, color="grey30", scale = 7, size=0.6) +
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1) +
  scale_fill_manual(values=c("#20F1F5","#B326E2"))+ xlab(paste("selection intensity K (log2)"))+ ylab("")+
  theme_ridges() + scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0.01, 0))

r3 %>% filter(direction != "NA") %>% filter(p.adj <= 0.1) %>%group_by(Group) %>% summarize(mean(log2_K))

# Need other groups
ggplot(r3 %>% filter(direction != "NA") %>% filter(p.adj <= 0.1), aes(x=Group, y=log2_K, fill=Group)) + 
  geom_violin() +
  theme_classic()

