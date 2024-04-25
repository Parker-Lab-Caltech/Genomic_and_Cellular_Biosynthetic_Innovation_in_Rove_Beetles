library(tidyr)
library(dplyr)
tab<-read.table("Eda_inact_L_UL_2.txt", sep ="\t", header=F, fill=T)
head(tab)

tab2<-separate(tab, col=V3, into=c('category', 'value'), sep=' ')

tab3<-tab2 %>% filter(is.na(V4))

tab3b<-spread(tab3, key = category, value = value) %>% select(-V4,-V5,-V6,-V7,-V8)

tab4<-tab2 %>% filter(!is.na(V4))

tab4b<-spread(tab4, key = V5, value = V8) %>% select(V1,V2,category, V6, COMPENSATION, 	`Deleted exon`, FS_DEL, FS_INS, `Missing exon`, SSM, START_MISSING, STOP) %>% rename(exon_num=category) %>% group_by(V1,V2) %>% mutate(exon_num2 = paste0(`exon_num`, collapse = ","), COMPENSATION2 = paste0(`COMPENSATION`, collapse = ","), FS_DEL2 = paste0(`FS_DEL`, collapse = ","),  FS_INS2 = paste0(`FS_INS`, collapse = ","), SSM2 = paste0(`SSM`, collapse = ","), START_MISSING2 = paste0(`START_MISSING`, collapse = ","), STOP2 = paste0(`STOP`, collapse = ","), V62 = paste0(`V6`, collapse = ","), del_exon = paste0(`Deleted exon`, collapse = ","), miss_exon = paste0(`Missing exon`, collapse = ",")) %>% select(-`Deleted exon`, -`Missing exon`) %>% select(V1,V2,exon_num2,COMPENSATION2, FS_DEL2, FS_INS2, SSM2, START_MISSING2, STOP2, V62, del_exon, miss_exon) %>% distinct()

head(as.data.frame(tab4b), n=40)

tab6<-read.table("./dcor_Ecitodaemon_TOGA_Qhox/temp/transcript_quality.tsv", sep ="\t", header=T) %>% separate(., col=Projection_ID, into=c('V1', 'V2'), sep="\\.(?=[^.]+$)")

tab6$V2<-as.integer(tab6$V2)

tab7<-read.table("./dcor_Ecitodaemon_TOGA_Qhox/temp/orthology_scores.tsv", sep ="\t", header=T) %>% rename(V1=gene, V2=chain)

tab7$V2<-as.integer(tab7$V2)

tab8 <- read.table("eda_loss_ul.txt",sep ="\t", header=F)

tab9<- read.table("./dcor_Ecitodaemon_TOGA_Qhox/loss_summ_data.tsv", sep ="\t", header=F) %>% select(V2,V3) %>% separate(., col=V2, into=c('V1', 'V2'), sep="\\.(?=[^.]+$)")

tab9$V2<-as.integer(tab9$V2)

tab5<- tab8 %>% left_join(tab9, by=c("V1"="V1")) %>% left_join(tab4b, by=c("V1"="V1", "V2"="V2")) %>% left_join(tab3b, by=c("V1"="V1", "V2"="V2")) %>% left_join (tab6, by=c("V1"="V1", "V2"="V2")) %>% left_join (tab7, by=c("V1"="V1", "V2"="V2"))

write.table(tab5, "mutations_ecitodaemon.txt", sep="\t", row.names=F)