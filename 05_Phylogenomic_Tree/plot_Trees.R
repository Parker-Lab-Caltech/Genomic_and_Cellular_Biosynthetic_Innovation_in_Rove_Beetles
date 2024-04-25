###########################
## Plot MCMCtree results ##
###########################
if (!any(rownames(installed.packages()) == "MCMCtreeR")) install.packages("MCMCtreeR")
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)

phy <- readMCMCtree("FigTree_run3_used_in_paper.tre", from.file = TRUE)

phy

par(mfrow=c(1,1))
par(mar = c(3, 1, 3, 3))
MCMC.tree.plot(phy, analysis.type = "MCMCtree", cex.tips = 1, time.correction = 100, scale.res = c("Period"), 
               plot.type = "phylogram", cex.age = 0.7, cex.labels = 0.6, edge.width = 2,node.ages = TRUE,
               relative.height = 0.03,  label.offset = 4, grey.bars = FALSE, show.tip.label = TRUE,
                no.margin = FALSE,col.age ="#00000080", add.time.scale=TRUE,add.abs.time = FALSE,
               lwd.bar = 4, node.method="bar")
axisPhylo(side=1, cex=0.3, pos = -1.5)

######################
## Plot ASTRAL Tree ##
######################
library("ggtree")

outgroups = c("PPYR", "Apla")

tree.name = "ASTRAL_phyloInfo_t2.tree"

astral.data <- read.astral(tree.name)

astral.data@phylo$tip.label <- c("Photinus","Agrilus","Tribolium","Aethina","Dendroctonus", "Leptinotarsa",
                                 "Anoplophora", "Onthophagus","Nicrophorus","Coproporus","Sepediphilus","Tachinus",
                                 "Gymnusa","Adinopsis","Deinopsis","Cypha","Oligota","Holobus",
                                 "Aleochara sp3","Aleochara sp1","Aleochara sp2", "Leptusa","Oxypoda","Liometoxenus",
                                 "Myllaena","Lissagria","Falagria","Drusilla", "Earota","Geostiba","Dalotia","Atheta",
                                 "Ecitodaemon","Ecitomorpha","Ecitophya")


quartets<-data.frame(astral.data@data[,9:11],node=astral.data@data$node)
#quartets$node <- 1:astral.data@phylo$Nnode+Ntip(astral.data@phylo)

pies <- nodepie(quartets, cols = 1:3)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = c("black","grey60","grey90")))

p<-ggtree(astral.data) + 
  geom_tiplab(fontface = 'italic',align=TRUE, linesize=.5) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+ ggplot2::xlim(0, 20)+
  scale_fill_manual(values = c(q1="black",q2="grey60",q3="grey90"))+
  geom_inset(pies, width = 0.1, height = 0.1) +
  theme(legend.position = "bottom") +
  geom_treescale(x=0, y=-0.5, width=2)
p