##############################
## ancestral reconstruction ##
##############################
library(phytools)
library(treeio)
library(dplyr)

phy <- read.tree("FigTree_run3_newickFormat.nwk")

phy <- tree_subset(phy, "Nves", levels_back = 1)

par(mfrow=c(1,1))

plotTree(ladderize(phy),ftype="i")

phy$tip.label

phy$tip.label <- c("Adinopsis","Deinopsis","Gymnusa","Aleochara sp3","Aleochara sp1","Aleochara sp2", 
                   "Atheta","Dalotia","Ecitophya","Ecitomorpha","Ecitodaemon","Earota","Geostiba",
                   "Drusilla","Falagria","Lissagria", "Myllaena","Leptusa","Liometoxenus","Oxypoda",
                   "Cypha","Holobus", "Oligota","Coproporus","Sepediphilus","Tachinus","Nicrophorus")

phy$edge.length<-phy$edge.length*100

chem<- read.table("PA_chemical_compounds_022123.txt",
                  sep="\t", header=T, row.names=2, na.strings="-9")

head(chem)


# loop through tree with different chemical classes
pdf("ancestralReconstruction_022123.pdf",width=7, height=7)
par(mfrow=c(3,3),lwd = 0.2)
y  <- NULL
for(i in 2:9){
  tryCatch({
    print(colnames(chem)[i])
    svl<-to.matrix(chem[,i],levels(factor(chem[,i])))
    #svl<-setNames(chem[,i],chem$Genus)
    #svl<-to.matrix(chem[,i],levels(factor(svl)))
    #svl<-as.matrix(as.numeric(chem[,i]))
    rownames(svl)<-chem$Genus
    svl
    
    fitER<-ace(svl[,2],phy,model="ER",type="discrete", CI = TRUE)
    print(fitER)
    
    # plot tree with ancestral state reconstruction
    #cols<-setNames(palette()[1:length(unique(svl))],sort(unique(svl)))
    cols<-c(`0`="black",`1`="orange")
    plotTree(ladderize(phy),cols,type="phylogram",fsize=0.8,ftype="i",mar=c(2,2,2,2), lwd=0.8, offset=0.2)
    title(main=colnames(chem)[i])
    nodelabels(node=1:phy$Nnode+Ntip(phy),
               pie=fitER$lik.anc,piecol=cols,cex=1)
    add.simmap.legend(colors=cols,prompt=FALSE,x=7,y=5,fsize=0.2,vertical=TRUE,shape="circle")
    axisPhylo(side=1, cex=0.3, pos = 0.5)
    #tiplabels(pie=to.matrix(svl,sort(unique(svl))),piecol=cols,cex=0.7, pch=19, lwd=0.2)
    
    Pr<-svl
    Pr["Earota",]<-c(0.5,0.5)
    Pr["Aleochara sp1",]<-c(0.5,0.5)
    Pr["Ecitodaemon",]<-c(0.9,0.1)
    Pr["Ecitomorpha",]<-c(0.9,0.1)
    Pr["Falagria",]<-c(0.5,0.5)
    QQ<-rerootingMethod(phy,Pr)
    QQ
    
    tmp <- data.frame(QQ$marginal.anc) %>% mutate(run=print(colnames(chem)[i]))
    y <- rbind(y, tmp)
    
    tiplabels(pie=QQ$marginal.anc[phy$tip.label,], piecol=cols,cex=0.6)
    nodelabels(pie=QQ$marginal[as.character(1:phy$Nnode +length(phy$tip)),],piecol=cols,cex=1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
dev.off()

write.table(y,"ace_ER_probabilites.txt", sep="\t", quote=FALSE)
