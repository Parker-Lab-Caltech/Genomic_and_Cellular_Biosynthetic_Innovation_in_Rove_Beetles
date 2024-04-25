install.packages("devtools")
devtools::install_github("schneebergerlab/findGSE")
install.packages("pracma")
install.packages("fGarch")
library("findGSE")

# Dalotia
findGSE(histo="./21mer_reads1.histo", sizek=21, outdir="./1/het_21mer",exp_hom=70)
# Coproporus
findGSE(histo="./21mer_reads2.histo", sizek=21, outdir="./2/het_21mer",exp_hom=40)
# Ecitophya
findGSE(histo="./21mer_reads3.histo", sizek=21, outdir="./3/het_21mer",exp_hom=45)
# Ecitomorpha
findGSE(histo="./21mer_reads11.histo", sizek=21, outdir="./11/het_21mer",exp_hom=75)
# Earota
findGSE(histo="./21mer_reads13.histo", sizek=21, outdir="./13/het_21mer",exp_hom=60)
# Deinopsis
findGSE(histo="./21mer_reads14.histo", sizek=21, outdir="./14/het_21mer",exp_hom=60)
# Ecitodaemon
findGSE(histo="./21mer_reads25_v2.histo", sizek=21, outdir="./25/het_21mer",exp_hom=50)#no fit
# Oxypoda
findGSE(histo="./21mer_reads33_v2.histo", sizek=21, outdir="./33/het_21mer",exp_hom=30)
# Drusilla
findGSE(histo="./21mer_reads35.histo", sizek=21, outdir="./35/het_21mer_v2",exp_hom=50)
# Geostiba
findGSE(histo="./21mer_reads38.histo", sizek=21, outdir="./38/het_21mer_v2",exp_hom=50)
# Liometoxenus
findGSE(histo="./21mer_reads39.histo", sizek=21, outdir="./39/het_21mer",exp_hom=80)
# Myllaena
findGSE(histo="./21mer_reads40.histo", sizek=21, outdir="./40/het_21mer_v2",exp_hom=100)
# Atheta
findGSE(histo="./21mer_reads41.histo", sizek=21, outdir="./41/het_21mer_v2",exp_hom=80)
# Leptusa
findGSE(histo="./21mer_reads42.histo", sizek=21, outdir="./42/het_21mer_v2",exp_hom=35)
# Falagria
findGSE(histo="./21mer_reads43.histo", sizek=21, outdir="./43/het_21mer_v2",exp_hom=120)
#Aleochara sp 2
findGSE(histo="./21mer_reads44.histo", sizek=21, outdir="./44/het_21mer")
# Lissagria
findGSE(histo="./21mer_reads45.histo", sizek=21, outdir="./45/het_21mer",exp_hom=40)
# Holobus
findGSE(histo="./21mer_reads46.histo", sizek=21, outdir="./46/het_21mer")
# Aleochara sp1
findGSE(histo="./21mer_reads47.histo", sizek=21, outdir="./47/het_21mer",exp_hom=140)
# Adinopsis
findGSE(histo="./21mer_reads48.histo", sizek=21, outdir="./48/het_21mer_2",exp_hom=200)
# Gymnusa
findGSE(histo="./21mer_reads52.histo", sizek=21, outdir="./52/het_21mer",exp_hom=140)
# Cypha
findGSE(histo="./21mer_reads55.histo", sizek=21, outdir="./55/het_21mer_2",exp_hom=200)
# Aleochara sp3
findGSE(histo="./21mer_reads56.histo", sizek=21, outdir="./56/het_21mer_4")
# Tachinus
findGSE(histo="./21mer_reads57.histo", sizek=21, outdir="./57/het_21mer")
# Sepediphilus
findGSE(histo="./21mer_reads58.histo", sizek=21, outdir="./58/het_21mer_2",exp_hom=200)
# Oligota
findGSE(histo="./21mer_reads59.histo", sizek=21, outdir="./59/het_21mer_4",exp_hom=140)