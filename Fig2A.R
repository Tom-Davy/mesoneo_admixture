rm(list=ls())

library(data.table)
library(RColorBrewer)
library(plotrix)
options(scipen = 999)
source("~/EvoFunc/code/functionsSelSweep.R")
poschr<-fread("/Users/davyt/EvoFunc/local_ancestry/store/chrpos.allchr.txt")
#etwd("/Users/davyt/slurmp/")
#setwd("/Users/davyt/slurmp/newpanel2_HMM/")
#setwd("/Users/davyt/slurmp/newpanel2_HMM/whg")
#setwd("/Users/davyt/slurmp/newpanel2_HMM/ehg")   
setwd("/Users/davyt/slurmp/rebuttal/redo/hmm/")
HG.ancestry<-fread("HG.means.txt")

NEO.ancestry<-fread("NEO.means.txt")

HET.ancestry<-fread("HET.means.txt")

LAdf<-cbind(poschr,HG.ancestry,HET.ancestry,NEO.ancestry)
colnames(LAdf)<-c("V1","V2","HG.ancestry","HET.ancestry","NEO.ancestry")
setwd("~/EvoFunc/ms_scripts")
rm(HG,NEO,HET,HG.ancestry,NEO.ancestry,HET.ancestry)
LAdf$HGall<-LAdf$HG.ancestry+(LAdf$HET.ancestry/2)
HGsd=sd(LAdf$HGall)
HGnsd=sd(LAdf$HG.ancestry)
LAdf$HGallZ<-(LAdf$HGall-mean(LAdf$HGall))/HGsd
LAdf$HGZ<-(LAdf$HG.ancestry-mean(LAdf$HG.ancestry))/HGnsd

LAdf$HGZcomp<-(LAdf$HGallZ - LAdf$HGZ)

### Null approach for Z scores ####
LAdf$row.names=LAdf$V1
LAdf$marker=LAdf$V2
LAdf<-as.data.frame(LAdf)
nullwind=5000000
nullsnps<-windowsampler(nullwind=5000000, LAdf, 22, 1)

null=LAdf[nullsnps,]
nullsd=sd(null$HGall)
null$HGallZ<-(null$HGall-mean(null$HGall))/nullsd
#plot(null$HGallZ)

LAdf$nullZ<-(LAdf$HGall - mean(null$HGall))/nullsd
#plot(LAdf$nullZ)


df_FUT<-LAdf[which(LAdf$row.names==19 & LAdf$marker > 49199228 & LAdf$marker < 49209208),]
df_HERC2<-LAdf[which(LAdf$row.names==2 & LAdf$marker > 197058796 & LAdf$marker < 197458278),]
mean(df_HERC2$nullZ)
mean(df_HERC2$HGall)

plot(df_FUT$HGall, x=df_FUT$marker)
abline(v=49206674)
abline(h=mean(LAdf$HGall), lty=2)

#iain.df<-cbind(LAdf$V1, LAdf$V2, LAdf$HGall, LAdf$nullZ)
#colnames(iain.df)<-c("chr","pos","HGanc","nullZ")
#write.table(x=iain.df, file="Zscores_iain_newpanel3_nullZ.tsv", quote=F, sep = "\t")

delimit_MHC<-LAdf[which(LAdf$V1==6 & LAdf$V2 > 28477797 & LAdf$V2 < 33448354),]
delimit_classII<-delimit_MHC[delimit_MHC$V2 > 32000000,]

plot(delimit_classII$V2, y=delimit_classII$nullZ)


outlier_z<-3.5

#layout(matrix(c(1, 1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE), heights = c(1))
df_filt_markerprox<-as.data.frame(LAdf)
df_filt_markerprox$marker<-df_filt_markerprox$V2
df_filt_markerprox$row.names<-df_filt_markerprox$V1
SPACER=100000
for(CHR in 1:21){
  df_current<-which(df_filt_markerprox$row.names == CHR)
  df_succ<-which(df_filt_markerprox$row.names == CHR+1)
  current_position<-max(df_filt_markerprox[df_current,]$marker)
  print(current_position)
  df_filt_markerprox[df_succ,]$marker <- df_filt_markerprox[df_succ,]$marker + SPACER + current_position
}

df_filt_markerprox$markerprox<-df_filt_markerprox$marker
df_filt_markerprox$marker<-df_filt_markerprox$V2

LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
LAcol_grey2<-LAcol(9)
LAcol_grey2[4]<-LAcol_grey2[5]
LAcol_grey2[6]<-LAcol_grey2[5]

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist


HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3, plot_filter=1) 
#HG_outliers<-HG_outliers[which(HG_outliers$gene=="HLA-E"),]




NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3, plot_filter=1) 
NEO_outliers=NEO_outliers[c(5,10),]


tiff(filename="davyetal2022_fig2A.tiff", res = 150, units="in", width = 14, height = 7, compression="lzw") 
LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
#layout(matrix(c(1,1,1,1,2), nrow=1, byrow=TRUE), heights = c(1))


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$nullZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.6), type='p', axes=F, ylab="Mesolithic Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$nullZ))+5)]), pch=20, main = " Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox)))


df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$nullZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$nullZ))+5)]))

box()

axis(2)
chrindex=NULL
for (x in 1:22){
  chrindexT<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindexH<-head(df_filt_markerprox[(which   (df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindex[x]<-(chrindexH+chrindexT)/2
}
chrmidpoint<-chrindex


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1)


axis(3, at = NULL, labels = FALSE, tick = FALSE)
axis(4, at = NULL, labels = FALSE, tick = FALSE)




HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3.5, plot_filter=1) 
#HG_outliers<-HG_outliers[which(HG_outliers$gene=="HLA-E"),]


# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

topoutliers<-HG_outliers
topoutliers<-HG_outliers

alloutliers<-rbind(HG_outliers, NEO_outliers)
#write.table(alloutliers, file="localancestry_z3outliers.txt",sep="\t", quote = F)

topoutliers<-HG_outliers
topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #if(x == 4){
  #  text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #}else{
  # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
  
  text(x = xpoints, y = ypoints+0.015, cex=1, col = "black", labels = topoutliers_unique[x,]$gene) }
#points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)

#}


## NEO genes

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3, plot_filter=1) 
#NEO_outliers=NEO_outliers[c(5,10),]

topoutliers<-NEO_outliers


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #points(x = xpoints, y = ypoints-0.045, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints-0.015, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}


xmax=max(df_filt_markerprox$markerprox)
ymax=0.65





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.97*ymax,
  yt = ymax - 0.87*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=((xmax-0.9*xmax)+(xmax-0.5*xmax))/2, y=(ymax - 0.97*ymax + ymax - 0.87*ymax )/2, labels = "Z-score")

dev.off()
###### end plot #
###

tiff(filename="nogene_33rebut_redo_z3_spacer_bnohist_np4_localancestry_newpanel4_wholegenome_nullzHG3p5NEO3_grey2_rep_namefix.tiff", res = 150, units="in", width = 14, height = 7, compression="lzw") 
LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
#layout(matrix(c(1,1,1,1,2), nrow=1, byrow=TRUE), heights = c(1))


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$nullZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.6), type='p', axes=F, ylab="Mesolithic Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$nullZ))+5)]), pch=20, main = " Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox)))


df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$nullZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$nullZ))+5)]))

box()

axis(2)
chrindex=NULL
for (x in 1:22){
  chrindexT<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindexH<-head(df_filt_markerprox[(which   (df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindex[x]<-(chrindexH+chrindexT)/2
}
chrmidpoint<-chrindex


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1)


axis(3, at = NULL, labels = FALSE, tick = FALSE)
axis(4, at = NULL, labels = FALSE, tick = FALSE)




HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3, plot_filter=1) 
HG_outliers<-NULL

# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

topoutliers<-HG_outliers
topoutliers<-HG_outliers

alloutliers<-rbind(HG_outliers, NEO_outliers)
#write.table(alloutliers, file="localancestry_z3outliers.txt",sep="\t", quote = F)

topoutliers<-HG_outliers
topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #if(x == 4){
  #  text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #}else{
  # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
  if(topoutliers_unique[x,]$gene == "DNM1L"){
    text(x = xpoints, y = ypoints+0.03, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
  } else
    if(topoutliers_unique[x,]$gene == "HECW2"){
      text(x = xpoints-75000000, y = ypoints+0.025, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
    } else {
      text(x = xpoints, y = ypoints+0.015, cex=1, col = "black", labels = topoutliers_unique[x,]$gene) }
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}
#}


## NEO genes

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3, plot_filter=1) 
#NEO_outliers=NEO_outliers[c(5,10),]

topoutliers<-NULL


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #points(x = xpoints, y = ypoints-0.045, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints-0.015, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}


xmax=max(df_filt_markerprox$markerprox)
ymax=0.65





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.97*ymax,
  yt = ymax - 0.87*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=((xmax-0.9*xmax)+(xmax-0.5*xmax))/2, y=(ymax - 0.97*ymax + ymax - 0.87*ymax )/2, labels = "Z-score")

dev.off()



tiff(filename="pres_spacer_bnohist_np4_localancestry_newpanel4_wholegenome_nullzHG3p5NEO3_grey2_rep_namefix.tiff", res = 150, units="in", width = 14, height = 7, compression="lzw") 
LAcol<-colorRampPalette(colors = c('#e40959','red','grey','blue','#412ca8'))
#layout(matrix(c(1,1,1,1,2), nrow=1, byrow=TRUE), heights = c(1))

outlier_z<-3.5

#layout(matrix(c(1, 1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE), heights = c(1))
df_filt_markerprox<-as.data.frame(LAdf)
df_filt_markerprox$marker<-df_filt_markerprox$V2
df_filt_markerprox$row.names<-df_filt_markerprox$V1
SPACER=100000
for(CHR in 1:21){
  df_current<-which(df_filt_markerprox$row.names == CHR)
  df_succ<-which(df_filt_markerprox$row.names == CHR+1)
  current_position<-max(df_filt_markerprox[df_current,]$marker)
  print(current_position)
  df_filt_markerprox[df_succ,]$marker <- df_filt_markerprox[df_succ,]$marker + SPACER + current_position
}

df_filt_markerprox$markerprox<-df_filt_markerprox$marker
df_filt_markerprox$marker<-df_filt_markerprox$V2


LAcol_grey2<-LAcol(9)
LAcol_grey2[4]<-LAcol_grey2[5]
LAcol_grey2[6]<-LAcol_grey2[5]


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$nullZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.6), type='p', axes=F, ylab="Mesolithic Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$nullZ))+5)]), pch=20, main = " Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox)), cex.lab=1.5)


df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$nullZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$nullZ))+5)]))

box()

axis(2)
chrindex=NULL
for (x in 1:22){
  chrindexT<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindexH<-head(df_filt_markerprox[(which   (df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindex[x]<-(chrindexH+chrindexT)/2
}
chrmidpoint<-chrindex


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1.5, cex.lab=1.5)


axis(3, at = NULL, labels = FALSE, tick = FALSE)
axis(4, at = NULL, labels = FALSE, tick = FALSE)

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist



HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3.5, plot_filter=1) 
HG_outliers<-HG_outliers[which(HG_outliers$gene=="HLA-E"),]


# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

topoutliers<-HG_outliers
topoutliers<-HG_outliers

alloutliers<-rbind(HG_outliers, NEO_outliers)
#write.table(alloutliers, file="localancestry_z3outliers.txt",sep="\t", quote = F)

topoutliers<-HG_outliers
topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #if(x == 4){
  #  text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #}else{
  # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
  if(topoutliers_unique[x,]$gene == "DNM1L"){
    text(x = xpoints, y = ypoints+0.03, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
  } else
    if(topoutliers_unique[x,]$gene == "HECW2"){
      text(x = xpoints-75000000, y = ypoints+0.025, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
    } else {
      text(x = xpoints, y = ypoints+0.025, cex=1.5, col = "black", labels = topoutliers_unique[x,]$gene) }
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}
#}


## NEO genes

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = 3.4, plot_filter=1) 


topoutliers<-NEO_outliers


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  
  text(x = xpoints, y = ypoints-0.025, cex=1.5, col = "black", labels = topoutliers_unique[x,]$gene)
  
}


xmax=max(df_filt_markerprox$markerprox)
ymax=0.65





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.97*ymax,
  yt = ymax - 0.87*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)









text(x=((xmax-0.9*xmax)+(xmax-0.5*xmax))/2, y=(ymax - 0.97*ymax + ymax - 0.87*ymax )/2, labels = "Z-score")


dev.off()


############

tiff(filename="spacer_allgenesfix_nonames_bnohist_np4_localancestry_newpanel4_wholegenome_nullzHG3p5NEO3_grey2_rep_namefix.tiff", res = 150, units="in", width = 14, height = 7, compression="lzw") 
LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
#layout(matrix(c(1,1,1,1,2), nrow=1, byrow=TRUE), heights = c(1))

outlier_z<-3

#layout(matrix(c(1, 1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE), heights = c(1))
df_filt_markerprox<-as.data.frame(LAdf)
df_filt_markerprox$marker<-df_filt_markerprox$V2
df_filt_markerprox$row.names<-df_filt_markerprox$V1
SPACER=100000
for(CHR in 1:21){
  df_current<-which(df_filt_markerprox$row.names == CHR)
  df_succ<-which(df_filt_markerprox$row.names == CHR+1)
  current_position<-max(df_filt_markerprox[df_current,]$marker)
  print(current_position)
  df_filt_markerprox[df_succ,]$marker <- df_filt_markerprox[df_succ,]$marker + SPACER + current_position
}

df_filt_markerprox$markerprox<-df_filt_markerprox$marker
df_filt_markerprox$marker<-df_filt_markerprox$V2


LAcol_grey2<-LAcol(9)
LAcol_grey2[4]<-LAcol_grey2[5]
LAcol_grey2[6]<-LAcol_grey2[5]


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$nullZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.6), type='p', axes=F, ylab="Mesolithic Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$nullZ))+5)]), pch=20, main = " Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox)))


#plot(x=0.01, y=0.01, ylim=c(0.0,0.6), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$nullZ))+5)]), pch=20, main = " Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox)))




df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$nullZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$nullZ))+5)]))

box()

axis(2)
chrindex=NULL
for (x in 1:22){
  chrindexT<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindexH<-head(df_filt_markerprox[(which   (df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindex[x]<-(chrindexH+chrindexT)/2
}
chrmidpoint<-chrindex


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1)


axis(3, at = NULL, labels = FALSE, tick = FALSE)
axis(4, at = NULL, labels = FALSE, tick = FALSE)

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist



HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = outlier_z, plot_filter=1) 
#HG_outliers<-HG_outliers[which(HG_outliers$gene=="HLA-E"),]


# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

topoutliers<-HG_outliers

alloutliers<-rbind(HG_outliers, NEO_outliers)
#write.table(alloutliers, file="localancestry_z3outliers.txt",sep="\t", quote = F)

topoutliers<-HG_outliers
topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")
}
coli<-x
if(coli%%15 == 0){
  coli<-coli+1
}
# if(x == 4){
# text(x = xpoints-250000000, y = ypoints+0.03, cex=0.5, col = "black", labels = topoutliers_unique[x,]$gene)
# }else{
# points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
# if(topoutliers_unique[x,]$gene == "DNM1L"){
#   text(x = xpoints, y = ypoints+0.03, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
# } else
#   if(topoutliers_unique[x,]$gene == "HECW2"){
#     text(x = xpoints-75000000, y = ypoints+0.025, cex=1, col = "black", labels = topoutliers_unique[x,]$gene)
#   } else {
# text(x = xpoints, y = ypoints+0.015, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene) }
#points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)




## NEO genes

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$nullZ), bed = glist0, outlier_threshold = outlier_z, plot_filter=1) 


topoutliers<-NEO_outliers


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #points(x = xpoints, y = ypoints-0.045, pty=20, cex=0.7, pch = 16, col = "green")
  # text(x = xpoints, y = ypoints-0.015, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}


xmax=max(df_filt_markerprox$markerprox)
ymax=0.65





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.97*ymax,
  yt = ymax - 0.87*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=((xmax-0.9*xmax)+(xmax-0.5*xmax))/2, y=(ymax - 0.97*ymax + ymax - 0.87*ymax )/2, labels = "Z-score")

dev.off()
###### end plot #





# 
# #plot.window(ylim=c(0.0,0.65), xlim=c(0,1))
# #xhist <- hist(df_filt_markerprox$HGall, breaks = 20, plot = FALSE)
# #barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", ylim = c(0.0,0.65))
# 
# 
# dev.off()
# 
# # null Z ######
# 
# ########


# Homozygosity plot....



#tiff(filename="localancestry_newpanel3_wholegenome_z_grey2.tiff", res = 200, width = 1500, height = 1000)
LAcol<-colorRampPalette(colors = c('red','grey','royalblue4')) ########
layout(matrix(c(1,1,1,1,2), nrow=1, byrow=TRUE), heights = c(1))

#layout(matrix(c(1, 1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE), heights = c(1))
df_filt_markerprox<-as.data.frame(LAdf)
df_filt_markerprox$marker<-df_filt_markerprox$V2
df_filt_markerprox$row.names<-df_filt_markerprox$V1
SPACER=80000000
for(CHR in 1:21){
  df_current<-which(df_filt_markerprox$row.names == CHR)
  df_succ<-which(df_filt_markerprox$row.names == CHR+1)
  current_position<-max(df_filt_markerprox[df_current,]$marker)
  print(current_position)
  df_filt_markerprox[df_succ,]$marker <- df_filt_markerprox[df_succ,]$marker + SPACER + current_position
}

df_filt_markerprox$markerprox<-df_filt_markerprox$marker
df_filt_markerprox$marker<-df_filt_markerprox$V2


LAcol_grey2<-LAcol(9)
LAcol_grey2[4]<-LAcol_grey2[5]
LAcol_grey2[6]<-LAcol_grey2[5]


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$HGZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HG, ylim=c(0.0,0.6), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$HGZ))+5)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox))) 





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.95*ymax,
  yt = ymax - 0.85*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)



df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$HGZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HG, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$HGZ))+5)]))

box()

axis(2)
chrindex=NULL
for (x in 1:22){
  chrindexT<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindexH<-head(df_filt_markerprox[(which   (df_filt_markerprox$row.names==x)),]$markerprox, 1)
  chrindex[x]<-(chrindexH+chrindexT)/2
}
chrmidpoint<-chrindex


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1)


axis(3, at = NULL, labels = FALSE, tick = FALSE)
axis(4, at = NULL, labels = FALSE, tick = FALSE)

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist


HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = df_filt_markerprox$HGZ, bed = glist0, outlier_threshold = outlier_z, plot_filter=0) 
# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$HGZ), bed = glist0, outlier_threshold = outlier_z, plot_filter=1) 

topoutliers<-HG_outliers

alloutliers<-rbind(HG_outliers, NEO_outliers)
#write.table(alloutliers, file="localancestry_z3outliers.txt",sep="\t", quote = F)


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HG
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #if(x == 4){
  #  text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #}else{
  # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}
#}


## NEO genes


topoutliers<-NEO_outliers


topoutliers_store=NULL
for (x in 1:length(unique(topoutliers$gene))){
  GENE<-unique(topoutliers$gene)[x]
  store<-topoutliers[topoutliers$gene==GENE,]
  maxpval<-store[which(store$pval == max(store$pval)),]
  topoutliers_store[[x]]<-maxpval
}

topoutliers_unique<-do.call(rbind, topoutliers_store)


df_smol<-df_filt_markerprox



for (x in 1:nrow(topoutliers_unique)){
  print(nrow(topoutliers_unique)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  chrstore<-topoutliers_unique[x,]$chr
  markerstore<-topoutliers_unique[x,]$markerprox
  #mbpeaks<-df_smol[which(df_smol$row.names==chrstore  & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  xpoints<-topoutliers_unique[x,]$markerprox
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HG
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #points(x = xpoints, y = ypoints-0.045, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints-0.045, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}


xmax=max(df_filt_markerprox$markerprox)
ymax=0.65





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.97*ymax,
  yt = ymax - 0.87*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=(xmax-0.9*xmax)+(xmax-0.85*xmax), y=(ymax - 0.83*ymax), labels = "Z-score")


### adding histogram in next panel

#(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HG, ylim=c(0.0,0.65), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol(10)[(floor((df_filt_markerprox_pos$HGZ))+6)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox+(0.4*df_filt_markerprox$markerprox)))) 

par(mar=c(5, 1, 4, 2) + 0.1)


A<-hist(df_filt_markerprox$HG.ancestry, breaks = 40, plot = F)


#rect(0, A$breaks[1:(length(A$breaks) - 1)], A$counts, A$breaks[2:length(A$breaks)])



horiz.hist <- function(Data, breaks=60, col="transparent", las=1, ylim=range(HBreaks), labelat=pretty(ylim), labels=labelat, border=par("fg"), ... )
{a <- hist(Data, plot=FALSE, breaks=breaks)
HBreaks <- a$breaks
HBreak1 <- a$breaks[1]
hpos <<- function(Pos) (Pos-HBreak1)*(length(HBreaks)-1)/ diff(range(HBreaks))
barplot(a$counts, space=0, horiz=T, ylim=hpos(ylim), col=LAcol(3)[2], border="black",axes=F )      
#axis(2, at=hpos(labelat), labels=labels, las=las, ...) 
#print("use hpos() to address y-coordinates") 
}
horiz.hist(df_filt_markerprox$HG, ylim = c(0,0.6))



#plot.window(ylim=c(0.0,0.65), xlim=c(0,1))
#xhist <- hist(df_filt_markerprox$HG, breaks = 20, plot = FALSE)
#barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", ylim = c(0.0,0.65))


dev.off()


# output for iain part 1 ...LA
iain.df<-cbind(LAdf$V1, LAdf$V2, LAdf$HGall, LAdf$nullZ)
colnames(iain.df)<-c("chr","pos","HGanc","nullZ")
write.table(x=iain.df,file="Zscores_iain_newpanel3_nullZ.tsv", quote=F, sep = "\t", colnames=F)
# LA - CCDS analysis



library(data.table)
library(fitdistrplus)
poschr<-fread("/Users/davyt/EvoFunc/local_ancestry/store/chrpos.allchr.txt")
colnames(poschr)<-c("chrom","position")
#posteriors<-fread("/Users/davyt/EvoFunc/local_ancestry/all.sites.posterior")
darken <- function(color, factor=1.4){ # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

### Carp-----

colnames(LAdf)[c(1,2)]<-c("chrom","position")
LAdf$NEO<-LAdf$NEO.ancestry
LAdf$HG<-LAdf$HG.ancestry
LAdf$NEOall<-LAdf$NEO+0.5*LAdf$HET
LAdf$HGall<-LAdf$HG+0.5*LAdf$HET
LAdf$diffhet<-abs(LAdf$NEOall-LAdf$HET)
LAdf$diff<-abs(LAdf$NEO-LAdf$HG)
LAdf$hetdiff<-abs(LAdf$HET-(LAdf$NEO+LAdf$HG))

df<-LAdf
glist<-read.table("/Users/davyt/EvoFunc/glist.txt")
#colnames(df)<-c("row.names","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
df_filt<-df
df_filt$uid<-seq_along(df_filt[,1])

glist_sub<-glist[-which(glist$V1=="X" | glist$V1=="Y" | glist$V1 =="XY"),]
glist<-glist_sub
glist$V5<-make.unique(glist$V4)
glist$V6<-abs(glist$V2-glist$V3)
df_filt$marker<-df_filt$position
df_filt$row.names<-df_filt$chrom
# new method - indexing!

df_split<-split(df_filt, f = df_filt$chrom) # Split by Chr

n=100

df_splat<-lapply(df_split, function(x){
  nr=nrow(x)
  split(x, rep(1:ceiling(nr/n), each=n, length.out=nr))
})

df_splat_flat= unlist(df_splat, recursive = FALSE)

# 
# split.genomic<-function(df,distance=100000,marker="marker"){
#   matchmarker=which(match(colnames(x), marker) == 1)
#   x$interval<-findInterval(x$marker, seq(from=0, to=max(x[, ..matchmarker]), by=distance))
#   xx<-split(x, f = x$interval)
#   return(xx)
# }
# 

# xy<-split.genomic(df_split)
#   

# distance=100000

df_splot<-lapply(df_split, function(x){
  matchmarker=which(match(colnames(x), "marker") == 1)
  x$interval<-findInterval(x$marker, seq(from=0, to=max(x[,matchmarker]), by=100000))
  return(split(x, f = x$interval))
})

df_splot<-unlist(df_splot, recursive = FALSE)

splot_SNPs<-lapply(df_splot, function(x){
  return(nrow(x))
})

splot_means_NEOall<-lapply(df_splot, function(x){
  return(mean(x$NEOall))
})

plot(x=splot_means_NEOall, y=splot_SNPs)


abline(v=mean(df_filt$NEOall), col="grey", lty=5)
abline(v=mean(df_filt$NEOall)-(1*sd(df_filt$NEOall)), col="red", lty=5)
abline(v=mean(df_filt$NEOall)-(2*sd(df_filt$NEOall)), col="red1", lty=5)
abline(v=mean(df_filt$NEOall)-(3*sd(df_filt$NEOall)), col="red2", lty=5)
abline(v=mean(df_filt$NEOall)-(4*sd(df_filt$NEOall)), col="red3", lty=5)

abline(v=mean(df_filt$NEOall)+(1*sd(df_filt$NEOall)), col="royalblue", lty=5)
abline(v=mean(df_filt$NEOall)+(2*sd(df_filt$NEOall)), col="royalblue1", lty=5)
abline(v=mean(df_filt$NEOall)+(3*sd(df_filt$NEOall)), col="royalblue2", lty=5)
abline(v=mean(df_filt$NEOall)+(4*sd(df_filt$NEOall)), col="royalblue3", lty=5)



#

splat_index<-lapply(df_splat, function(x){
  lapply(x, function(y){
    index_a<-unlist(head(y, 1))[c(1,2)]
    index_b<-unlist(tail(y, 1))[2]
    print(index_a)
    return(c(index_a,index_b))
  })
})

ptm <- proc.time()
df_ccds<-apply(glist, 1, function(x){
  #print(x) # print glist row
  z<-t(as.data.frame(splat_index[[as.numeric(x[1])]])) # as.numeric(x[1]) -> Current chr;splat_index is a list, and each list item is an list of positions such as to index df_splat, which is a nested list itself.
  z.min<-max(which(z[,2] < as.numeric(x[2])))
  if(z.min == "-Inf"){
    z.min<-1
  }
  z.max<-min(which(z[,3] > as.numeric(x[3])))
  
  out<-df_splat[[as.numeric(x[1])]][[z.min]]
  if(z.max != z.min){
    out<-rbind(out,df_splat[[as.numeric(x[1])]][[z.max]] )} # We concatenate the relevant look-up segments to finally delineate our genic co-ordinates.
  out$GENE<-as.character(x[5])
  out[which(out$position > as.numeric(x[2]) & out$position < as.numeric(x[3])),]
})
proc.time() - ptm

df_ccds2<-lapply(df_ccds, function(x){
  x$NEO<-x$NEO.ancestry
  x$HG<-x$HG.ancestry
  x$HET<-x$HET.ancestry
  return(x)
})



df_ccds<-df_ccds2


names(df_ccds)<-glist$V4

gene_means_HET<-lapply(df_ccds, function(x){
  return(mean(x$HET))
})

gene_means_NEO<-lapply(df_ccds, function(x){
  return(mean(x$NEO))
})

gene_means_HG<-lapply(df_ccds, function(x){
  return(mean(x$HG))
})

gene_means_NEOall<-lapply(df_ccds, function(x){
  return(mean(x$NEOall))
})

gene_means_HET<-lapply(df_ccds, function(x){
  return(mean(x$HET))
})


gene_means_HGall<-lapply(df_ccds, function(x){
  return(mean(x$HGall))
})



gene_SNPs<-lapply(df_ccds, function(x){
  return(nrow(x))
})
# 
# 
# 
# legend(x="topleft", legend=c("Proximal to HL
LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
col.HG<-"#412ca8"
col.NEO<-"#e40959"
col.IBM<-c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000")


gene_SNPs_scaled<-as.numeric(gene_SNPs)/as.numeric(glist$V6)
loggene_SNPs_scaled<--log(gene_SNPs_scaled)
tiff(filename="v50_hues_panel3_nullLAvsGENESNPS_volcano_logscaled_NEO42_HG_HLAdemarked_publication_res150.tiff",  units="in", width=8, height=8, res=150, compression = "lzw")
plot(y=as.vector(loggene_SNPs_scaled), x=as.vector(gene_means_HGall), pch = 20, col = "grey", ylab = "Number of SNPs per gene (ln)", xlab = "Mesolithic Ancestry", main = "Mesolithic Ancestry per-gene in the Admixed Mid-Neolithic", xlim=c(0.05, 0.6))


# 
# abline(v=mean(df_filt$NEOall), col="grey10", lty=5)
# abline(v=mean(df_filt$NEOall)-(1*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)-(2*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)-(3*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)-(4*sd(df_filt$NEOall)), col="grey", lty=5)
# 
# abline(v=mean(df_filt$NEOall)+(1*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)+(2*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)+(3*sd(df_filt$NEOall)), col="grey", lty=5)
# abline(v=mean(df_filt$NEOall)+(4*sd(df_filt$NEOall)), col="grey", lty=5)
cols<-colorRampPalette(brewer.pal(9,"YlOrRd"))
ccds<-do.call(rbind, df_ccds)

abline(v=mean(ccds$HGall), col="grey10", lty=5)
abline(v=mean(ccds$HGall)+(1*sd(ccds$HGall)), col=add.alpha(col.HG, 0.5), lty=5)
abline(v=mean(ccds$HGall)+(2*sd(ccds$HGall)), col=add.alpha(col.HG, 0.5), lty=5)
abline(v=mean(ccds$HGall)+(3*sd(ccds$HGall)), col=add.alpha(col.HG, 0.5), lty=5)
abline(v=mean(ccds$HGall)+(4*sd(ccds$HGall)), col=add.alpha(col.HG, 0.5), lty=5)

abline(v=mean(ccds$HGall)-(1*sd(ccds$HGall)), col=add.alpha(col.NEO, 0.5), lty=5)
abline(v=mean(ccds$HGall)-(2*sd(ccds$HGall)), col=add.alpha(col.NEO, 0.5), lty=5)
abline(v=mean(ccds$HGall)-(3*sd(ccds$HGall)), col=add.alpha(col.NEO, 0.5), lty=5)
abline(v=mean(ccds$HGall)-(4*sd(ccds$HGall)), col=add.alpha(col.NEO, 0.5), lty=5)



# sapply(seq_along(MHCgenes),function(x){
#    print(x)
#    z<-which(glist$V4 == MHCgenes[x])
#    points(y=as.vector(gene_SNPs)[z], x=as.vector(gene_means_NEOall)[z], pch=12, cex=0.5, col="black")
#   points(y=as.vector(gene_SNPs)[z], x=as.vector(gene_means_NEOall)[z], col = colfunc(length(MHCgenes))[x], pch=20, cex=1)
# } )
# HLAgenes<-glist[grep(glist$V4, pattern = "^HLA"),4]
# sapply(HLAgenes,function(x){
#   z<-which(glist$V4 == x)
#   points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_NEOall)[z], col = "black", pch=15, cex=1)
# } )

proximal=300000
z<-which(glist$V4 =="HLA-E")
nearHLAE<-glist[which(glist$V1==6 & glist$V2 > glist[z,]$V2-proximal & glist$V3 <  glist[z,]$V3 + proximal),]$V4
s<-which(glist$V4 == "SLC24A5")
nearSLC24A5<-glist[which(glist$V1 == glist[s,]$V1 & glist$V2 > glist[s,]$V2-proximal & glist$V3 <  glist[s,]$V3 + proximal),]$V4
x<-which(glist$V4 == "IQCJ")
nearIQCJ<-glist[which(glist$V1 == glist[x,]$V1 & glist$V2 > glist[x,]$V2-proximal & glist$V3 <  glist[x,]$V3 + proximal),]$V4
c<-which(glist$V4 == "DNM1L")
nearDNM1L<-glist[which(glist$V1 == glist[c,]$V1 & glist$V2 > glist[c,]$V2-proximal & glist$V3 <  glist[c,]$V3 + proximal),]$V4
d<-which(glist$V4 == "YARS2")
nearYARS<-glist[which(glist$V1 == glist[d,]$V1 & glist$V2 > glist[,]$V2-proximal & glist$V3 <  glist[d,]$V3 + proximal),]$V4


sapply(nearDNM1L,function(x){
  z<-which(glist$V4 ==x)
  points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = col.IBM[4])
} )


# sapply(nearYARS,function(x){
#    z<-which(glist$V4 == x)
#    points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = "purple4")
#  } )




sapply(nearHLAE,function(x){
  z<-which(glist$V4 == x)
  points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = col.IBM[1])
} )


sapply(nearSLC24A5,function(x){
  z<-which(glist$V4 == x)
  points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = col.IBM[2])
} )


sapply(nearIQCJ,function(x){
  z<-which(glist$V4 == x)
  points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = col.IBM[3])
} )


sapply(nearYARS2,function(x){
  z<-which(glist$V4 == x)
  points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.2,  pch=20, col = col.IBM[5])
} )



z<-which(glist$V4 =="HLA-E")
points(y=as.vector(loggene_SNPs_scaled)[z], x=as.vector(gene_means_HGall)[z], cex=1.5,  pch=18, col = col.IBM[1])
text(y=as.numeric(as.vector(loggene_SNPs_scaled)[z]), x=as.numeric(as.vector(gene_means_HGall)[x])+0.035, cex=1.1,  labels = "HLA-E", col = col.IBM[1])

s<-which(glist$V4 == "SLC24A5")
points(y=as.vector(loggene_SNPs_scaled)[s], x=as.vector(gene_means_HGall)[s], cex=1.5,  pch=18, col = col.IBM[2])
text(y=as.numeric(as.vector(loggene_SNPs_scaled)[s]), x=as.numeric(as.vector(gene_means_HGall)[s])-0.055, cex=1,  labels = "SLC24A5", col = col.IBM[2])


x<-which(glist$V4 == "IQCJ")
points(y=as.vector(loggene_SNPs_scaled)[x], x=as.vector(gene_means_HGall)[x], cex=1.5,  pch=18, col = col.IBM[3])
text(y=as.numeric(as.vector(loggene_SNPs_scaled)[x]), x=as.numeric(as.vector(gene_means_HGall)[x])+0.035, cex=1,  labels = "IQCJ", col = col.IBM[3])



x<-which(glist$V4 == "DNM1L")
points(y=as.vector(loggene_SNPs_scaled)[x], x=as.vector(gene_means_HGall)[x], cex=1.5,  pch=18, col = col.IBM[4])
text(y=as.numeric(as.vector(loggene_SNPs_scaled)[x]), x=as.numeric(as.vector(gene_means_HGall)[x])+0.055, cex=1,  labels = "DNM1L", col = col.IBM[4])


x<-which(glist$V4 == "YARS2")
points(y=as.vector(loggene_SNPs_scaled)[x], x=as.vector(gene_means_HGall)[x], cex=1.5,  pch=18, col = col.IBM[5])
text(y=as.numeric(as.vector(loggene_SNPs_scaled)[x]), x=as.numeric(as.vector(gene_means_HGall)[x])+0.055, cex=1,  labels = "YARS2", col = col.IBM[5])


legend(x="topleft", legend=c("Proximal to HLA-E","Proximal to SLC24A5","Proximal to IQCJ","Proximal to DNM1L","Proximal To YARS2"),  col=c(col.IBM[1:5]), pch=c(20))


dev.off()




peak1<-delimit_MHC[which(delimit_MHC$marker > 30000000 & delimit_MHC$marker < 31000000),]
plot(peak1$nullZ, x=peak1$marker)
peak1[min(which(peak1$nullZ > 2)),]


peak2<-delimit_MHC[which(delimit_MHC$marker > 32000000 & delimit_MHC$marker < 32999330+100000),]
plot(peak2$nullZ, x=peak2$marker)
peak2[min(which(peak2$nullZ > 2)),]
peak2[max(which(peak2$nullZ > 2)),]


glist[glist$V4=="HERC2",]
df<-df_filt_markerprox
HERC<-df[which(df$row.names==15 & df$marker > 28356182 & df$marker < 28567298),]
