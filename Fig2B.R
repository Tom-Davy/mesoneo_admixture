# zoom-in of AGDP MW, 1240k Mw p-val & LA.
rm(list=ls())
setwd("~/EvoFunc/ms_scripts")
## AGDP -
source("~/EvoFunc/code/functionsSelSweep.R")
library(stringr)
library(plotrix)
library(data.table)
library(scales)

mapp<-fread("/Users/davyt/Downloads/wetransfer_heng-li-mapp-tracks_2022-10-26_1431/hs37m_filt35_99.bed")
colnames(mapp)<-c("chr","marker1","marker2")


####

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist
glist2<-glist[grep(glist$V4, pattern="HLA"),]
glist<-glist2
glist3<-glist
glist15<-glist0[glist0$V1==15,]

# gwasplot<-fread("/Users/davyt/slurmp/allgwas_standardized.6.tsv")
# colnames(gwasplot)<-c("chrom","pos","alt","ref","test","low_confidence","beta","pval","study","trait")




df<-as.data.frame(fread("/Users/davyt/slurmp/rebuttal/redo/REF_MNEO_HG_NEO_freqstats_mincount1_intersect.tsv"))
fread("~/slurmp/iain_sent/final/")
colnames(df)<-c("row.names","marker1","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
#df<-df[,-11]
#df<-df[-which(df$row.names==23),]
#df<-df[-which(df$row.names==24),]
#df$rawF4<-df$datawg
df$df.uid<-seq_along(df[,1])



df$expected_p2<-(df$p3*0.3)+(df$p4*0.7)
y<-df$expected_p2
df$datawg<-((df$p2-y)^2)/(2*(1-y^2-(1-y)^2))




df$datawg_nonsliding<-df$datawg
df<-na.omit(df)
df<-df[-which(df$datawg == "Inf"),]
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0

filename=paste0("freqtest_fix4", snp,"mincount20_5kyafix_nullfix",paste0("NEO",NEOFILT,"NOTNEOFILT",NOTNEOFILT))









#### GT lengths -----
df$GT1L<-str_length(df$GT1)/2
df$GT2L<-str_length(df$GT2)/2
df$GT3L<-str_length(df$GT3)/2
df$GT4L<-str_length(df$GT4)/2

# Mincount Hist
GT1<-hist(str_length(df$GT1)/2, breaks = 100, plot = FALSE)
GT2<-hist(str_length(df$GT2)/2, breaks = 100, plot = FALSE)
GT3<-hist(str_length(df$GT3)/2, breaks = 100, plot = FALSE)
GT4<-hist(str_length(df$GT4)/2, breaks = 50, plot = FALSE)



# Mincount Hist

plot(GT1)
plot(GT2)
plot(GT3)
plot(GT4)
# top fads e.g. rs174538
df<-df[which(df$GT4L >= 20),]

#df<-df[,-11]
#df<-df[-which(df$row.names==23),]
#df<-df[-which(df$row.names==24),]
#df$rawF4<-df$datawg
df$df.uid<-seq_along(df[,1])
#df<-df[-which(is.na(df$p2)),]
#df$expected_p2<-(df$p3*0.3)+(df$p4*0.7)

df$F2_NEOREF<-(df$p4-df$p1)^2
df$F2_MNEOREF<-(df$p2-df$p1)^2
df$F2_HGREF<-(df$p3-df$p1)^2


df<-na.omit(df)
#df<-df[-which(df$datawg == "Inf"),]
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0

filename=paste0("freqtest_fix4", snp,"mincount20_5kyafix_nullfix",paste0("NEO",NEOFILT,"NOTNEOFILT",NOTNEOFILT))

#### GT lengths -----
df$GT1L<-str_length(df$GT1)/2
df$GT2L<-str_length(df$GT2)/2
df$GT3L<-str_length(df$GT3)/2
df$GT4L<-str_length(df$GT4)/2

# Mincount Hist
GT1<-hist(str_length(df$GT1)/2, breaks = 100, plot = FALSE)
GT2<-hist(str_length(df$GT2)/2, breaks = 100, plot = FALSE)
GT3<-hist(str_length(df$GT3)/2, breaks = 100, plot = FALSE)
GT4<-hist(str_length(df$GT4)/2, breaks = 50, plot = FALSE)
# top fads e.g. rs174538

rs174538.GT1<-str_length(df[which(df$rs.number=="rs174538"),]$GT1)
rs174538.GT2<-str_length(df[which(df$rs.number=="rs174538"),]$GT2)
rs174538.GT3<-str_length(df[which(df$rs.number=="rs174538"),]$GT3)
rs174538.GT4<-str_length(df[which(df$rs.number=="rs174538"),]$GT4)


rs174546.GT1<-str_length(df[which(df$rs.number=="rs174546"),]$GT1)
rs174546.GT2<-str_length(df[which(df$rs.number=="rs174546"),]$GT2)
rs174546.GT3<-str_length(df[which(df$rs.number=="rs174546"),]$GT3)
rs174546.GT4<-str_length(df[which(df$rs.number=="rs174546"),]$GT4)

#Single Plots
plot(GT1, col = "gray")
abline(v=rs174538.GT1, lty = 4, col = "red", lwd =2)
abline(v=rs174546.GT1, lty = 4, col = "blue", lwd =2)

plot(GT2, col = "gray")
abline(v=rs174538.GT2, lty = 4, col = "red", lwd =2)
abline(v=rs174546.GT2, lty = 4, col = "blue", lwd =2)

plot(GT3, col = "gray")
abline(v=rs174538.GT3, lty = 4, col = "red", lwd =2)
abline(v=rs174546.GT3, lty = 4, col = "blue", lwd =2)

plot(GT4, col = "gray")
abline(v=rs174538.GT4, lty = 4, col = "red", lwd =2)
abline(v=rs174546.GT4, lty = 4, col = "blue", lwd =2)

#Multi-plot
plot(GT1, col  = "green")
plot(GT2, add =TRUE, col = "purple")
plot(GT3, add =TRUE, col = "blue")
plot(GT4, add = TRUE, col = "red")
abline(v=15, col ="black", lwd =2)
abline(v=40, col="black", lwd = 2)
legend(x="topleft", legend = c("MSL","MNEO","HG","NEO"), lty = 1, lwd = 3, col = c("green","purple","blue","red"))



glist<-read.table("/Users/davyt/EvoFunc/glist.txt")

# Nei's Gene Diversity Test

df$NGDMSL<-1-( ( df$p1^2) + ((1-df$p1)^2))
df$NGDMNEO<-1-( ( df$p2^2) + ((1-df$p2)^2))
df$NGDHG<-1-( ( df$p3^2) + ((1-df$p3)^2))
df$NGDNEO<-1-( ( df$p4^2) + ((1-df$p4)^2))



df$NGDMSL2<-1-df$p1^2
df$NGDMNEO2<-1-df$p2^2
df$NGDHG2<-1-df$p3^2
df$NGDNEO2<-1-df$p4^2


df$HZdiff<-abs(df$NGDNEO - df$NGDMNEO)
df$HZdiff2<-abs(df$NGDHG - df$NGDMNEO)

##MEDIAN CORRECTION
MSLmed = median(df$GT1L)
MNEOmed = median(df$GT2L)
HGmed = median(df$GT3L)
NEOmed = median(df$GT4L)

##

df$NGDMSLcorr<-1-( ( df$p1^2) + ((1-df$p1)^2)) * ((df$GT1L/ (df$GT1L - 1)))
df$NGDMNEOcorr<-1-( ( df$p2^2) + ((1-df$p2)^2))* ((df$GT2L / (df$GT2L - 1)))
df$NGDHGcorr<-1-( ( df$p3^2) + ((1-df$p3)^2))* ((df$GT3L / (df$GT3L - 1)))
df$NGDNEOcorr<-1-( ( df$p4^2) + ((1-df$p4)^2))* ((df$GT4L  / (df$GT4L  - 1)))

## MEANCORRECTION


MSLmean = mean(df$GT1L)
MNEOmean = mean(df$GT2L)
HGmean = mean(df$GT3L)
NEOmean = mean(df$GT4L)


##
df$NGDMSLcorrmean<-1-( ( df$p1^2) + ((1-df$p1)^2)) * ((MSLmean / (MSLmean - 1)))
df$NGDMNEOcorrmean<-1-( ( df$p2^2) + ((1-df$p2)^2))* ((MNEOmean / (MNEOmean - 1)))
df$NGDHGcorrmean<-1-( ( df$p3^2) + ((1-df$p3)^2))* ((HGmean / (HGmean - 1)))
df$NGDNEOcorrmean<-1-( ( df$p4^2) + ((1-df$p4)^2))* ((NEOmean / (NEOmean - 1)))



#df_filt<-df[which(df$datawg != 0),]
df_filt<-df
datawg<-df_filt$datawg


## filtering based on coverage -----
datawg<-datawg
df_filt_temp<-df_filt[which(df_filt$GT4L > NEOFILT & df_filt$GT1 > 0 & df_filt$GT2 > NOTNEOFILT & df_filt$GT3 > NOTNEOFILT),]
df_filt<-df_filt_temp
#df_filt<-df_filt[which(df_filt$datawg>0),]
##inverse
#datawg<-abs(df_filt$datawg)
datawg<-df_filt$datawg
df_filt$datawg<-datawg
##
df_filt$uid<-seq_along(df_filt[,1])
##

df_prime<-df_filt
#datawg<-datawg+abs(min(datawg))+0.00001
df_filt$datawg<-datawg
plot(density(df_filt$datawg))
row.names<-22
nullwind=5000000
df_mixedweight_keep<-df_filt

#### window-distance null ------ -----
nullwind=5000000
ptm <- proc.time()
nullsnps<-windowsampler(nullwind, df_filt, row.names, snp)
proc.time() - ptm
nullsnps_dist<-df_filt$marker[nullsnps[snp,]]-df_filt$marker[nullsnps[1,]]
if(snp>1){
  datanull<-get_null(datawg, snp, nullsnps)
  #datanull<-datanull/nullsnps_dist
}
# Null Classic -----
row.names<-22
nullsnps<-windowsampler(nullwind, df_filt, row.names, snp)

datanull<-get_null(datawg, snp, nullsnps)
#datanull<-datanull[-which.max(datanull)]


## Null SNPs plot . . . 

SNP_plot=df_filt[seq(from=1, to=nrow(df_filt), by = 800),]

plot(x=SNP_plot$row.names, y=SNP_plot$marker, col = colfunc(1500), pch =20, cex=5)




### Show null values acorss the genome
nulls<-nullsnps[snp,]
df_nulls<-df_filt[nulls,]


points(x=df_nulls$row.names, y=df_nulls$marker, col="red", pch=20, cex = 1)

# 
options(scipen=999)
message("Generating Sliding Windows")
slidingwindow<-windowseeker(snp, df_filt, row.names)
windiff<-df_filt[slidingwindow[,51],]$marker - df_filt[slidingwindow[,1],]$marker
if(snp > 1){
  datawg<-rowMeans(get_null(datawg, snp, slidingwindow))
  p2_sliding<-rowMeans(get_null(df_filt$p2, snp, slidingwindow))
  p3_sliding<-rowMeans(get_null(df_filt$p3, snp, slidingwindow))
  p4_sliding<-rowMeans(get_null(df_filt$p4, snp, slidingwindow))
  exp2_sliding<-rowMeans(get_null(df_filt$expected_p2, snp, slidingwindow))
  NGDHG_sliding<-rowMeans(get_null(df_filt$NGDHG, snp, slidingwindow))
  NGDNEO_sliding<-rowMeans(get_null(df_filt$NGDNEO, snp, slidingwindow))
  NGDMSL_sliding<-rowMeans(get_null(df_filt$NGDMSL, snp, slidingwindow))
  NGDMNEO_sliding<-rowMeans(get_null(df_filt$NGDMNEO, snp, slidingwindow))
  # rawF4_sliding<-rowMeans(get_null(df_filt$rawF4, snp, slidingwindow))
  HZdiff_sliding<-rowMeans(get_null(df_filt$HZdiff, snp, slidingwindow))
  slidingwindow_dist<-df_filt$marker[slidingwindow[,snp]]-df_filt$marker[slidingwindow[,1]]
  #datawg<-datawg/slidingwindow_dist
  
  df_filt$F2_MNEO_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_MNEOREF, snp, slidingwindow)), df_filt, 22)
  df_filt$F2_HG_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_HGREF, snp, slidingwindow)), df_filt, 22)
  df_filt$F2_NEO_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_NEOREF, snp, slidingwindow)), df_filt, 22)
  
  
  datawg2<-slidingwindow_NAccountant_test(datawg, df_filt, 22)
  HZdiff_sliding<-slidingwindow_NAccountant_test(df_filt$HZdiff, df_filt, 22)
  HZdiff2_sliding<-slidingwindow_NAccountant_test(df_filt$HZdiff2, df_filt, 22)
  #rawF4_sliding<-slidingwindow_NAccountant_test(rawF4_sliding, df_filt, 22)
  
  p2_sliding<-slidingwindow_NAccountant_test(p2_sliding, df_filt, 22)
  p3_sliding<-slidingwindow_NAccountant_test(p3_sliding, df_filt, 22)
  p4_sliding<-slidingwindow_NAccountant_test(p4_sliding, df_filt, 22)
  exp2_sliding<-slidingwindow_NAccountant_test(exp2_sliding, df_filt, 22)
  
  
  NGDHG_sliding<-slidingwindow_NAccountant_test(NGDHG_sliding, df_filt, 22)
  NGDNEO_sliding<-slidingwindow_NAccountant_test(NGDNEO_sliding, df_filt, 22)
  NGDMNEO_sliding<-slidingwindow_NAccountant_test(NGDMNEO_sliding, df_filt, 22)
  NGDMSL_sliding<-slidingwindow_NAccountant_test(NGDMSL_sliding, df_filt, 22) 
  
  
  df_filt$datawg2<-datawg2
  df_filt$sliding_p2<-p2_sliding
  df_filt$sliding_p3<-p3_sliding
  df_filt$sliding_p4<-p4_sliding
  df_filt$exp2_sliding<-exp2_sliding
  
  # df_filt$rawF4_sliding<-rawF4_sliding
  
  df_filt$NGDHG_sliding<-NGDHG_sliding
  df_filt$NGDNEO_sliding<-NGDNEO_sliding
  df_filt$NGDMNEO_sliding<-NGDMNEO_sliding
  df_filt$NGDMSL_sliding<-NGDMSL_sliding
  
  
  df_filt$HZdiff_sliding<-HZdiff_sliding
  df_filt$HZ2diff_sliding<-HZdiff2_sliding  
  
  df_filt_store<-df_filt
  df_filt<-df_filt[which(df_filt$datawg2!="NA"),]
  df_filt$datawg<-df_filt$datawg2
  df_filt$slidingwindowdist<-slidingwindow_dist
  datawg<-df_filt$datawg
  
}





WID=14

##### Pval derivation ------
G<-fitdist(datanull, "gamma", method = "mle", keepdata = T)
message("Fitting Genome-wide Null Distribution")
nullfitG<-fitdist(datanull, "gamma", keepdata=T, method ="mle")[1]



pval<-pgamma(datawg, shape=as.vector(nullfitG$estimate[1]), rate = as.vector(nullfitG$estimate[2]), lower.tail = F, log.p = T)
datapvalG<--pval/log(10)
max(datapvalG)


df_filt$datawg<-datawg
df_filt$pval<-datapvalG
df_filt$filtUID<-seq_along(df_filt$pval)
df_filt$row.names<-as.integer(df_filt$row.names)


df_filt2<-df_filt[,c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,58)]
x<-colnames(df_filt2)
x[3]<-"GT_REF"
x[4]<-"GT_MNEO"
x[5]<-"GT_HG"
x[6]<-"GT_NEO"
x[7]<-"p_REF"
x[8]<-"p_MNEO"
x[9]<-"p_HG"
x[10]<-"p_NEO"

colnames(df_filt2)<-x
#write.table(df_filt2, file="1240k_fadm.tsv",sep="\t", row.names=F)


max(datapvalG)

k1240<-df_filt[df_filt$row.names==6,]



# 

glist3<-glist2[-c(22:29),]
glist_arti<-glist3
glist_arti<-glist_arti[-c(15:21),]


glist_arti[1,]<-c(6,29691116,31478901, "Class I")
glist_arti[2,]<-c(6,32407618,32731330, "Class II")
glist_arti[3,]<-c(6,31539875,32003195, "Class III")
glist_arti<-glist_arti[-c(4:nrow(glist_arti)),]

glist_HLA<-glist3[-c(1:3),]
glist_HLA$uid<-seq_along(glist_HLA$V1)
glist_HLA$bool<-is.odd(glist_HLA$uid)#class II = 10120 -> 10081
#class I = 9989 -> 9916
#class III = 10057 -> 9998

k1240<-df_filt
k1240<-k1240[k1240$row.names==6,]

mappa<-fread("~/EvoFunc/CRGalign24mer_MHC.tsv")
colnames(mappa)<-c("M","Mstart","Mend","Map")
mappa$row.names<-6
df_mappa<-merge(df_filt, mappa, by = c("a","b"))



df_filt_markerprox<-df_filt
SPACER=100000
for(row.names in 1:21){
  df_current<-which(df_filt_markerprox$row.names == row.names)
  df_succ<-which(df_filt_markerprox$row.names == row.names+1)
  current_marker<-max(df_filt_markerprox[df_current,]$marker)
  print(current_marker)
  df_filt_markerprox[df_succ,]$marker <- df_filt_markerprox[df_succ,]$marker + SPACER + current_marker
  
}
df_filt$markerprox<-df_filt_markerprox$marker


plot_filter<-0
plot_max<-60
df<-df_filt
genomic_values<-df_filt$pval
bed=glist
outlier_threshold<--log10(5e-8)
row.names
chromosomes<-22
topoutliers<-return_outliers(df_filt, df_filt$pval, glist, outlier_threshold, plot_filter)
#topoutliers<-topoutliers[which(topoutliers$gene=="HLA-DQB1"),]
df_smol<-df[which(genomic_values > plot_filter),]
df_smol$uid<-seq_along(df_smol[,1])
GVsmol<-genomic_values[which(genomic_values > plot_filter)]
df_smol$smoluid<-seq_along(df_smol$uid)
df_smol$GVsmol<-GVsmol






mean(df$NGDMSL)
mean(df$NGDNEO)
mean(df$NGDHG)
mean(df$NGDMNEO)


coolsnp<-df_filt[which(df_filt$row.names==11 & df_filt$marker > 0 & df_filt$marker < 12337556),]
plot(coolsnp$datawg_nonsliding)


WID=14

tiff(filename="davyetal2022_Fig2B.tiff", units = "in", width=14, height=7, res=150, compression="lzw")


plot.new()
box()
plot.window(xlim = c(-5, max(df_smol$markerprox+10000)), ylim = c(plot_filter, plot_max))
#title(main="5M null distribution")
df_smol<-df_smol
abline(h= -log10(5e-8), col = "#e40959", lty = 2, lwd=2)
mbpeaks=NULL
for (x in 1:nrow(topoutliers)){
  print(nrow(topoutliers)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  row.namesstore<-topoutliers[x,]$row.names 
  markerstore<-topoutliers[x,]$marker
  mbpeaks<-df_smol[which(df_smol$row.names==row.namesstore & df_smol$GVsmol >= outlier_threshold & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  ypoints<-as.vector(unlist(mbpeaks$GVsmol))
  xpoints<-as.vector(unlist(mbpeaks$markerprox))
  #if(is.na(topoutliers[x,]$locus)){
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  points(x = xpoints, y = ypoints, pty=20, cex=0.7, pch = 16, col = "blue")
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}

#topoutlierG$datawg.sliding<-df_smol[topoutlierG$smolluid,]$datawg.sliding
#df$slidingorder = findInterval(df$slidingwindowdist, sort(df$slidingwindowdist))
#df$slidingwindowdistlog<-log10(df$slidingwindowdist)
#df$slidingorderlog = findInterval(df$slidingwindowdistlog, sort(df$slidingwindowdistlog))
#pal = colorRampPalette(c("blue", "red"))
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="grey",cex=as.numeric(is.odd(as.numeric(df_smol$row.names))), pch = 20)
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="goldenrod2",cex=as.numeric(is.even(as.numeric(df_smol$row.names))), pch = 20)

#points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col=df_filt$Col)
#legend("topright", col=pal(2), pch=19,legend=c(round(range(df$slidingorderlog10), 1)))


#### Add Axes ------



### Add Gene Names ----
if(max(genomic_values)>=plot_filter){
  coli<-0
  for (unigene in unique(topoutliers[,3]))
  {
    {
      toGuid<-which(topoutliers[,3]==unigene)
      tophits<-topoutliers[toGuid,]
      topgene<-unigene
      toppval<-tophits[which.max(tophits$pval),]$pval
      locus<-tophits[which.max(tophits$pval),]$markerprox
      coli<-coli+1
      if(coli%%15 == 0){
        coli<-coli+1
      }
      ### COL
      text(x=locus, y=(toppval), labels = as.character(topgene), cex=0.6, pos = 3)
    }
  }}

## Finish Manhattan -----
row.namesindex=NULL
for (x in 1:22){
  row.namesindex[x]<-tail(df_smol[(which(df_smol$row.names==x)),]$markerprox, 1)
}
row.namesindex2<-append(0, row.namesindex)
row.namesmidpoint<-row.namesindex2[-length(row.namesindex2)] + diff(row.namesindex2)/2

axis(1, at = row.namesmidpoint, labels=c(1:22), tick = FALSE, cex.axis=0.8)
axis(1, at = row.namesmidpoint, labels = FALSE, tick = TRUE)
axis(2, at=seq(from = 0, to = 60, by = 5))
title(xlab = "Chromosome", ylab = "P-value (-log10)", cex.main=0.8)
title(main="F_adm, V50, snp = 51")
#legend(x="topleft",legend=c(paste("0.05 Bonferroni Correction, n = ", nrow(df_filt)),"5e-8"), col = c("blue","red"), lty=2, cex =0.6)


dev.off()


#png(filename = paste0(filename, ".QQ.png"), width=700, height=700, res=300, pointsize = 4)
tiff(filename = paste0(filename, ".rebut.redomapp99.box.QQ.tiff"), width=5, height=5, units = "in", res=200, compression = "lzw")
shape = as.vector(head(nullfitG$estimate))[1]
rate = as.vector(head(nullfitG$estimate))[2]
data<-pgamma(datawg, shape = shape, rate=rate, lower.tail = F)
data<-data[data>0]
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- ceiling(max(logobspval))+3
expmax <- ceiling(max(logexppval))+3
plot(c(0,expmax), c(0,expmax), col="grey", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l")
box()
points(logexppval, logobspval, pch=20, col="goldenrod2", cex=1)

abline(a = 0, b = 1, col = "grey", lwd =4)
shape = as.vector(head(nullfitG$estimate))[1]
rate = as.vector(head(nullfitG$estimate))[2]
data<-pgamma(datanull, shape = shape, rate=rate, lower.tail=F)
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
points(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)
#title(main= "5M observed vs Expected pvals", xlab = "")
lines(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)
dev.off()



topoutliers=NULL

tiff(filename=paste0(filename,".w",WID, "redo.rebutmapp99.60.empty.hue.xqc.redo.berisa.selscan.tiff"), units = "in", width=14, height=7, res=150, compression="lzw")


plot.new()
box()
plot.window(xlim = c(-5, max(df_smol$markerprox+10000)), ylim = c(plot_filter, plot_max))
#title(main="5M null distribution")
df_smol<-df_smol
abline(h= -log10(5e-8), col = "#e40959", lty = 2, lwd=2)
mbpeaks=NULL
for (x in 1:nrow(topoutliers)){
  print(nrow(topoutliers)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  row.namesstore<-topoutliers[x,]$row.names 
  markerstore<-topoutliers[x,]$marker
  mbpeaks<-df_smol[which(df_smol$row.names==row.namesstore & df_smol$GVsmol >= outlier_threshold & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  ypoints<-as.vector(unlist(mbpeaks$GVsmol))
  xpoints<-as.vector(unlist(mbpeaks$markerprox))
  #if(is.na(topoutliers[x,]$locus)){
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  points(x = xpoints, y = ypoints, pty=20, cex=0.7, pch = 16, col = "blue")
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}

#topoutlierG$datawg.sliding<-df_smol[topoutlierG$smolluid,]$datawg.sliding
#df$slidingorder = findInterval(df$slidingwindowdist, sort(df$slidingwindowdist))
#df$slidingwindowdistlog<-log10(df$slidingwindowdist)
#df$slidingorderlog = findInterval(df$slidingwindowdistlog, sort(df$slidingwindowdistlog))
#pal = colorRampPalette(c("blue", "red"))
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="grey",cex=as.numeric(is.odd(as.numeric(df_smol$row.names))), pch = 20)
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="goldenrod2",cex=as.numeric(is.even(as.numeric(df_smol$row.names))), pch = 20)

#points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col=df_filt$Col)
#legend("topright", col=pal(2), pch=19,legend=c(round(range(df$slidingorderlog10), 1)))


#### Add Axes ------



### Add Gene Names ----
if(max(genomic_values)>=plot_filter){
  coli<-0
  for (unigene in unique(topoutliers[,3]))
  {
    {
      toGuid<-which(topoutliers[,3]==unigene)
      tophits<-topoutliers[toGuid,]
      topgene<-unigene
      toppval<-tophits[which.max(tophits$pval),]$pval
      locus<-tophits[which.max(tophits$pval),]$markerprox
      coli<-coli+1
      if(coli%%15 == 0){
        coli<-coli+1
      }
      ### COL
      text(x=locus, y=(toppval), labels = as.character(topgene), cex=0.6, pos = 3)
    }
  }}

## Finish Manhattan -----
row.namesindex=NULL
for (x in 1:22){
  row.namesindex[x]<-tail(df_smol[(which(df_smol$row.names==x)),]$markerprox, 1)
}
row.namesindex2<-append(0, row.namesindex)
row.namesmidpoint<-row.namesindex2[-length(row.namesindex2)] + diff(row.namesindex2)/2

axis(1, at = row.namesmidpoint, labels=c(1:22), tick = FALSE, cex.axis=0.8)
axis(1, at = row.namesmidpoint, labels = FALSE, tick = TRUE)
axis(2, at=seq(from = 0, to = 60, by = 5))
title(xlab = "Chromosome", ylab = "P-value (-log10)", cex.main=0.8)
title(main="F_adm, V50, snp = 51")
#legend(x="topleft",legend=c(paste("0.05 Bonferroni Correction, n = ", nrow(df_filt)),"5e-8"), col = c("blue","red"), lty=2, cex =0.6)


dev.off()




plot_filter<-0
plot_max<-60
df<-df_filt
genomic_values<-df_filt$pval
bed=glist
outlier_threshold<--log10(5e-8)
row.names
chromosomes<-22
topoutliers<-return_outliers(df_filt, df_filt$pval, glist, outlier_threshold, plot_filter)
topoutliers<-topoutliers[which(topoutliers$gene=="HLA-DQB1"),]
df_smol<-df[which(genomic_values > plot_filter),]
df_smol$uid<-seq_along(df_smol[,1])
GVsmol<-genomic_values[which(genomic_values > plot_filter)]
df_smol$smoluid<-seq_along(df_smol$uid)
df_smol$GVsmol<-GVsmol


tiff(filename=paste0(filename,".w",WID, "redo.rebutmapp99.60.fig1.hue.xqc.redo.berisa.selscan.tiff"), units = "in", width=14, height=7, res=150, compression="lzw")


plot.new()
box()
plot.window(xlim = c(-5, max(df_smol$markerprox+10000)), ylim = c(plot_filter, plot_max))
#title(main="5M null distribution")
df_smol<-df_smol
abline(h= -log10(5e-8), col = "#e40959", lty = 2, lwd=2)
mbpeaks=NULL
for (x in 1:nrow(topoutliers)){
  print(nrow(topoutliers)-x)
  #ypoints<-as.vector(unlist(mbpeaks$pval))
  #xpoints<-as.vector(unlist(mbpeaks$smoluid))
  row.namesstore<-topoutliers[x,]$row.names 
  markerstore<-topoutliers[x,]$marker
  mbpeaks<-df_smol[which(df_smol$row.names==row.namesstore & df_smol$GVsmol >= outlier_threshold & df_smol$marker>=(markerstore-500000) & df_smol$marker <= (markerstore+500000)),]
  
  ypoints<-as.vector(unlist(mbpeaks$GVsmol))
  xpoints<-as.vector(unlist(mbpeaks$markerprox))
  #if(is.na(topoutliers[x,]$locus)){
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  points(x = xpoints, y = ypoints, pty=20, cex=0.7, pch = 16, col = "blue")
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}

#topoutlierG$datawg.sliding<-df_smol[topoutlierG$smolluid,]$datawg.sliding
#df$slidingorder = findInterval(df$slidingwindowdist, sort(df$slidingwindowdist))
#df$slidingwindowdistlog<-log10(df$slidingwindowdist)
#df$slidingorderlog = findInterval(df$slidingwindowdistlog, sort(df$slidingwindowdistlog))
#pal = colorRampPalette(c("blue", "red"))
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="grey",cex=as.numeric(is.odd(as.numeric(df_smol$row.names))), pch = 20)
points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col="goldenrod2",cex=as.numeric(is.even(as.numeric(df_smol$row.names))), pch = 20)

#points(df_smol$markerprox,genomic_values[which(genomic_values > plot_filter)],col=df_filt$Col)
#legend("topright", col=pal(2), pch=19,legend=c(round(range(df$slidingorderlog10), 1)))


#### Add Axes ------



### Add Gene Names ----
if(max(genomic_values)>=plot_filter){
  coli<-0
  for (unigene in unique(topoutliers[,3]))
  {
    {
      toGuid<-which(topoutliers[,3]==unigene)
      tophits<-topoutliers[toGuid,]
      topgene<-unigene
      toppval<-tophits[which.max(tophits$pval),]$pval
      locus<-tophits[which.max(tophits$pval),]$markerprox
      coli<-coli+1
      if(coli%%15 == 0){
        coli<-coli+1
      }
      ### COL
      text(x=locus, y=(toppval), labels = as.character(topgene), cex=0.6, pos = 3)
    }
  }}

## Finish Manhattan -----
row.namesindex=NULL
for (x in 1:22){
  row.namesindex[x]<-tail(df_smol[(which(df_smol$row.names==x)),]$markerprox, 1)
}
row.namesindex2<-append(0, row.namesindex)
row.namesmidpoint<-row.namesindex2[-length(row.namesindex2)] + diff(row.namesindex2)/2

axis(1, at = row.namesmidpoint, labels=c(1:22), tick = FALSE, cex.axis=0.8)
axis(1, at = row.namesmidpoint, labels = FALSE, tick = TRUE)
axis(2, at=seq(from = 0, to = 60, by = 5))
title(xlab = "Chromosome", ylab = "P-value (-log10)", cex.main=0.8)
title(main="F_adm, V50, snp = 51")
#legend(x="topleft",legend=c(paste("0.05 Bonferroni Correction, n = ", nrow(df_filt)),"5e-8"), col = c("blue","red"), lty=2, cex =0.6)


dev.off()



