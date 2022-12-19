setwd("~/EvoFunc/ms_scripts")
rm(list=ls())
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
df<-as.data.frame(fread("/Users/davyt/slurmp/rebuttal/redo/agdp_fadm/v2/6.rebut_redo.fadm.HG.NEO.MNEO.MSL.agdp.nogenewiz.freqstats.minmaf005.notransversions.tsv.map99.bed"))

colnames(df)<-c("row.names","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
#df<-df[,-11]
#df<-df[-which(df$row.names==23),]
#df<-df[-which(df$row.names==24),]
#df$rawF4<-df$datawg
df$df.uid<-seq_along(df[,1])



df$expected_p2<-(df$p3*0.3)+(df$p4*0.7)
y<-df$expected_p2
df$datawg<-((df$p2-y)^2)/(2*(1-y^2-(1-y)^2))




df$datawg_nonsliding<-df$datawg
#df<-na.omit(df)
df<-df[-which(df$datawg == "Inf"),]
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0
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

df$n=df$GT1L+df$GT2L+df$GT3L

df$maf<-apply(df, 1, function(x){
  GTS=str_c(x[3], x[4], x[5])
  y=rawToChar(unique(charToRaw(GTS)))
  GT1=str_split(y, pattern="")[[1]][1]
  GT2=str_split(y, pattern="")[[1]][2]
  
  minor=min(c(str_count(GTS, GT1),str_count(GTS, GT2)))
  maf=minor/as.numeric(x[20])
  return(maf)
  
})



df$df.uid<-seq_along(df[,1])
#df<-df[-which(is.na(df$p2)),]
#df$expected_p2<-(df$p3*0.3)+(df$p4*0.7)

df$F2_NEOREF<-(df$p4-df$p1)^2
df$F2_MNEOREF<-(df$p2-df$p1)^2
df$F2_HGREF<-(df$p3-df$p1)^2


#df<-na.omit(df)
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

# plot(x=df_nulls$markerprox, y=datanull, axes=F, pch=20)
# 
# row.namesindex=NULL
# for (x in 1:22){
#   row.namesindex[x]<-tail(df_nulls[(which(df_nulls$row.names==x)),]$markerprox, 1)
# }
# row.namesindex2<-append(0, row.namesindex)
# row.namesmidpoint<-row.namesindex2[-length(row.namesindex2)] + diff(row.namesindex2)/2
# 
# axis(1, at = row.namesmidpoint, labels=c(1:22), tick = FALSE, cex.axis=0.8)
# axis(1, at = row.namesmidpoint, labels = FALSE, tick = TRUE)
# axis(2)
# #title(xlab = "Chromosome", ylab = "P-value (-log10)", cex.main=0.8)
# title(main="F_adm, V50, fixed NEO, no Greece_Peloponnesse_N, EHG")
# plot(density(datawg))
# plot(density(datanull))
# dev.off()
# 
# # Sliding Window -----
df_filt$row.names=1
options(scipen=999)
message("Generating Sliding Windows")
slidingwindow<-windowseeker(snp, df_filt, 1)

datawg<-rowMeans(get_null(datawg, snp, slidingwindow))
datawg<-c(rep("NA",25), datawg, c(rep("NA", 25)))
df_filt$datawg<-datawg
df_filt=df_filt[-which(datawg=="NA"),]
df_filt$datawg=as.numeric(df_filt$datawg)
AGDP_593<-df_filt
AGDP_593$row.names=6

glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist
glist2<-glist[grep(glist$V4, pattern="HLA"),]
glist<-glist2
glist3<-glist
glist15<-glist0[glist0$V1==15,]


AGDP_noMHC<-AGDP_593[which(AGDP_593$marker < 29500000 | AGDP_593$marker > 33000000),]
mean(AGDP_noMHC$datawg)

# gwasplot<-fread("/Users/davyt/slurmp/allgwas_standardized.6.tsv")
# colnames(gwasplot)<-c("chrom","pos","alt","ref","test","low_confidence","beta","pval","study","trait")


## LA
poschr<-fread("/Users/davyt/EvoFunc/local_ancestry/store/chrpos.allchr.txt")
setwd("/Users/davyt/slurmp/rebuttal/redo/hmm")

HG.ancestry<-fread("HG.means.txt")

NEO.ancestry<-fread("NEO.means.txt")

HET.ancestry<-fread("HET.means.txt")

LAdf<-cbind(poschr,HG.ancestry,HET.ancestry,NEO.ancestry)
colnames(LAdf)<-c("V1","V2","HG.ancestry","HET.ancestry","NEO.ancestry")

rm(HG,NEO,HET,HG.ancestry,NEO.ancestry,HET.ancestry)
LAdf$HGall<-LAdf$HG.ancestry+(LAdf$HET.ancestry/2)
HGsd=sd(LAdf$HGall)
HGnsd=sd(LAdf$HG.ancestry)
LAdf$HGallZ<-(LAdf$HGall-mean(LAdf$HGall))/HGsd
LAdf$HGZ<-(LAdf$HG.ancestry-mean(LAdf$HG.ancestry))/HGnsd

LAdf$HGZcomp<-(LAdf$HGallZ - LAdf$HGZ)


LAdf$row.names=LAdf$V1
LAdf$marker=LAdf$V2
LAdf<-as.data.frame(LAdf)
nullwind=5000000
nullsnps<-windowsampler(nullwind=5000000, LAdf, 22, 1)
null=LAdf[nullsnps,]
LAdf_iain_null<-null[,-c(10,11)]
colnames(LAdf_iain_null)[c(1,2)]<-c("row.names","marker")
#write.table(LAdf_iain_null, file="LAdf_iain_null_rebuttal.tsv", sep="\t", quote=F, col.names=T, row.names=F)
nullsd=sd(null$HGall)
null$HGallZ<-(null$HGall-mean(null$HGall))/nullsd
#plot(null$HGallZ)

LAdf$nullZ<-(LAdf$HGall - mean(null$HGall))/nullsd

LAdf_iain<-LAdf[,-c(10,11),]
#(LAdf_iain, file="LAdf_iain_rebuttal.tsv", quote=F, col.names=T, row.names=F, sep="\t")
colnames(LAdf_iain)[c(1,2)]<-c("row.names","marker")


LAdf6<-LAdf[LAdf$V1==6,]
LAdf15<-LAdf[LAdf$V1==15,]


delimit_MHC<-LAdf6
delimit_SLC<-LAdf15
#colnames(delimit_MHC)[c(1,2)]<-c("row.names","marker")

## 1240k 
# zoom-in of AGDP MW, 1240k Mw p-val & LA.
#rm(list=ls())

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


plotmax=30
LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
color.legend(
  xl=33950000,
  xr=34190000,
  yb = plotmax - 0.7*plotmax,
  yt = plotmax - 0.2*plotmax-1, gradient = "y" ,
  legend = seq(from= 0, to = 10, by=1)-5,
  rect.col = LAcol(11), cex = 0.6)
dev.off()


df<-as.data.frame(fread("/Users/davyt/slurmp/rebuttal/redo/REF_MNEO_HG_NEO_freqstats_mincount1_intersect.tsv"))
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

# # Sliding Window -----
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




datapvalG[which(datapvalG=="Inf")]<-20
df_filt$datawg<-datawg
df_filt$pval<-datapvalG
df_filt$filtUID<-seq_along(df_filt$pval)
df_filt$row.names<-as.integer(df_filt$row.names)

max(datapvalG)

k1240<-df_filt[df_filt$row.names==6,]



# 
glist<-read.table("~/EvoFunc/glist.txt")
glist0<-glist
glist2<-glist[grep(glist$V4, pattern="HLA"),]
glist<-glist2
glist3<-glist
glist15<-glist0[glist0$V1==15,]

glist3<-glist2[-c(22:29),]
glist_arti<-glist3
glist_arti<-glist_arti[-c(15:21),]


glist_arti[1,]<-c(6,29691116,31478901, "Class I")
glist_arti[2,]<-c(6,32407618,32731330, "Class II")
glist_arti[3,]<-c(6,31539875,32003195, "Class III")
glist_arti<-glist_arti[-c(4:nrow(glist_arti)),]

glist_HLA<-glist3[-c(1:3),]
glist_HLA$uid<-seq_along(glist_HLA$V1)
glist_HLA$bool<-is.odd(glist_HLA$uid)
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
plot_max<-55
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


glist3<-glist2[-c(22:29),]
glist_arti<-glist3
glist_arti<-glist_arti[-c(15:21),]


glist_arti[1,]<-c(6,29691116,31478901, "Class I")
glist_arti[2,]<-c(6,32407618,32731330, "Class II")
glist_arti[3,]<-c(6,31539875,32003195, "Class III")
glist_arti<-glist_arti[-c(4:nrow(glist_arti)),]

glist_HLA<-glist3[-c(1:3),]
glist_HLA$uid<-seq_along(glist_HLA$V1)
glist_HLA$bool<-is.odd(glist_HLA$uid)

k1240<-df_filt
k1240<-k1240[k1240$row.names==6,]

mappa<-fread("~/EvoFunc/CRGalign24mer_MHC.tsv")
colnames(mappa)<-c("M","Mstart","Mend","Map")
mappa$row.names<-6
df_mappa<-merge(df_filt, mappa, by = c("a","b"))




setwd("~/EvoFunc/ms_scripts")
tiff(filename="davyetal2022_fig2Cbeta.tiff", height= 7, units = "in", width = 14, res = 150, compression="lzw")
options(scipen=999)
plotlim<-c(29500000, 33000000)
par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
box()
plot.window(xlim = plotlim, ylim = c(-0.6,0.25))
grid(nx = NULL, ny=NULL)

lines(x=AGDP_593$marker, y=AGDP_593$datawg, col="darkblue", lwd=3,xlim = plotlim, pch=20)
lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)

apply(glist_arti, 1,function(x){
  rect(xleft=x[2], xright=x[3], ytop=-0.3, ybot=-0.40, col="grey85",border = NA)
})

apply(glist_HLA, 1, function(x){
  rect(xleft=as.numeric(x[2]), xright=as.numeric(x[3]), ybot = -0.40, ytop= -0.30, col="rosybrown3", border="rosybrown3")
  if(x[6]=="FALSE"){
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.285, labels = x[4], cex=0.75, adj=0)
  }else{
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.405,labels = x[4],cex=0.75, adj=1)
  }
})


pc<-0.0005*range(df_filt$marker)[2]
apply(glist_arti, 1,function(x){
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325, ybot=-0.375, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35, label=x[4], col="grey30")
})


legend(x="topleft", legend=c("Fadm (1240k)","Fadm (Shotgun)"), col=c("goldenrod2","darkblue"), lwd=2)

scale_limits<-c(-0.18,-0.045)
delimit_MHC$scale<-scales::rescale(delimit_MHC$HGall, to=c(scale_limits[1],scale_limits[2]))



points(y=delimit_MHC$scale, x=delimit_MHC$marker, type='l')

points(y=delimit_MHC$scale, x=delimit_MHC$marker)

points(y=delimit_MHC$scale, x=delimit_MHC$marker, col=(LAcol(9)[(trunc((delimit_MHC$nullZ)))+5]), type='p', lwd=1, pch=20)


meanline=scale_limits[1]-(((min(delimit_MHC$HGall))-mean(LAdf$HGall))/2.688357)



options(scipen=999)
axis(side=2, at =seq(from=0,to=0.5, by=0.05))
axis(1)

title(xlab="Marker (Chromosome 6)")
mtext(side = 2, text ='Fadm' , line=3, adj=0.85)

options(scipen=999)

max(delimit_MHC$HGall)-min(delimit_MHC$HGall) / scale_limits[2]- scale_limits[1]

start=scale_limits[1]-(((min(delimit_MHC$HGall))-0.15)/2.688357)
end=scale_limits[2]+((0.55-(max(delimit_MHC$HGall)))/2.688357)
axis(side=2, at=seq(from=start, to=end, by = (end-start)/10), labels=seq(from=0, to=1, by=0.1), cex.axis=0.65)
mtext(side=2, line=3, text="Mesolithic Ancestry (%)")



# box to sep Fadm/ LAD . .

LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
color.legend(
  xl=plotlim[1],
  xr=3e+07,
  yb = -0.065,
  yt = -0.025,
  gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

text(x=(plotlim[1]+3e+07)/2, y=(-0.045), labels = "Z-score", col="grey15")



dev.off()


##### REF PLOT ######## -----
## LA
poschr<-fread("/Users/davyt/EvoFunc/local_ancestry/store/chrpos.allchr.txt")
setwd("/Users/davyt/slurmp/rebuttal/redo/hmm_ref/")

HG.ancestry<-fread("localancestry.HG.allchr.txt")
HG.ancestry<-HG.ancestry[,3]

NEO.ancestry<-fread("localancestry.NEO.allchr.txt")
NEO.ancestry<-NEO.ancestry[,3]

HET.ancestry<-fread("localancestry.HET.allchr.txt")
HET.ancestry<-HET.ancestry[,3]

REFLAdf<-cbind(poschr,HG.ancestry,HET.ancestry,NEO.ancestry)
colnames(REFLAdf)<-c("V1","V2","HG.ancestry","HET.ancestry","NEO.ancestry")

rm(HG,NEO,HET,HG.ancestry,NEO.ancestry,HET.ancestry)
REFLAdf$HGall<-REFLAdf$HG.ancestry+(REFLAdf$HET.ancestry/2)
HGsd=sd(REFLAdf$HGall)
HGnsd=sd(REFLAdf$HG.ancestry)
REFLAdf$HGallZ<-(REFLAdf$HGall-mean(REFLAdf$HGall))/HGsd
REFLAdf$HGZ<-(REFLAdf$HG.ancestry-mean(REFLAdf$HG.ancestry))/HGnsd

REFLAdf$HGZcomp<-(REFLAdf$HGallZ - REFLAdf$HGZ)
colnames(REFLAdf)[c(1,2)]<-c("row.names","marker")
REFLAdf<-REFLAdf[REFLAdf$row.names==6,]

setwd("~/EvoFunc/ms_scripts")


tiff(filename="davyetal2022_figSF8.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
options(scipen=999)
plotlim<-c(29500000, 33000000)
par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
box()
plot.window(xlim = plotlim, ylim = c(-1.2,1.25))
grid(nx = NULL, ny=NULL)

lines(x=AGDP_593$marker, y=AGDP_593$datawg, col="darkblue", lwd=3)
lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)

apply(glist_arti, 1,function(x){
  rect(xleft=x[2], xright=x[3], ytop=-0.3-0.5, ybot=-0.40-0.5, col="grey85",border = NA)
})

apply(glist_HLA, 1, function(x){
  rect(xleft=as.numeric(x[2]), xright=as.numeric(x[3]), ybot = -0.40-0.5, ytop= -0.30-0.5, col="rosybrown3", border="rosybrown3")
  if(x[6]=="FALSE"){
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.285-0.5, labels = x[4], cex=0.75, adj=0)
  }else{
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.405-0.5,labels = x[4],cex=0.75, adj=1)
  }
})


pc<-0.0005*range(df_filt$marker)[2]
apply(glist_arti, 1,function(x){
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325-0.5, ybot=-0.375-0.5, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35-0.5, label=x[4], col="grey30")
})


legend(x=29350000,y=0.4, legend=c("Fadm (1240k)","Fadm (Shotgun)"), col=c("goldenrod2","darkblue"), lwd=2)

scale_limits<-c(-0.18,-0.045)
delimit_MHC$scale<-scales::rescale(delimit_MHC$HGall, to=c(scale_limits[1],scale_limits[2]))


options(scipen=999)
axis(side=2, at =seq(from=0,to=0.3, by=0.1))
axis(1)

title(xlab="Marker (Chromosome 6)")
mtext(side = 2, text =expression("F"[adm]) , line=3, adj=0.6)

options(scipen=999)

max(delimit_MHC$HGall)-min(delimit_MHC$HGall) / scale_limits[2]- scale_limits[1]

lines(k1240$F2_MNEO_sliding+0.5, x=k1240$marker, lty=1, lwd=2, pch=20)

lines(k1240$F2_HG_sliding+0.5, x=k1240$marker, col="blue")

lines(k1240$F2_NEO_sliding+0.5,x=k1240$marker, col="red")

axis(side=2, at =seq(from=0+0.5,to=1+0.5, by=0.2), labels = seq(from=0,to=1,by=0.2))
mtext(side = 2, text =expression("F"[2]) , line=3, adj=0.8)



legend(x="topleft", legend=c(expression('F'[2]*' (Neolithic:Reference)'),expression("F"[2]* ' (Mesolithic:Reference)'),expression("F"[2]*' (Admixed Neolithic:Reference)')), col=c("red","blue","black"), lty=1, lwd=2)



# REF LAD

REFLAdf$rescaled<-scales::rescale(x = REFLAdf$HGall, to=c(-0.45, -0.15))

points(x=REFLAdf$marker, y=REFLAdf$rescaled, type='l')

points(y=REFLAdf$rescaled, x=REFLAdf$marker)

points(y=REFLAdf$rescaled, x=REFLAdf$marker, col=(LAcol(11)[(round(REFLAdf$HGall*10, digits = 0)+1)]), type='p', lwd=1, pch=20)


axis(side=2, at=seq(from=-0.45, to=-0.15, by=0.3/4), labels=seq(from=0, to = 100 , by =100/4), cex.lab=0.8)

mtext(side = 2, text = "Mesolithic Ancestry %" , line=3, adj=0.32)

dev.off()

k1240<-df_filt[df_filt$row.names==15,]


AGDP_593_SLC<-fread("/Users/davyt/slurmp/rebuttal/redo/agdp_refilt/15.fadm.HG.NEO.MNEO.MSL.agdp.nogenewiz.freqstats.minmaf005.notransversions.tsv.map99.bed")

colnames(AGDP_593_SLC)<-c("row.names","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
AGDP_593_SLC$df.uid<-seq_along(AGDP_593_SLC[,1])



AGDP_593_SLC$expected_p2<-(df$p3*0.3)+(df$p4*0.7)
y<-AGDP_593_SLC$expected_p2
AGDP_593_SLC$datawg<-((df$p2-y)^2)/(2*(1-y^2-(1-y)^2))




AGDP_593_SLC$datawg_nonsliding<-df$datawg
AGDP_593_SLC<-na.omit(df)
AGDP_593_SLC<-AGDP_593_SLC[-which(AGDP_593_SLC),]
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0





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

df$n=df$GT1L+df$GT2L+df$GT3L

df$maf<-apply(df, 1, function(x){
  GTS=str_c(x[3], x[4], x[5])
  y=rawToChar(unique(charToRaw(GTS)))
  GT1=str_split(y, pattern="")[[1]][1]
  GT2=str_split(y, pattern="")[[1]][2]
  
  minor=min(c(str_count(GTS, GT1),str_count(GTS, GT2)))
  maf=minor/as.numeric(x[20])
  return(maf)
  
})


AGDP_593_SLC<-na.omit(AGDP_593_SLC)
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0
glist<-read.table("/Users/davyt/EvoFunc/glist.txt")

# # Sliding Window -----
df$row.names==1
df_filt$row.names=1
options(scipen=999)
message("Generating Sliding Windows")
slidingwindow<-windowseeker(snp, df_filt, 1)

datawg<-rowMeans(get_null(datawg, snp, slidingwindow))
datawg<-c(rep("NA",25), datawg, c(rep("NA", 25)))
df_filt$datawg<-datawg



mean(df$NGDMSL)
mean(df$NGDNEO)
mean(df$NGDHG)
mean(df$NGDMNEO)





glist0<-glist
glist<-glist0


glist_SLC<-glist0[which(glist0$V4=="SLC24A5"),]
glist00<-glist0[sample(size=10000, x=1:nrow(glist0), replace = F),]

glist000<-rbind(glist00,glist_SLC)

glist00<-glist00[order(glist00[,1], glist00[,2] ),]
glist000$uid<-seq_along(glist000$V1)
glist000$bool<-is.odd(glist000$uid)


glist0$uid<-seq_along(glist0$V1)
glist0$bool<-is.odd(glist0$uid)
glist0<-glist0[-which(glist0$V4=="NDUFAF4P1"),]
setwd("~/EvoFunc/ms_scripts")
tiff(filename="davyetal2022_figS5.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")

plotlim<-c(48413168-1000000, 48513168+1000000)
par(mar=c(6.1,5.1,4.1,5.1))
plot.new()
plot.window(xlim = plotlim, ylim = c(-0.5,0.3))
grid(nx = NULL, ny=NULL)

lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)



pc<-0.0005*range(df_filt$marker)[2]
apply(glist_arti, 1,function(x){
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325, ybot=-0.375, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35, label=x[4], col="grey30")
})


legend(x="topleft", legend=c("Fadm (1240k)"), col=c("goldenrod2"), lwd=2)

scale_limits<-c(-0.18,-0.045)
delimit_SLC$scale<-scales::rescale(delimit_SLC$HGall, to=c(scale_limits[1],scale_limits[2]))


points(y=delimit_SLC$scale, x=delimit_SLC$marker, type='l')

points(y=delimit_SLC$scale, x=delimit_SLC$marker)
points(y=delimit_SLC$scale, x=delimit_SLC$marker, col=(LAcol(9)[(trunc((delimit_SLC$nullZ)))+5]), type='p', lwd=1, pch=20)

meanline=scale_limits[1]-(((min(LAdf$HGall))-mean(LAdf$HGall))/2.688357)

abline(h=meanline, lty=2, col="Indianred4", lwd=2)

options(scipen=999)
axis(side=2, at =seq(from=0,to=0.3, by=0.05))
axis(1, at=seq(from=47500000, to=47500000+(5*500000), by=500000))

title(xlab="Marker (Chromosome 15)")
mtext(side = 2, text ='Fadm' , line=3, adj=0.85)

options(scipen=999)

max(delimit_SLC$HGall)-min(delimit_SLC$HGall) / scale_limits[2]- scale_limits[1]


start=scale_limits[1]-(((min(delimit_SLC$HGall))-0.15)/2.688357)
end=scale_limits[2]+((0.55-(max(delimit_SLC$HGall)))/2.688357)
axis(side=2, at=seq(from=start, to=end, by = (end-start)/8), labels=seq(from=0.15, to=0.55, by=0.05), cex.axis=0.65)
mtext(side=2, line=3, text="Mesolithic Ancestry (%)")


LAcol<-colorRampPalette(colors = c('#e40959','grey','#412ca8'))
color.legend(
  xl=plotlim[1],
  xr=plotlim[1]+((plotlim[2]-plotlim[1])/8),
  yb = -0.065,
  yt = -0.025,
  gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

x.text=(plotlim[1] + ((plotlim[1]+((plotlim[2]-plotlim[1])/8))))/2




text(x=x.text, y=(-0.045), labels = "LAD Z-score", col="grey15")

apply(glist0, 1, function(x){
  if(x[1]==15){
    rect(xleft=as.numeric(x[2]), xright=as.numeric(x[3]), ybot = -0.40, ytop= -0.30, col="grey90", border="black")
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.285, labels = x[4], cex=0.75, adj=0)

  }
})



dev.off()

