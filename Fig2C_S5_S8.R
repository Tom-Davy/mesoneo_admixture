# zoom-in of AGDP MW, 1240k Mw p-val & LA.
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


a=0.822222222
b=0.929824561
c=0.75

y<-(a*0.3)+(b*0.7)
x<-((c-y)^2)/(2*(1-y^2-(1-y)^2))



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
df_filt$row.names==1
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
x<-df$expected_p2
df$datawg<-((df$p2-x)^2)/(2*x*(1-x))
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
  
  # df_filt$rawF4_sliding<-rawF4_sliding
  
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
glist_HLA$bool<-is.odd(glist_HLA$uid)#class II = 10120 -> 10081
#class I = 9989 -> 9916
#class III = 10057 -> 9998

k1240<-df_filt
k1240<-k1240[k1240$row.names==6,]

mappa<-fread("~/EvoFunc/CRGalign24mer_MHC.tsv")
colnames(mappa)<-c("M","Mstart","Mend","Map")
mappa$row.names<-6
df_mappa<-merge(df_filt, mappa, by = c("a","b"))





## AGDP/1240k x-ver





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





#fix glist; 
#fix delimit_MHC
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





tiff(filename="v2AGDPFIXrebut_map99_agdp51_sortfixLAREF_253_ref_pres4_p4_golden_hues_1240kvsAGDPsnpwindows_pruned51_HLA_box_v2.tiff", height= 7, units = "in", width = 14, res = 150, compression="lzw")
options(scipen=999)
plotlim<-c(29500000, 33000000)
par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
box()
plot.window(xlim = plotlim, ylim = c(-0.6,0.6))
grid(nx = NULL, ny=NULL)
#rect(xleft = 0, xright=30300000000, ytop=end+0.035, ybot=start, col=add.alpha("grey20", 0.5))
#abline(h=c(end+0.035), lty=1)

lines(x=AGDP_593$marker, y=AGDP_593$datawg, col="darkblue", lwd=3,xlim = plotlim, pch=20)
lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)
#abline(h=mean(AGDP_noMHC$datawg), col="red", lty=2, lwd=2)
#abline(h=mean(AGDP_MHC$datawg), col="black")

apply(glist_arti, 1,function(x){
  rect(xleft=x[2], xright=x[3], ytop=-0.3, ybot=-0.40, col="grey85",border = NA)
  #text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.3, label=x[4], col="grey85")
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
  #rect(xleft=(as.numeric(x[2])+as.numeric(x[3]))/3, xright=2*((as.numeric(x[2])+as.numeric(x[3]))/3), ytop=-0.3, ybot=-0.4, col="grey50",border = NA)
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325, ybot=-0.375, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35, label=x[4], col="grey30")
})


legend(x="topleft", legend=c("Fadm (1240k)","Fadm (Shotgun)"), col=c("goldenrod2","darkblue"), lwd=2)

# LA but rescaled - scale_limits[1] to scale_limits[2]?
scale_limits<-c(-0.18,-0.045)
delimit_MHC$scale<-scales::rescale(delimit_MHC$HGall, to=c(scale_limits[1],scale_limits[2]))



points(y=delimit_MHC$scale, x=delimit_MHC$marker, type='l')

points(y=delimit_MHC$scale, x=delimit_MHC$marker)

points(y=delimit_MHC$scale, x=delimit_MHC$marker, col=(LAcol(9)[(trunc((delimit_MHC$nullZ)))+5]), type='p', lwd=1, pch=20)




#abline(h=(mean(LAdf$HGall)/(mean(LAdf$HGall)+min(delimit_MHC$HGall)))*(scale_limits[1]-scale_limits[2]), lty=2, col="indianred4", lwd=2)

meanline=scale_limits[1]-(((min(delimit_MHC$HGall))-mean(LAdf$HGall))/2.688357)

#abline(h=meanline, lty=2, col="Indianred4", lwd=2)

options(scipen=999)
axis(side=2, at =seq(from=0,to=0.5, by=0.05))
axis(1)

title(xlab="Marker (Chromosome 6)")
mtext(side = 2, text ='Fadm' , line=3, adj=0.85)

options(scipen=999)
#axis(side=4, at=(seq(from=scale_limits[1],to=scale_limits[2], by=-(scale_limits[1]-scale_limits[2])/5)), labels = seq(from=min(delimit_MHC$HGall), to=max(delimit_MHC$HGall), by=(max(delimit_MHC$HGall)-min(delimit_MHC$HGall))/5))

max(delimit_MHC$HGall)-min(delimit_MHC$HGall) / scale_limits[2]- scale_limits[1]
#so we change by 0.3629282 LAD in the space of 0.135 axis units.
# #so one point of change in axis is 2.688357 in LAD
# abline(h=scale_limits[1]-(0.03/2.688357), col="green")
# abline(h=scale_limits[2]+(0.01/2.688357), col="green")

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




# df_mappa$scale<-rescale(df_mappa$Map, to=c(0.3,0.4))
# points(x=df_mappa$marker, y=df_mappa$scale, type='l')
# 
# 
# axis(2, at=c(0.3,0.35,0.4), labels = c(0,0.5,1))
# mtext(side=2, line=3, adj=0.965, text = "Mappability")

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



tiff(filename="agdp51_redo2_REBUTTAL_REF_res4_p4_golden_hues_1240kvsAGDPsnpwindows_pruned5_HLA_box_v2.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
options(scipen=999)
plotlim<-c(29500000, 33000000)
par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
box()
plot.window(xlim = plotlim, ylim = c(-1.2,1.25))
grid(nx = NULL, ny=NULL)
#rect(xleft = 0, xright=30300000000, ytop=end+0.035, ybot=start, col=add.alpha("grey20", 0.5))
#abline(h=c(end+0.035), lty=1)

lines(x=AGDP_593$marker, y=AGDP_593$datawg, col="darkblue", lwd=3)
lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)

apply(glist_arti, 1,function(x){
  rect(xleft=x[2], xright=x[3], ytop=-0.3-0.5, ybot=-0.40-0.5, col="grey85",border = NA)
  #text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.3, label=x[4], col="grey85")
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
  #rect(xleft=(as.numeric(x[2])+as.numeric(x[3]))/3, xright=2*((as.numeric(x[2])+as.numeric(x[3]))/3), ytop=-0.3, ybot=-0.4, col="grey50",border = NA)
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325-0.5, ybot=-0.375-0.5, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35-0.5, label=x[4], col="grey30")
})


legend(x=29350000,y=0.4, legend=c("Fadm (1240k)","Fadm (Shotgun)"), col=c("goldenrod2","darkblue"), lwd=2)

# LA but rescaled - scale_limits[1] to scale_limits[2]?
scale_limits<-c(-0.18,-0.045)
delimit_MHC$scale<-scales::rescale(delimit_MHC$HGall, to=c(scale_limits[1],scale_limits[2]))

#points(y=delimit_MHC$scale, x=delimit_MHC$marker, col=(LAcol(9)[(trunc((delimit_MHC$nullZ)))+5]), type='p', lwd=1, pch=20)

#abline(h=(mean(LAdf$HGall)/(mean(LAdf$HGall)+min(delimit_MHC$HGall)))*(scale_limits[1]-scale_limits[2]), lty=2, col="indianred4", lwd=2)

#meanline=scale_limits[1]-(((min(delimit_MHC$HGall))-mean(LAdf$HGall))/2.688357)

#abline(h=meanline, lty=2, col="Indianred4", lwd=2)

options(scipen=999)
axis(side=2, at =seq(from=0,to=0.3, by=0.1))
axis(1)

title(xlab="Marker (Chromosome 6)")
mtext(side = 2, text =expression("F"[adm]) , line=3, adj=0.6)

options(scipen=999)
#axis(side=4, at=(seq(from=scale_limits[1],to=scale_limits[2], by=-(scale_limits[1]-scale_limits[2])/5)), labels = seq(from=min(delimit_MHC$HGall), to=max(delimit_MHC$HGall), by=(max(delimit_MHC$HGall)-min(delimit_MHC$HGall))/5))

max(delimit_MHC$HGall)-min(delimit_MHC$HGall) / scale_limits[2]- scale_limits[1]

lines(k1240$F2_MNEO_sliding+0.5, x=k1240$marker, lty=1, lwd=2, pch=20)

lines(k1240$F2_HG_sliding+0.5, x=k1240$marker, col="blue")

lines(k1240$F2_NEO_sliding+0.5,x=k1240$marker, col="red")
#rect(xleft = 0, xright=30300000000, ytop=end+0.035, ybot=start, col=add.alpha("grey20", 0.5))
#abline(h=c(end+0.035), lty=1)


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

#
#
## SLC tiff
k1240<-df_filt[df_filt$row.names==15,]


AGDP_593_SLC<-fread("/Users/davyt/slurmp/rebuttal/redo/agdp_refilt/15.fadm.HG.NEO.MNEO.MSL.agdp.nogenewiz.freqstats.minmaf005.notransversions.tsv.map99.bed")

colnames(AGDP_593_SLC)<-c("row.names","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
#df<-df[,-11]
#df<-df[-which(df$row.names==23),]
#df<-df[-which(df$row.names==24),]
#df$rawF4<-df$datawg
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
#df<-df[-which(df$datawg == "Inf"),]
snp<-51
NEOFILT<-0
row.names<-22
NOTNEOFILT<-0

#### GT lengths -----
glist<-read.table("/Users/davyt/EvoFunc/glist.txt")

# Nei's Gene Diversity Test


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
#glist0=glist0[-which(glist0$V4=="SLC24A5"),]
glist00<-glist0[sample(size=10000, x=1:nrow(glist0), replace = F),]

glist000<-rbind(glist00,glist_SLC)

glist00<-glist00[order(glist00[,1], glist00[,2] ),]
glist000$uid<-seq_along(glist000$V1)
glist000$bool<-is.odd(glist000$uid)


glist0$uid<-seq_along(glist0$V1)
glist0$bool<-is.odd(glist0$uid)
glist0<-glist0[-which(glist0$V4=="NDUFAF4P1"),]

tiff(filename="rebutt_map_minmaf_notvsion_p4_mapp_hues_1240kvsAGDPsnpwindows_pruned5_SLC.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")

plotlim<-c(48413168-1000000, 48513168+1000000)
par(mar=c(6.1,5.1,4.1,5.1))
plot.new()
plot.window(xlim = plotlim, ylim = c(-0.5,0.3))
grid(nx = NULL, ny=NULL)

#lines(x=AGDP_593_SLC$marker, y=AGDP_593_SLC$datawg, col="darkblue", lwd=3)
lines(x=k1240$marker, y=k1240$datawg, col="goldenrod2",lwd=3)

#apply(glist_arti, 1,function(x){
#  rect(xleft=x[2], xright=x[3], ytop=-0.3, ybot=-0.40, col="grey85",border = NA)
#text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.3, label=x[4], col="grey85")
#})



pc<-0.0005*range(df_filt$marker)[2]
apply(glist_arti, 1,function(x){
  #rect(xleft=(as.numeric(x[2])+as.numeric(x[3]))/3, xright=2*((as.numeric(x[2])+as.numeric(x[3]))/3), ytop=-0.3, ybot=-0.4, col="grey50",border = NA)
  rect(xleft=((as.numeric(x[2])+as.numeric(x[3]))/2)-pc, xright=((as.numeric(x[2])+as.numeric(x[3]))/2)+pc, ytop=-0.325, ybot=-0.375, col="grey87", border = "grey55")
  text(x=(as.numeric(x[2])+as.numeric(x[3]))/2, y=-0.35, label=x[4], col="grey30")
})


legend(x="topleft", legend=c("Fadm (1240k)"), col=c("goldenrod2"), lwd=2)

# LA but rescaled - scale_limits[1] to scale_limits[2]?
scale_limits<-c(-0.18,-0.045)
delimit_SLC$scale<-scales::rescale(delimit_SLC$HGall, to=c(scale_limits[1],scale_limits[2]))


points(y=delimit_SLC$scale, x=delimit_SLC$marker, type='l')

points(y=delimit_SLC$scale, x=delimit_SLC$marker)
points(y=delimit_SLC$scale, x=delimit_SLC$marker, col=(LAcol(9)[(trunc((delimit_SLC$nullZ)))+5]), type='p', lwd=1, pch=20)

#points(y=REFLAdf$rescaled, x=REFLAdf$marker, col=(LAcol(11)[(round(REFLAdf$HGall*10, digits = 0)+1)]), type='p', lwd=1, pch=20)

#abline(h=(mean(LAdf$HGall)/(mean(LAdf$HGall)+min(delimit_MHC$HGall)))*(scale_limits[1]-scale_limits[2]), lty=2, col="indianred4", lwd=2)

meanline=scale_limits[1]-(((min(LAdf$HGall))-mean(LAdf$HGall))/2.688357)

abline(h=meanline, lty=2, col="Indianred4", lwd=2)

options(scipen=999)
axis(side=2, at =seq(from=0,to=0.3, by=0.05))
axis(1, at=seq(from=47500000, to=47500000+(5*500000), by=500000))

title(xlab="Marker (Chromosome 15)")
mtext(side = 2, text ='Fadm' , line=3, adj=0.85)

options(scipen=999)
#axis(side=4, at=(seq(from=scale_limits[1],to=scale_limits[2], by=-(scale_limits[1]-scale_limits[2])/5)), labels = seq(from=min(delimit_MHC$HGall), to=max(delimit_MHC$HGall), by=(max(delimit_MHC$HGall)-min(delimit_MHC$HGall))/5))



max(delimit_SLC$HGall)-min(delimit_SLC$HGall) / scale_limits[2]- scale_limits[1]
#so we change by 0.3629282 LAD in the space of 0.135 axis units.
# #so one point of change in axis is 2.688357 in LAD
# abline(h=scale_limits[1]-(0.03/2.688357), col="green")
# abline(h=scale_limits[2]+(0.01/2.688357), col="green")

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
    #if(x[6]=="FALSE"){
    text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.285, labels = x[4], cex=0.75, adj=0)
    # }else{
    #  text(x = (as.numeric(x[2]) +as.numeric(x[3]) )/2, srt=90, y=-0.405,labels = x[4],cex=0.75, adj=1)
  }
})


# 
# df_mappa$scale<-rescale(df_mappa$Map, to=c(0.3,0.4))
# points(x=df_mappa$marker, y=df_mappa$scale, type='l')
# 
# 
# axis(2, at=c(0.3,0.35,0.4), labels = c(0,0.5,1))
# mtext(side=2, line=3, adj=0.965, text = "Mappability")

dev.off()
## SLC tiff


k1240<-df_filt[df_filt$row.names==15,]
glist<-glist0

###############################
#pdf("1240kvsAGDPsnpwindows_pruned2.pdf",width=12, height = 8)

#pdf("1240kvsAGDPsnpwindows_pruned2.pdf",width=12, height = 8)
AGDP_593_15<-fread("/Users/davyt/slurmp/autosomes.MW.HG.NEO.MNEO.MSL.agdp.tsv.minMAC3.minpopNEO30.nogenewiz551_SLC24A5zoomin.tsv")
AGDP_593<-AGDP_593_15
delimit_MHC<-LAdf15


tiff(filename="1240kvsAGDPsnpwindows_pruned2_SLC.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
par(mar=c(5, 4, 4, 4) + 0.1)
m=48434869
plot(AGDP_593$datawg, x =AGDP_593$marker, ylim =c(-0.15,0.3), xlim=c(m-3e+06,m+3e+06), xlab = "Marker (Chr 15)", ylab="Fadm", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp SLC24A5 signal", type='p', lwd=5, col="#d29618")


#plot(k1240$datawg, x =k1240$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight raw value", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")

points(k1240$datawg, x = k1240$marker, col = "goldenrod2", type='l', lwd=3, cex=0.5)
points(k1240$datawg, x = k1240$marker, col = "black", type='p', lwd=3, cex=0.5)

legend(x="topleft", legend=c("1240k Array - 21 snp window","Whole-genome shotgun - 551 snp window"), col = c("black","grey"), lty=c(1), lwd=3, cex=0.8)


axis(1)
axis(2, at=c(0,0.05,0.1,0.15,0.2,0.25,0.3), cex.axis=0.7)

#dev.off()

axis(2, at=c(-0.06125, -0.1125,-(0.15+0.125)/2), labels = c("GENES","1240k SNPs","Shotgun SNPs"), cex.axis=0.5, las=1)

abline(h=scale_limits[2], lty=1) #genes upper
abline(h=-0.1, lty=1) # genes lower
abline(h=-0.125) # 1240k
abline(h=-0.15) #

TRACKSPACERS=c(scale_limits[2],-0.1,-0.125,-0.15)
#abline(h=TRACKSPACERS[1], col="orange")
#abline(h=TRACKSPACERS[2], col="green")
#abline(h=TRACKSPACERS[4], col="red")
#abline(h=TRACKSPACERS[3], col="blue")
## SNPS: 1-2
rect(xleft=AGDP_593$marker,ybottom=TRACKSPACERS[4]+0.0015, ytop=TRACKSPACERS[3]-0.0015, xright=AGDP501$marker+100)
abline(h=-1, lty=2)



rect(xleft=k1240$marker,ybottom=TRACKSPACERS[3]+0.0015, ytop=TRACKSPACERS[2]-0.0015, xright=k1240$marker+100)
abline(h=-1, lty=2)
## GWAS: 2-3
## GWAS: 2-3
gwascol<-colorRampPalette(colors = c('red','blue'))

#genes
glistdiff<-TRACKSPACERS[1]-TRACKSPACERS[2]
glistjitter<-6
glistbool.tmp<-rep_len(TRACKSPACERS[1], length.out  = glistjitter-2)
tmp2<-c(1:(glistjitter-2))
glistbool.tmp2<-glistbool.tmp-(tmp2*(glistdiff/glistjitter))


glist15$BOOL<-rep_len(x = glistbool.tmp2,length.out = nrow(glist15))
#abline(h=unique(glist$BOOL), lty=2, col="red")



apply(glist15, MARGIN = 1, function(x){
  if(x[1]==15){
    rect(xleft=x[2], xright=x[3], ybottom=(as.numeric(x[5])-(glistdiff/glistjitter*2)), ytop=(as.numeric(x[5])-(2*(glistdiff/(glistjitter*2)))), col = "gray50", border = "gray37")
    text(x = ((as.numeric(x[2])+as.numeric(x[3]))/2), y=as.numeric(x[5]), labels = as.character(x[4]), cex=0.4, col = "grey5")
  }
})

par(new=TRUE)



delimit_marker<-delimit_MHC[,2]
delimit_MHC$marker<-delimit_MHC$V2


delimit_MHC_pos<-delimit_MHC[which(delimit_MHC$HGall >= 0),]
plot(x=delimit_MHC_pos$marker, y=delimit_MHC_pos$HGall, ylim=c(-1.6,0.9), xlim=c(m-3e+06,m+3e+06), type='b', axes=F, ylab=NA, xlab=NA, col=(LAcol(10)[(floor((delimit_MHC_pos$HGallZ))+6)]), pch=20) 

delimit_MHC_neg<-delimit_MHC[which(delimit_MHC$HGall < 0),]
points(x=delimit_MHC_neg$marker, y=delimit_MHC_neg$HGall, pch=20, col=(LAcol(10)[(ceiling((delimit_MHC_neg$HGallZ))+6)]))


plotmax=0.9
LAcol<-colorRampPalette(colors = c('red','grey','royalblue4'))
color.legend(
  xl=5.1e+07,
  xr=m+3e+06,
  yb = 0.2,
  yt = 0.85,
  gradient = "y" ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

text(x=(33950000+34190000)/2, y=(plotmax), labels = "Z-score")



axis(4, at = c(0,0.2,0.4,0.6,0.8,1))
mtext(side = 4, text = "% Neolithic Ancestry", line=2, adj=0.875)


axis(4, at=c(0.3,0.35,0.4), labels = c(0,0.5,1))
mtext(side=4, line=3, adj=0.965, text = "Mappability")

dev.off()

## round negative and positive sepe
tiff("1240kvsAGDPsnpwindows_pruned2.pdf",width=12, height = 8)
tiff(filename="1240kvsAGDPsnpwindows_pruned2.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
par(mar=c(5, 4, 4, 4) + 0.1)
plot(AGDP_593$datawg, x =AGDP_593$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight statistic", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")


#plot(k1240$datawg, x =k1240$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight raw value", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")

points(k1240$datawg, x = k1240$marker, col = "black", type='l', lwd=3, cex=0.5)
points(k1240$datawg, x = k1240$marker, col = "black", type='p', lwd=3, cex=0.5)

legend(x="topleft", legend=c("1240k Array - 21 snp window","Whole-genome shotgun - 551 snp window"), col = c("black","grey"), lty=c(1), lwd=3, cex=0.8)


axis(1)
axis(2, at=c(0,0.05,0.1,0.15,0.2,0.25,0.3), cex.axis=0.7)

#dev.off()

axis(2, at=c(-0.06125, -0.1125,-(0.15+0.125)/2), labels = c("GENES","1240k SNPs","Shotgun SNPs"), cex.axis=0.5, las=1)

abline(h=scale_limits[2], lty=1) #genes upper
abline(h=-0.1, lty=1) # genes lower
abline(h=-0.125) # 1240k
abline(h=-0.15) #

TRACKSPACERS=c(scale_limits[2],-0.1,-0.125,-0.15)
#abline(h=TRACKSPACERS[1], col="orange")
#abline(h=TRACKSPACERS[2], col="green")
#abline(h=TRACKSPACERS[4], col="red")
#abline(h=TRACKSPACERS[3], col="blue")
## SNPS: 1-2
rect(xleft=AGDP_593$marker,ybottom=TRACKSPACERS[4]+0.0015, ytop=TRACKSPACERS[3]-0.0015, xright=AGDP501$marker+100)
abline(h=-1, lty=2)



rect(xleft=k1240$marker,ybottom=TRACKSPACERS[3]+0.0015, ytop=TRACKSPACERS[2]-0.0015, xright=k1240$marker+100)
abline(h=-1, lty=2)
## GWAS: 2-3
## GWAS: 2-3
gwascol<-colorRampPalette(colors = c('red','blue'))

#genes
glistdiff<-TRACKSPACERS[1]-TRACKSPACERS[2]
glistjitter<-6
glistbool.tmp<-rep_len(TRACKSPACERS[1], length.out  = glistjitter-2)
tmp2<-c(1:(glistjitter-2))
glistbool.tmp2<-glistbool.tmp-(tmp2*(glistdiff/glistjitter))


glist3$BOOL<-rep_len(x = glistbool.tmp2,length.out = nrow(glist3))
#abline(h=unique(glist$BOOL), lty=2, col="red")



apply(glist3, MARGIN = 1, function(x){
  if(x[1]==6){
    rect(xleft=x[2], xright=x[3], ybottom=(as.numeric(x[5])-(glistdiff/glistjitter*2)), ytop=(as.numeric(x[5])-(2*(glistdiff/(glistjitter*2)))), col = "gray50", border = "gray37")
    text(x = ((as.numeric(x[2])+as.numeric(x[3]))/2), y=as.numeric(x[5]), labels = as.character(x[4]), cex=0.4, col = "grey5")
  }
})

par(new=TRUE)



delimit_MHC$marker<-delimit_MHC[,2]


delimit_MHC_pos<-delimit_MHC[which(delimit_MHC$HGall >= 0),]
plot(x=delimit_MHC_pos$marker, y=delimit_MHC_pos$HGallz, ylim=c(-1.6,0.7), xlim=c(28000000,34000000), type='b', axes=F, ylab=NA, xlab=NA, col=(LAcol(10)[(floor((delimit_MHC_pos$HGall))+6)]), pch=20) 

delimit_MHC_neg<-delimit_MHC[which(delimit_MHC$HGall < 0),]
points(x=delimit_MHC_neg$marker, y=delimit_MHC_neg$HGallz, pch=20, col=(LAcol(10)[(ceiling((delimit_MHC_neg$HGall))+6)]))


plotmax=0.9
LAcol<-colorRampPalette(colors = c('red','grey','royalblue4'))
color.legend(
  xl=28000000,
  xr=29500000,
  yb = -1.3,
  yt = -1.1 ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

text(x=(33950000+34190000)/2, y=(plotmax), labels = "Z-score")



axis(4, at = c(0,0.2,0.4,0.6,0.8,1))
mtext(side = 4, text = "% Neolithic Ancestry", line=2, adj=0.875)


dev.off()



## HERC tiff


#pdf("1240kvsAGDPsnpwindows_pruned2.pdf",width=12, height = 8)

AGDP_593<-AGDP_593_SLC
delimit_MHC<-LAdf15


tiff(filename="1240kvsAGDPsnpwindows_pruned2_HERC2.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
par(mar=c(5, 4, 4, 4) + 0.1)
m=28356182
plot(AGDP_593$datawg, x =AGDP_593$marker, ylim =c(-0.15,0.3), xlim=c(m-3e+06,m+3e+06), xlab = "Marker (Chr 15)", ylab="Mixed-weight statistic", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HERC2 signal", type='b', lwd=5, col="grey50")


#plot(k1240$datawg, x =k1240$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight raw value", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")

points(k1240$datawg, x = k1240$marker, col = "black", type='l', lwd=3, cex=0.5)
points(k1240$datawg, x = k1240$marker, col = "black", type='p', lwd=3, cex=0.5)

legend(x="topleft", legend=c("1240k Array - 21 snp window","Whole-genome shotgun - 551 snp window"), col = c("black","grey"), lty=c(1), lwd=3, cex=0.8)


axis(1)
axis(2, at=c(0,0.05,0.1,0.15,0.2,0.25,0.3), cex.axis=0.7)

#dev.off()

axis(2, at=c(-0.06125, -0.1125,-(0.15+0.125)/2), labels = c("GENES","1240k SNPs","Shotgun SNPs"), cex.axis=0.5, las=1)

abline(h=scale_limits[2], lty=1) #genes upper
abline(h=-0.1, lty=1) # genes lower
abline(h=-0.125) # 1240k
abline(h=-0.15) #

TRACKSPACERS=c(scale_limits[2],-0.1,-0.125,-0.15)
#abline(h=TRACKSPACERS[1], col="orange")
#abline(h=TRACKSPACERS[2], col="green")
#abline(h=TRACKSPACERS[4], col="red")
#abline(h=TRACKSPACERS[3], col="blue")
## SNPS: 1-2
rect(xleft=AGDP_593$marker,ybottom=TRACKSPACERS[4]+0.0015, ytop=TRACKSPACERS[3]-0.0015, xright=AGDP501$marker+100)
abline(h=-1, lty=2)



rect(xleft=k1240$marker,ybottom=TRACKSPACERS[3]+0.0015, ytop=TRACKSPACERS[2]-0.0015, xright=k1240$marker+100)
abline(h=-1, lty=2)
## GWAS: 2-3
## GWAS: 2-3
gwascol<-colorRampPalette(colors = c('red','blue'))

#genes
glistdiff<-TRACKSPACERS[1]-TRACKSPACERS[2]
glistjitter<-6
glistbool.tmp<-rep_len(TRACKSPACERS[1], length.out  = glistjitter-2)
tmp2<-c(1:(glistjitter-2))
glistbool.tmp2<-glistbool.tmp-(tmp2*(glistdiff/glistjitter))


glist15$BOOL<-rep_len(x = glistbool.tmp2,length.out = nrow(glist15))
#abline(h=unique(glist$BOOL), lty=2, col="red")



apply(glist15, MARGIN = 1, function(x){
  if(x[1]==15){
    rect(xleft=x[2], xright=x[3], ybottom=(as.numeric(x[5])-(glistdiff/glistjitter*2)), ytop=(as.numeric(x[5])-(2*(glistdiff/(glistjitter*2)))), col = "gray50", border = "gray37")
    text(x = ((as.numeric(x[2])+as.numeric(x[3]))/2), y=as.numeric(x[5]), labels = as.character(x[4]), cex=0.4, col = "grey5")
  }
})

par(new=TRUE)



delimit_marker<-delimit_MHC[,2]
delimit_MHC$marker<-delimit_MHC$V2


delimit_MHC_pos<-delimit_MHC[which(delimit_MHC$HGall >= 0),]
plot(x=delimit_MHC_pos$marker, y=delimit_MHC_pos$HGall, ylim=c(-1.6,0.9), xlim=c(m-3e+06,m+3e+06), type='b', axes=F, ylab=NA, xlab=NA, col=(LAcol(10)[(floor((delimit_MHC_pos$HGallZ))+6)]), pch=20) 

delimit_MHC_neg<-delimit_MHC[which(delimit_MHC$HGall < 0),]
points(x=delimit_MHC_neg$marker, y=delimit_MHC_neg$HGall, pch=20, col=(LAcol(10)[(ceiling((delimit_MHC_neg$HGallZ))+6)]))


plotmax=0.9
LAcol<-colorRampPalette(colors = c('red','grey','royalblue4'))
color.legend(
  xl=3.1e+07,
  xr=m+3e+06,
  yb = 0.2,
  yt = 0.85,
  gradient = "y" ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

text(x=(33950000+34190000)/2, y=(plotmax), labels = "Z-score")



axis(4, at = c(0,0.2,0.4,0.6,0.8,1))
mtext(side = 4, text = "% Neolithic Ancestry", line=2, adj=0.875)


dev.off()

## round negative and positive sepe
tiff(filename="1240kvsAGDPsnpwindows_pruned2.tiff", height= 8, units = "in", width = 14, res = 150, compression="lzw")
par(mar=c(5, 4, 4, 4) + 0.1)
plot(AGDP_593$datawg, x =AGDP_593$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight statistic", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")


#plot(k1240$datawg, x =k1240$marker, ylim =c(-0.15,0.3), xlim=c(28000000,34000000), xlab = "Marker (Chr 6)", ylab="Mixed-weight raw value", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal", type='b', lwd=5, col="grey50")

points(k1240$datawg, x = k1240$marker, col = "black", type='l', lwd=3, cex=0.5)
points(k1240$datawg, x = k1240$marker, col = "black", type='p', lwd=3, cex=0.5)

legend(x="topleft", legend=c("1240k Array - 21 snp window","Whole-genome shotgun - 551 snp window"), col = c("black","grey"), lty=c(1), lwd=3, cex=0.8)


axis(1)
axis(2, at=c(0,0.05,0.1,0.15,0.2,0.25,0.3), cex.axis=0.7)

#dev.off()

axis(2, at=c(-0.06125, -0.1125,-(0.15+0.125)/2), labels = c("GENES","1240k SNPs","Shotgun SNPs"), cex.axis=0.5, las=1)

abline(h=scale_limits[2], lty=1) #genes upper
abline(h=-0.1, lty=1) # genes lower
abline(h=-0.125) # 1240k
abline(h=-0.15) #

TRACKSPACERS=c(scale_limits[2],-0.1,-0.125,-0.15)
#abline(h=TRACKSPACERS[1], col="orange")
#abline(h=TRACKSPACERS[2], col="green")
#abline(h=TRACKSPACERS[4], col="red")
#abline(h=TRACKSPACERS[3], col="blue")
## SNPS: 1-2
rect(xleft=AGDP_593$marker,ybottom=TRACKSPACERS[4]+0.0015, ytop=TRACKSPACERS[3]-0.0015, xright=AGDP501$marker+100)
abline(h=-1, lty=2)



rect(xleft=k1240$marker,ybottom=TRACKSPACERS[3]+0.0015, ytop=TRACKSPACERS[2]-0.0015, xright=k1240$marker+100)
abline(h=-1, lty=2)
## GWAS: 2-3
## GWAS: 2-3
gwascol<-colorRampPalette(colors = c('red','blue'))

#genes
glistdiff<-TRACKSPACERS[1]-TRACKSPACERS[2]
glistjitter<-6
glistbool.tmp<-rep_len(TRACKSPACERS[1], length.out  = glistjitter-2)
tmp2<-c(1:(glistjitter-2))
glistbool.tmp2<-glistbool.tmp-(tmp2*(glistdiff/glistjitter))


glist3$BOOL<-rep_len(x = glistbool.tmp2,length.out = nrow(glist3))
#abline(h=unique(glist$BOOL), lty=2, col="red")



apply(glist3, MARGIN = 1, function(x){
  if(x[1]==6){
    rect(xleft=x[2], xright=x[3], ybottom=(as.numeric(x[5])-(glistdiff/glistjitter*2)), ytop=(as.numeric(x[5])-(2*(glistdiff/(glistjitter*2)))), col = "gray50", border = "gray37")
    text(x = ((as.numeric(x[2])+as.numeric(x[3]))/2), y=as.numeric(x[5]), labels = as.character(x[4]), cex=0.4, col = "grey5")
  }
})

par(new=TRUE)



delimit_MHC$marker<-delimit_MHC[,2]


delimit_MHC_pos<-delimit_MHC[which(delimit_MHC$HGall >= 0),]
plot(x=delimit_MHC_pos$marker, y=delimit_MHC_pos$HGallz, ylim=c(-1.6,0.9), xlim=c(28000000,34000000), type='b', axes=F, ylab=NA, xlab=NA, col=(LAcol(10)[(floor((delimit_MHC_pos$HGall))+6)]), pch=20) 

delimit_MHC_neg<-delimit_MHC[which(delimit_MHC$HGall < 0),]
points(x=delimit_MHC_neg$marker, y=delimit_MHC_neg$HGallz, pch=20, col=(LAcol(10)[(ceiling((delimit_MHC_neg$HGall))+6)]))


plotmax=0.9
LAcol<-colorRampPalette(colors = c('red','grey','royalblue4'))
color.legend(
  xl=33950000,
  xr=34190000,
  yb = plotmax - 0.8*plotmax,
  yt = plotmax - 0.05*plotmax, gradient = "y" ,
  legend = seq(from= -4, to = 4, by=1),
  rect.col = LAcol(9), cex = 0.6)

text(x=(33950000+34190000)/2, y=(plotmax), labels = "Z-score")



axis(4, at = c(0,0.2,0.4,0.6,0.8,1))
mtext(side = 4, text = "% Neolithic Ancestry", line=2, adj=0.875)


dev.off()


#####

delimit_MHC<-LAdf6
colnames(delimit_MHC)[c(1,2)]<-c("row.names","marker")
pdf(file="LA_allchroms.tmp.pdf", width=20, height=40)
par(mfrow=c(11,2))
for(x in 1:22){
  LA<-LAdf[which(LAdf$V1==x),]
  #plot(x=LA$V2, y=LA$HGall, ylim=c(0,0.8),type='p', col=(LAcol(10)[(floor((LA$HGallZ))+6)]), pch=20, main = paste0("Chromosome ",x), xlab="Marker",ylab="HG Ancestry (%)")
  
  LA_pos<-LA[which(LA$HGallZ >=0),]
  LA_neg<-LA[which(LA$HGallZ < 0),]
  
  plot(x=LA_pos$V2, xlim=c(min(LA$V2), max(LA$V2)), y=LA_pos$HGall, ylim=c(0,0.8), type='p', col=(LAcol(10)[(floor((LA_pos$HGallZ))+6)]), pch=20, main = paste0("Chromosome ",x), xlab="Marker",ylab="Mesolithic Ancestry (%)", cex=2)
  
  
  points(x=LA_neg$V2, xlim=c(0, max(LA$V2)), y=LA_neg$HGall, ylim=c(0,0.8), type='p', col=(LAcol(10)[(ceiling((LA_neg$HGallZ))+6)]), pch=20, cex=2)
}
dev.off()



# Whole genome Z-plotting


tiff(filename="localancestry_wholegenome_z_grey2_panel2.tiff", res = 200, width = 1500, height = 1000)


layout(matrix(c(1,1,1,1,2
), nrow=1, byrow=TRUE), heights = c(1))

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


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$HGallZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.6), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$HGallZ))+5)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox))) 





color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.95*ymax,
  yt = ymax - 0.85*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)



df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$HGallZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$HGallZ))+5)]))

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


HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = df_filt_markerprox$HGallZ, bed = glist0, outlier_threshold = 4, plot_filter=0) 
# manual edit to remove clipping names ]
HG_outliers<-HG_outliers[-11,]

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$HGallZ), bed = glist0, outlier_threshold = 4, plot_filter=0) 

topoutliers<-HG_outliers



alloutliers<-cbind(HG_outliers, NEO_outliers)
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
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  if(x == 4){
    text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  }else{
    # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
    text(x = xpoints, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
    #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
  }}


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
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
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
  yb = ymax - 0.95*ymax,
  yt = ymax - 0.85*ymax, gradient = "x",
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=(xmax-0.9*xmax)+(xmax-0.85*xmax), y=(ymax - 0.83*ymax), labels = "Z-score")


### adding histogram in next panel

#(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.65), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol(10)[(floor((df_filt_markerprox_pos$HGallZ))+6)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox+(0.4*df_filt_markerprox$markerprox)))) 

par(mar=c(5, 1, 4, 2) + 0.1)


A<-hist(df_filt_markerprox$HGall, breaks = 40, plot = F)


#rect(0, A$breaks[1:(length(A$breaks) - 1)], A$counts, A$breaks[2:length(A$breaks)])



horiz.hist <- function(Data, breaks=60, col="transparent", las=1, ylim=range(HBreaks), labelat=pretty(ylim), labels=labelat, border=par("fg"), ... )
{a <- hist(Data, plot=FALSE, breaks=breaks)
HBreaks <- a$breaks
HBreak1 <- a$breaks[1]
hpos <<- function(Pos) (Pos-HBreak1)*(length(HBreaks)-1)/ diff(range(HBreaks))
barplot(a$counts, space=0, horiz=T, ylim=hpos(ylim), col=LAcol(3)[2], border="black",axes=F, )      
#axis(2, at=hpos(labelat), labels=labels, las=las, ...) 
print("use hpos() to address y-coordinates") }

horiz.hist(df_filt_markerprox$HGall, ylim = c(0,0.6))



#plot.window(ylim=c(0.0,0.65), xlim=c(0,1))
#xhist <- hist(df_filt_markerprox$HGall, breaks = 20, plot = FALSE)
#barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", ylim = c(0.0,0.65))


dev.off()

####

# delta Z score
tiff(filename="localancestry_wholegenome_z_delta_grey2.tiff", res = 200, width = 1500, height = 1000)


layout(matrix(c(1,1,1,1,2
), nrow=1, byrow=TRUE), heights = c(1))

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


df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$HGZcomp >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGZcomp, ylim=c(-2.5,2), type='p', axes=F, ylab="delta Z-score", xlab="Chromosome", col=(LAcol_grey2[(floor((df_filt_markerprox_pos$HGZcomp))+5)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox))) 





df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$HGZcomp < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGZcomp, pch=20, col=(LAcol_grey2[(ceiling((df_filt_markerprox_neg$HGZcomp))+5)]))

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
dev.off()

HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = df_filt_markerprox$HGZcomp, bed = glist0, outlier_threshold = 0.65*max(LAdf$HGZcomp), plot_filter=0) 
# manual edit to remove clipping names ]
#HG_outliers<-HG_outliers[-11,]

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$HGZcomp), bed = glist0, outlier_threshold = 0.99*min(-1* df_filt_markerprox$HGZcomp), plot_filter=0) 

topoutliers<-HG_outliers



alloutliers<-cbind(HG_outliers, NEO_outliers)
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
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  if(x == 4){
    text(x = xpoints-250000000, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  }else{
    # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
    text(x = xpoints, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
    #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
  }}


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
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
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
  yb = ymax - 0.95*ymax,
  yt = ymax - 0.85*ymax, gradient = "x" ,
  legend = seq(from= -4, to = 4, by=1),
  align = "rb",
  rect.col = LAcol_grey2, cex = 0.85)

text(x=(xmax-0.9*xmax)+(xmax-0.85*xmax), y=(ymax - 0.83*ymax), labels = "Z-score")


### adding histogram in next panel

#(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.65), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol(10)[(floor((df_filt_markerprox_pos$HGallZ))+6)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic", xlim=c(0,max(df_filt_markerprox$markerprox+(0.4*df_filt_markerprox$markerprox)))) 

par(mar=c(5, 1, 4, 2) + 0.1)


A<-hist(df_filt_markerprox$HGall, breaks = 40, plot = F)


#rect(0, A$breaks[1:(length(A$breaks) - 1)], A$counts, A$breaks[2:length(A$breaks)])



horiz.hist <- function(Data, breaks=60, col="transparent", las=1, ylim=range(HBreaks), labelat=pretty(ylim), labels=labelat, border=par("fg"), ... )
{a <- hist(Data, plot=FALSE, breaks=breaks)
HBreaks <- a$breaks
HBreak1 <- a$breaks[1]
hpos <<- function(Pos) (Pos-HBreak1)*(length(HBreaks)-1)/ diff(range(HBreaks))
barplot(a$counts, space=0, horiz=T, ylim=hpos(ylim), col=LAcol(3)[2], border="black",axes=F, )      
#axis(2, at=hpos(labelat), labels=labels, las=las, ...) 
print("use hpos() to address y-coordinates") }

horiz.hist(df_filt_markerprox$HGall, ylim = c(0,0.6))



#plot.window(ylim=c(0.0,0.65), xlim=c(0,1))
#xhist <- hist(df_filt_markerprox$HGall, breaks = 20, plot = FALSE)
#barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", ylim = c(0.0,0.65))


dev.off()


####
pdf(file="localancestry_wholegenome_z3.pdf", width=12, height=7)




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




df_filt_markerprox_pos<-df_filt_markerprox[which(df_filt_markerprox$HGallZ >= 0),]
plot(x=df_filt_markerprox_pos$markerprox, y=df_filt_markerprox_pos$HGall, ylim=c(0.0,0.65), type='p', axes=F, ylab="Hunter-Gatherer Ancestry (%)", xlab="Chromosome", col=(LAcol(10)[(floor((df_filt_markerprox_pos$HGallZ))+6)]), pch=20, main = "Local ancestry in the admixed middle-Neolithic") 

df_filt_markerprox_neg<-df_filt_markerprox[which(df_filt_markerprox$HGallZ < 0),]
points(x=df_filt_markerprox_neg$markerprox, y=df_filt_markerprox_neg$HGall, pch=20, col=(LAcol(10)[(ceiling((df_filt_markerprox_neg$HGallZ))+6)]))



axis(2)
chrindex=NULL
for (x in 1:22){
  chrindex[x]<-tail(df_filt_markerprox[(which(df_filt_markerprox$row.names==x)),]$markerprox, 1)
}
chrindex2<-append(0, chrindex)
chrmidpoint<-chrindex2[-length(chrindex2)] + diff(chrindex2)/2


axis(1, at = chrmidpoint, labels=c(1:22), tick = TRUE, cex.axis=1)



HG_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = df_filt_markerprox$HGallZ, bed = glist0, outlier_threshold = 3, plot_filter=0) 

NEO_outliers<-return_outliers(df = df_filt_markerprox, genomic_values = (-1* df_filt_markerprox$HGallZ), bed = glist0, outlier_threshold =3, plot_filter=0) 

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
  # points(x = xpoints, y = ypoints+0.03, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints+0.03, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}

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
  ypoints<-df_filt_markerprox[topoutliers_unique[x,]$uid,]$HGall
  
  
  #points(x=xpoints,y=ypoints, pch=4, col="red")}
  coli<-x
  if(coli%%15 == 0){
    coli<-coli+1
  }
  #points(x = xpoints, y = ypoints-0.045, pty=20, cex=0.7, pch = 16, col = "green")
  text(x = xpoints, y = ypoints-0.045, cex=0.7, col = "black", labels = topoutliers_unique[x,]$gene)
  #points(x = xpoints, y = ypoints, col=colfunc(15)[coli%%15], pty=20, cex=0.7)
}




color.legend(
  xl=xmax-0.9*xmax,
  xr=xmax-0.5*xmax,
  yb = ymax - 0.95*ymax,
  yt = ymax - 0.85*ymax, gradient = "x" ,
  legend = seq(from= 0, to = 10, by=1)-5,
  rect.col = LAcol(11), cex = 1)


text(x=(xmax-0.9*xmax)+(xmax-0.8*xmax), y=(ymax - 0.92*ymax), labels = "Z-score")


dev.off()


# abline(h=mean(LAdf$HGall), col="grey", lty=5)
# abline(h=mean(LAdf$HGall)-(1*sd(LAdf$HGall)), col="red", lty=5)
# abline(h=mean(LAdf$HGall)-(2*sd(LAdf$HGall)), col="red1", lty=5)
# abline(h=mean(LAdf$HGall)-(3*sd(LAdf$HGall)), col="red2", lty=5)
# abline(h=mean(LAdf$HGall)-(4*sd(LAdf$HGall)), col="red3", lty=5)
# 
# abline(h=mean(LAdf$HGall)+(1*sd(LAdf$HGall)), col="royalblue", lty=5)
# abline(h=mean(LAdf$HGall)+(2*sd(LAdf$HGall)), col="royalblue1", lty=5)
# abline(h=mean(LAdf$HGall)+(3*sd(LAdf$HGall)), col="royalblue2", lty=5)
# abline(h=mean(LAdf$HGall)+(4*sd(LAdf$HGall)), col="royalblue3", lty=5)
# abline(h=mean(LAdf$HGall)+(5*sd(LAdf$HGall)), col="royalblue3", lty=5)
# abline(h=mean(LAdf$HGall)+(6*sd(LAdf$HGall)), col="royalblue3", lty=5)

#axis(4, at=c(0,0.1,0.2,0.3,0.4))
axis(4, at=seq(from=0, to=0.7, by=0.1))
#axis(4,at=seq(from=0.6, to=1.2, by=0.1))
#axis(4, at=c(0,10,20,30,40,50,60,70,80,90,100))
mtext(side = 4, text = "Hunter-Gatherer Ancestry", line = 2.4, las = 3)
plot(y=LAdf6$HG.ancestry ,x=LAdf6$V2, ylim=c(0,1), xlim=c(32006018,33199864))



legend("topright", col=gwascol(5), pch=19)

pvallegend<-unique(cut(gwasplot$pval, breaks = 5))
pvallegend_manual<-c(-.94,2.71,3.48,4.24,5.01,5.78)
pvallegend_manual_pval<-10^-exp(gwasplot$pval)
x<-sort(10^-exp(gwasplot$pval))
abline(h=0.02)
abline(h=-log10(5e-8), lty=2, col = "red")






plot(x=AGDP501$marker, AGDP501$GT2L, ylim=c(-40, 150), type='p',axes=F, col = "cadetblue", ylab = NA, xlab=NA, lwd=2, pch=20)
#points(x=fullMHC$marker, fullMHC$GT2L, type='p', col = "darkorchid", lwd=2,pch=20)
#points(x=fullMHC$marker, fullMHC$GT3L, type='p',col = "cadetblue",lwd=2, pch=20)
#axis(4, at=c(0,0.1,0.2,0.3,0.4)) #axis(4,at=seq(from=0.6, to=1.2, by=0.1)) axis(4, at=c(0,10,20,30,40,50,60,70,80,90,100))
mtext(side =4, text = "SNP Coverage", line = 2.4, las = 3)
legend(x="topright", col=c("coral3","cadetblue","darkorchid"), legend = c("HG","NEO","MNEO"), title= "Coverage", pch =20,bg="white")
#abline(h=20, col="black", lty=1, lwd =2)
axis(side=4, at=seq(from=0,to=100,by=10))






##### swapping SNP coverage for 1240k MW raw value. . .

plot(AGDP501$datawg, x =AGDP501$marker, ylim =c(-0.03,0.25), xlab = "Marker (Chr 6)", ylab="Mixed-weight raw value", axes=F, pch=20, main = "Zoom-in on mixed-weight agdp HLA signal")
axis(1)
axis(2, at=c(0,0.05,0.1,0.15,0.2), cex.axis=0.7)
axis(2, at=c(-0.05,-0.1, -0.15), labels = c("GENES","GWAS","SNPs"), cex.axis=0.5, las=1)

TRACKSPACERS=c(0.02,0.005,-0.01,scale_limits[2])
TRACKSPACER=0.0015
GENEHEIGHT=0.3*(abs(TRACKSPACERS[4])-abs(TRACKSPACERS[3]))

abline(h=TRACKSPACERS, lty=2)
## SNPS: 1-2
rect(xleft=AGDP501$marker,ybottom=TRACKSPACERS[2]+0.0015, ytop=TRACKSPACERS[1]-0.0015, xright=AGDP501$marker+100)
abline(h=-1, lty=2)
## GWAS: 2-3
gwascol<-colorRampPalette(colors = c('red','blue'))
#gwasplot$pvalog<-(gwasplot$pval)
#gwasplot$col <- cut(gwasplot$pvalog, breaks = 100, )

gwasplot$col<-findInterval(gwasplot$beta, seq(from=min(gwasplot$beta), to=max(gwasplot$beta), length.out=10))
rect(xleft=gwasplot$pos,ybottom=TRACKSPACERS[3], ytop=TRACKSPACERS[2], xright=gwasplot$pos+2000, col = gwascol(10)[ceiling(gwasplot$col)] ,border=NA)

#genes
glist$BOOL<-rep_len(x=c(TRACKSPACERS[3]-((abs(TRACKSPACERS[4])-abs(TRACKSPACERS[3]))*1/3), (TRACKSPACERS[3]-((abs(TRACKSPACERS[4])-abs(TRACKSPACERS[3]*2/3))))), length.out = nrow(glist))
#abline(h=unique(glist$BOOL), lty=2, col="red")



apply(glist, MARGIN = 1, function(x){
  if(x[1]==6){
    rect(xleft=x[2], xright=x[3], ybottom=as.numeric(x[5])-GENEHEIGHT, ytop=as.numeric(x[5]), col = "gray50", border = "gray37")
    text(x = ((as.numeric(x[2])+as.numeric(x[3]))/2), y=as.numeric(x[5])-TRACKSPACER, labels = as.character(x[4]), cex=0.62, col = "grey5")
  }
})
legend("topright", col=gwascol(5), pch=19)

pvallegend<-unique(cut(gwasplot$pval, breaks = 5))
pvallegend_manual<-c(-.94,2.71,3.48,4.24,5.01,5.78)
pvallegend_manual_pval<-10^-exp(gwasplot$pval)
x<-sort(10^-exp(gwasplot$pval))
abline(h=0.02)
abline(h=-log10(5e-8), lty=2, col = "red")
par(new=TRUE)

plot(y=k1240$datawg, x=k1240$marker, ylim=c(-0.03,0.25), type='p',axes=F, col = "orange", ylab = NA, xlab=NA, lwd=2, pch=20, xlim=c(32000000,33200000))
#points(x=fullMHC$marker, fullMHC$GT2L, type='p', col = "darkorchid", lwd=2,pch=20)
#points(x=fullMHC$marker, fullMHC$GT3L, type='p',col = "cadetblue",lwd=2, pch=20)
#axis(4, at=c(0,0.1,0.2,0.3,0.4)) #axis(4,at=seq(from=0.6, to=1.2, by=0.1)) axis(4, at=c(0,10,20,30,40,50,60,70,80,90,100))
mtext(side =4, text = "1240k Mixed-Weight raw value", line = 2.4, las = 3)
#abline(h=20, col="black", lty=1, lwd =2)
axis(side=4, at=seq(from=0,to=100,by=10))









