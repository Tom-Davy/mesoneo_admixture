#### ---- Prereq------
options(scipen=-30)

source("~/EvoFunc/code/functionsSelSweep.R")x
library(fitdistrplus)
library(data.table)
library(scales)
library(stringr)


#### Parameters
snp<-51
row.names<-22
nullwind<-5000000

#### Carp
df<-as.data.frame(fread("/Users/davyt/slurmp/rebuttal/redo/REF_MNEO_HG_NEO_freqstats_mincount1_intersect.tsv"))
colnames(df)<-c("row.names","marker1","marker","GT1","GT2","GT3","GT4","p1","p2","p3","p4","datawg","rs.number")
df$df.uid<-seq_along(df[,1])

## F2s
df$F2_NEOHG<-(df$p3-df$p4)^2
df$F2_MNEOHG<-(df$p2-df$p3)^2
df$F2_NEOMNEO<-(df$p2-df$p4)^2


#### GT lengths
df$GT1L<-str_length(df$GT1)/2
df$GT2L<-str_length(df$GT2)/2
df$GT3L<-str_length(df$GT3)/2
df$GT4L<-str_length(df$GT4)/2

# Mincount Hist
GT1<-hist(str_length(df$GT1)/2, breaks = 100, plot = FALSE)
GT2<-hist(str_length(df$GT2)/2, breaks = 100, plot = FALSE)
GT3<-hist(str_length(df$GT3)/2, breaks = 100, plot = FALSE)
GT4<-hist(str_length(df$GT4)/2, breaks = 50, plot = FALSE)



# Filtering on >20 NEO
df<-df[which(df$GT4L >= 20),]
df_filt<-df
datawg<-df_filt$datawg



#### Generating Null Distribution -------

#### Generating Null Windows



nullsnps<-windowsampler(nullwind, df_filt, row.names, snp)
datanull<-get_null(datawg, snp, nullsnps)

#Generating per-statistic nulls
datanull.MN<-get_null(df_filt$F2_NEOMNEO, snp, nullsnps)
datanull.MH<-get_null(df_filt$F2_MNEOHG, snp, nullsnps)
datanull.HN<-get_null(df_filt$F2_NEOHG, snp, nullsnps)


setwd("~/EvoFunc/ms_scripts")

# Sliding Window -----
message("Generating Sliding Windows")
slidingwindow<-windowseeker(snp, df_filt, row.names)
if(snp > 1){
  df_filt$MN_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_NEOMNEO, snp, slidingwindow)), df_filt, 22)
  df_filt$MH_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_MNEOHG, snp, slidingwindow)), df_filt, 22)
  df_filt$HN_sliding<-slidingwindow_NAccountant_test(rowMeans(get_null(df_filt$F2_NEOHG, snp, slidingwindow)), df_filt, 22)
  df_filt_store<-df_filt
  df_filt=df_filt[-which(is.na(df_filt$HN_sliding)),]
}



##### Pval derivation ------
G.MN<-fitdist(datanull.MN, "gamma", method = "mle", keepdata = T)
G.MH<-fitdist(datanull.MH, "gamma", method = "mle", keepdata = T)
G.HN<-fitdist(datanull.HN, "gamma", method = "mle", keepdata = T)



pval.MN<--pgamma(df_filt$MN_sliding, shape=as.vector(G.MN$estimate[1]), rate = as.vector(G.MN$estimate[2]), lower.tail = F, log.p = T)/log(10)
pval.MH<--pgamma(df_filt$MH_sliding, shape=as.vector(G.MH$estimate[1]), rate = as.vector(G.MH$estimate[2]), lower.tail = F, log.p= T)/log(10)
pval.HN<--pgamma(df_filt$HN_sliding, shape=as.vector(G.HN$estimate[1]), rate = as.vector(G.HN$estimate[2]), lower.tail = F, log.p = T)/log(10)





### Plotting . . .

### Plotting prep/carp
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
#### Big Manhattan -----
options(scipen=999)

glist<-read.table("/Users/davyt/EvoFunc/glist.txt")


#### Classic ----
plot_filter<-0
plot_max<-25
df<-df_filt
genomic_values<-pval.MN
bed=glist
outlier_threshold<--log10(5e-8)
chromosomes<-22

topoutliers.HN<-return_outliers(df_filt, pval.HN, glist, outlier_threshold, plot_filter)
topoutliers.MH<-return_outliers(df_filt, pval.MH, glist, outlier_threshold, plot_filter)
topoutliers.MN<-return_outliers(df_filt, pval.MN, glist, outlier_threshold, plot_filter)
topoutliers=topoutliers.MN
#### Set a different plot_filter above for easier plotting. This handles that parameter.
df_smol<-df_filt[which(genomic_values > plot_filter),]
df_smol$uid<-seq_along(df_smol[,1])
GVsmol<-genomic_values[which(genomic_values > plot_filter)]
df_smol$smoluid<-seq_along(df_smol$uid)
df_smol$GVsmol<-GVsmol




tiff(filename="figS3A.tiff", units="in", width=14, height=7, res=150, compression = "lzw")


plot.new()
box()
plot.window(xlim = c(-5, max(df_smol$markerprox+10000)), ylim = c(plot_filter, plot_max))
#title(main="5M null distribution")
df_smol<-df_smol
abline(h= -log10(5e-8), col = "#e40959", lty = 2, lwd=2)
mbpeaks=NULL




col.HN<-"grey70"
col.MN<-"#bed357"
col.MH<-"#ee9c49ff"

points(df_filt$markerprox,pval.HN,col=col.HN, pch = 2)
points(df_filt$markerprox,pval.MH,col=col.MH, pch = 1)
points(df_filt$markerprox,pval.MN,col=col.MN, pch = 0)
#points(df_smol$markerpro

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


row.namesindex=NULL
for (x in 1:22){
  row.namesindex[x]<-tail(df_smol[(which(df_smol$row.names==x)),]$markerprox, 1)
}
row.namesindex2<-append(0, row.namesindex)
row.namesmidpoint<-row.namesindex2[-length(row.namesindex2)] + diff(row.namesindex2)/2

axis(1, at = row.namesmidpoint, labels=c(1:22), tick = FALSE, cex.axis=0.8)
axis(1, at = row.namesmidpoint, labels = FALSE, tick = TRUE)
axis(2, at=seq(from = 0, to = 50, by = 5))
title(xlab = "Chromosome", ylab = "P-value (-log10)", cex.main=0.8)
title(main="F2 comparisons, snp 51")



col.HN<-"grey70"
col.MN<-"#bed357"
col.MH<-"#ee9c49ff"

legend(x="topleft",legend=c("MNEO-NEO","MNEO-HG","HG-NEO"), col=c(col.MN, col.MH, col.HN), pch=c(0,1,2))

dev.off()


#png(filename = paste0(filename, ".QQ.png"), width=700, height=700, res=300, pointsize = 4)
tiff(filename = "figs3B.tiff", width=14, height=5, units = "in", res=150, compression = "lzw")

par(mfrow=c(1,3))



box()
shape = as.vector(head(G.MN$estimate))[1]
rate = as.vector(head(G.MN$estimate))[2]
data<-pgamma(df_filt$MN_sliding, shape = shape, rate=rate, lower.tail = F)
data<-data[data>0]
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- ceiling(max(logobspval))+3
expmax <- ceiling(max(logexppval))+3
plot(c(0,expmax), c(0,expmax), col="grey", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,24), las=1, xaxs="i", yaxs="i", bty="l")
box()
points(logexppval, logobspval, pch=20, col=col.MN, cex=1)

abline(a = 0, b = 1, col = "grey", lwd =4)
shape = as.vector(head(G.MN$estimate))[1]
rate = as.vector(head(G.MN$estimate))[2]
data<-pgamma(datanull.MN, shape = shape, rate=rate, lower.tail=F)
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
points(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)
#title(main= "5M observed vs Expected pvals - MN", xlab = "")
lines(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)





box()
shape = as.vector(head(G.HN$estimate))[1]
rate = as.vector(head(G.HN$estimate))[2]
data<-pgamma(df_filt$HN_sliding, shape = shape, rate=rate, lower.tail = F)
data<-data[data>0]
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- ceiling(max(logobspval))+3
expmax <- ceiling(max(logexppval))+3
plot(c(0,expmax), c(0,expmax), col="grey", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,24), las=1, xaxs="i", yaxs="i", bty="l")
box()
points(logexppval, logobspval, pch=20, col=col.HN, cex=1)

abline(a = 0, b = 1, col = "grey", lwd =4)
shape = as.vector(head(G.HN$estimate))[1]
rate = as.vector(head(G.HN$estimate))[2]
data<-pgamma(datanull.HN, shape = shape, rate=rate, lower.tail=F)
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
points(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)
#title(main= "5M observed vs Expected pvals", xlab = "")
lines(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)




box()
shape = as.vector(head(G.MH$estimate))[1]
rate = as.vector(head(G.MH$estimate))[2]
data<-pgamma(df_filt$MH_sliding, shape = shape, rate=rate, lower.tail = F)
data<-data[data>0]
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- ceiling(max(logobspval))+3
expmax <- ceiling(max(logexppval))+3
plot(c(0,expmax), c(0,expmax), col="grey", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,24), las=1, xaxs="i", yaxs="i", bty="l")
box()
points(logexppval, logobspval, pch=20, col=col.MH, cex=1)

abline(a = 0, b = 1, col = "grey", lwd =4)
shape = as.vector(head(G.MH$estimate))[1]
rate = as.vector(head(G.MH$estimate))[2]
data<-pgamma(datanull.MH, shape = shape, rate=rate, lower.tail=F)
obspval<-sort(data)
logobspval<--log(obspval)/log(10)
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
points(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)
#title(main= "5M observed vs Expected pvals", xlab = "")
lines(y = logobspval, x=logexppval, col = "black", pch=20, cex=1)



dev.off()