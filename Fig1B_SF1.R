library(maps)
setwd("~/EvoFunc/ms_scripts")
source("~/EvoFunc/code/functionsSelSweep.R")

metadata<-read.table("/Users/davyt/slurmp/rebuttal/redo/all_meta_cols_pops.tsv")
colnames(metadata)<-c("GROUP","IND","AGE","LAT","LONG","COV","POP")

df_split<-split(metadata, f = metadata$POP) # Split by Chr
z=NULL
y=0

lapply(df_split, function(x){
  x[,3:6]<-as.numeric(x[,3:6])
})

#HG minmax

tiff(filename="davyetal_2022_fig1B.tiff", res=250, compression = "lzw", height=2000, width=2000)

map(xlim=c(-12,40), ylim=c(35,70), resolution = 0, bg = "white", fill = T, col="grey90")


grey2<-add.alpha("grey50",0.3)
df<-do.call("rbind", df_split)
df<-df[order(df$LAT,df$POP,decreasing = c(T,T)),]
cols<-as.vector(NULL)
HG.col<-"#412CA8"
NEO.col<-"#e40959"
MNEO.col<-"goldenrod2"
cols[which(df$POP=="HG")]<-HG.col
cols[which(df$POP=="NEO")]<-NEO.col
cols[which(df$POP=="MNEO")]<-MNEO.col
df$cols<-cols
cols<-as.vector(NULL)

### SHADOW

rm(x.store)
multiplier=0
for(x in 1:nrow(df)){
  x.current<-as.numeric(df[x,4:5])
  x.col<-df[x,]$cols
  
  
  if(!exists("x.store")){
    x.store<-as.numeric(df[x,4:5])
    x.shadow<-x.store+c(-0.1, 0.1)
    points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    if(multiplier > 2){
      x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
      points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20 ,lwd=0)}

    } else {
    x.shadow<-x.current+c(-0.1, 0.1)
    points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    multiplier=0
  }
  
  x.store<-x.current
  
}



##POLE

rm(x.store)
multiplier=0


for(x in 1:nrow(df)){
  x.current<-as.numeric(df[x,4:5])
  x.col<-df[x,]$cols
  
  
  
  if(!exists("x.store")){
    x.store<-as.numeric(df[x,4:5])
    points(x=x.store[2], y=x.store[1], col = df[x,]$cols,pch=20, cex=1)
    points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    x.plot<-x.current+(multiplier*c(0.15,0))
    x.ring<-x.plot
    points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
    points(x=x.plot[2], y=x.plot[1], col = df[x,]$cols, pch=20)
    
    
  } else {
    multiplier=0
    points(x=x.current[2], y=x.current[1], col = df[x,]$cols, pch=20)
    points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
  }
  
  x.store<-x.current
  
}


df_blu<-df[df$cols=="blue",]

legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c(NEO.col,MNEO.col,HG.col), pch=20, cex=1.5)
legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c("black"), pch=c(1), bg = 'transparent', pt.cex = 0.9, cex=1.5)



df_NEO<-df[df$POP=="NEO_NEW",]
df_blu<-df[df$col=="blue",]



dev.off()


tiff(filename="davy2022etal_figS1B.tiff",width = 1000, height=1000, res=150, compression = "lzw")
plot.new()
box()
plot.window(xlim = c(0, 4), ylim = c(min(df$AGE)-1000,max(df$AGE)+1000))
grid(nx=F,ny=NULL)
par(mar=c(5, 4, 4, 2) + 0.1)
boxplot(df_split$HG$AGE, df_split$NEO$AGE, df_split$MNEO$AGE, col=c(add.alpha(HG.col, 1),add.alpha(NEO.col,1),add.alpha(MNEO.col,1)),yaxt=F,add=T, lwd=1.5, pch=21)
axis(1, labels = c("Mesolithic","Neolithic","Admixed Neolithic"), at =c(1,2,3))
axis(2)
title(ylab = "Time (kya)")
dev.off()

tiff(filename="davy2022etal_figS1A.tiff",width = 1000, height=1000, res=150, compression = "lzw")

par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
plot.window(xlim = c(0, 4), ylim = c(min(df$COV),  1240000))
grid(nx=F,ny=NULL)
boxplot(df_split$HG$COV, df_split$NEO$COV, df_split$MNEO$COV, col=c(add.alpha(HG.col, 1),add.alpha(NEO.col,1),add.alpha(MNEO.col,1)), add=T, lwd=1.5, pch=21)
axis(1, labels = c("Mesolithic","Neolithic","Admixed Neolithic"), at =c(1,2,3))
axis(2)
title(ylab = "Coverage (Number of SNPs)")
dev.off()





