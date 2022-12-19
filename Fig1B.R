rm(list=ls())
library(rworldmap)
library(terra)
library(tmaptools)
library(maptools)
library(maps)
mapGriddedData(mapRegion = "europe", dataset = "", colourPalette = "white2Black")
setwd("~/EvoFunc/v50/")
source("~/EvoFunc/code/functionsSelSweep.R")

old=add.alpha("green", alpha=0.4)
new=add.alpha("blue", alpha=0.4)

metadata<-read.table("/Users/davyt/slurmp/rebuttal/redo/all_meta_cols_pops.tsv")

colnames(metadata)<-c("GROUP","IND","AGE","LAT","LONG","COV","POP")
#metadata[,3:6]<-as.numeric(metadata[,3:6])



df_split<-split(metadata, f = metadata$POP) # Split by Chr
z=NULL
y=0
for(x in 1:3){
  
  name<-names(df_split[x])
  z[[x]]<-cbind(name, apply(X = df_split[[x]], 2, function(x){sd(x)}))
}

names(z)<-unique(as.data.frame(do.call(rbind, z))[,1])
meta_sd<-z

lapply(df_split, function(x){
  x[,3:6]<-as.numeric(x[,3:6])
})
#longlat
# #HG minmax
# tiff(filename="v50_v2_comparison_longlat.tiff", res=250, compression = "lzw", height=2000, width=2000)
# #par(mfrow=c(1,4))
# #
# #plot.new()
# map(xlim=c(-12,40), ylim=c(35,70), col = "grey50", resolution = 0)
# 
# 
# grey2<-add.alpha("grey50",0.3)
# 
# 
# NEO_NEW<-df_split$NEO_NEW[order(df_split$NEO_NEW$LAT, decreasing = T),]
# #points(x=df_split$NEO_NEW$LONG, y=df_split$NEO_NEW$LAT, col=add.alpha("blue", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(NEO_NEW)){
#   x.current<-NEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-NEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     if(multiplier > 2){
#     x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)}
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black")
#     points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="blue", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# HG_NEW<-df_split$HG_NEW[order(df_split$HG_NEW$LAT, decreasing = T),]
# #points(x=df_split$HG_NEW$LONG, y=df_split$HG_NEW$LAT, col=add.alpha("green", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(HG_NEW)){
#   x.current<-HG_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-HG_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "green",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black")
#     points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="green", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# 
# MNEO_NEW<-df_split$MNEO_NEW[order(df_split$MNEO_NEW$LAT, decreasing = T),]
# #points(x=df_split$MNEO_NEW$LONG, y=df_split$MNEO_NEW$LAT, col=add.alpha("red", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(MNEO_NEW)){
#   x.current<-MNEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-MNEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "red",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black")
#     points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="red", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# dev.off()
# 
# 
# 
# 
# ##### Split shadow and pole
# 
# 
# 
# 
# #HG minmax
# tiff(filename="v50_v2_comparison_longlat_newshad.tiff", res=250, compression = "lzw", height=2000, width=2000)
# #par(mfrow=c(1,4))
# #
# #plot.new()
# map(xlim=c(-12,40), ylim=c(35,70), resolution = 0, bg = "white", fill = T, col="grey90")
# 
# 
# grey2<-add.alpha(grey2,0.3)
# 
# ### SHADOW
# NEO_NEW<-df_split$NEO_NEW[order(df_split$NEO_NEW$LAT, decreasing = T),]
# #points(x=df_split$NEO_NEW$LONG, y=df_split$NEO_NEW$LAT, col=add.alpha("blue", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(NEO_NEW)){
#   x.current<-NEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-NEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     #points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
#     #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     if(multiplier > 2){
#       x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#       points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20 ,lwd=0)}
#     #x.plot<-x.current+(multiplier*c(0.15,0))
#     #x.ring<-x.plot
#     #points(x=x.ring[2], y=x.ring[1], col = "black")
#     #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     #points(x=x.current[2], y=x.current[1], col="blue", pch=20)
#     multiplier=0
#     #points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# HG_NEW<-df_split$HG_NEW[order(df_split$HG_NEW$LAT, decreasing = T),]
# #points(x=df_split$HG_NEW$LONG, y=df_split$HG_NEW$LAT, col=add.alpha("green", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(HG_NEW)){
#   x.current<-HG_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-HG_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     #points(x=x.store[2], y=x.store[1], col = "green",pch=20, cex=1)
#     #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20,lwd=0)
#     #x.plot<-x.current+(multiplier*c(0.15,0))
#     #x.ring<-x.plot
#     #points(x=x.ring[2], y=x.ring[1], col = "black")
#     #points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     #points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     #points(x=x.current[2], y=x.current[1], col="green", pch=20)
#     multiplier=0
#     #points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# 
# MNEO_NEW<-df_split$MNEO_NEW[order(df_split$MNEO_NEW$LAT, decreasing = T),]
# #points(x=df_split$MNEO_NEW$LONG, y=df_split$MNEO_NEW$LAT, col=add.alpha("red", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(MNEO_NEW)){
#   x.current<-MNEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-MNEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     #points(x=x.store[2], y=x.store[1], col = "red",pch=20, cex=1)
#     #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20,lwd=0)
#     #x.plot<-x.current+(multiplier*c(0.15,0))
#     #x.ring<-x.plot
#     #points(x=x.ring[c], y=x.ring[1], col = "black")
#     #points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     #points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     multiplier=0
#     #points(x=x.current[2], y=x.current[1], col="red", pch=20)
#     #points(x=x.current[2], y=x.current[1], col="black")
#   }
#   
#   x.store<-x.current
#   
# }
# 
# ##POLE
# 
# rm(x.store)
# multiplier=0
# 
# 
# for(x in 1:nrow(NEO_NEW)){
#   x.current<-NEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-NEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     #if(multiplier > 2){
#     #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)}
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
#     points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="blue", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# HG_NEW<-df_split$HG_NEW[order(df_split$HG_NEW$LAT, decreasing = T),]
# #points(x=df_split$HG_NEW$LONG, y=df_split$HG_NEW$LAT, col=add.alpha("green", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(HG_NEW)){
#   x.current<-HG_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-HG_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "green",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black" ,cex=0.95)
#     points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="green", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="green", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
#   }
#   
#   x.store<-x.current
#   
# }
# 
# 
# 
# MNEO_NEW<-df_split$MNEO_NEW[order(df_split$MNEO_NEW$LAT, decreasing = T),]
# #points(x=df_split$MNEO_NEW$LONG, y=df_split$MNEO_NEW$LAT, col=add.alpha("red", 0.3), pch=20, cex=4)
# rm(x.store)
# multiplier=0
# for(x in 1:nrow(MNEO_NEW)){
#   x.current<-MNEO_NEW[x,4:5]
#   
#   
#   
#   if(!exists("x.store")){
#     x.store<-MNEO_NEW[x,4:5]
#     #x.shadow<-x.store+c(-0.1, 0.1)
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
#     points(x=x.store[2], y=x.store[1], col = "red",pch=20, cex=1)
#     points(x=x.store[2], y=x.store[1], col = "black",cex=0.95)
#     
#   } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
#     multiplier=multiplier+1
#     #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
#     #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)
#     x.plot<-x.current+(multiplier*c(0.15,0))
#     x.ring<-x.plot
#     points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
#     points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     
#     
#   } else {
#     #x.shadow<-x.current+c(-0.1, 0.1)
#     #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
#     points(x=x.plot[2], y=x.plot[1], col="red", pch=20)
#     multiplier=0
#     points(x=x.current[2], y=x.current[1], col="red", pch=20)
#     points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
#   }
#   
#   x.store<-x.current
# }
# 
# 
# legend(x="topright", legend=c("NEO","MNEO","HG"), col=c("blue","red","green"), pch=20, cex=1.5)
# legend(x="topright", legend=c("NEO","MNEO","HG"), col=c("black"), pch=c(1), bg = 'transparent', pt.cex = 0.9, cex=1.5)
# 
# dev.off()
# 
# 


### Newshadow + Proper stacking between pops

#HG minmax
tiff(filename="v50_rebuttal_redo_comparison_longlat_newshad_stack.tiff", res=250, compression = "lzw", height=2000, width=2000)
#par(mfrow=c(1,4))
#plot.new()
map(xlim=c(-12,40), ylim=c(35,70), resolution = 0, bg = "white", fill = T, col="grey90")


grey2<-add.alpha("grey50",0.3)
df<-do.call("rbind", df_split)
df<-df[order(df$LAT,df$POP,decreasing = c(T,T)),]
cols<-as.vector(NULL)
cols[which(df$POP=="HG")]<-"blue"
cols[which(df$POP=="NEO")]<-"red"
cols[which(df$POP=="MNEO")]<-"green"
df$cols<-cols
cols<-as.vector(NULL)

### SHADOW

#points(x=df_split$NEO_NEW$LONG, y=df_split$NEO_NEW$LAT, col=add.alpha("blue", 0.3), pch=20, cex=4)
rm(x.store)
multiplier=0
for(x in 1:nrow(df)){
  x.current<-as.numeric(df[x,4:5])
  x.col<-df[x,]$cols
  
  
  if(!exists("x.store")){
    x.store<-as.numeric(df[x,4:5])
    x.shadow<-as.numeric((x.store))+c(-0.1, 0.1)
    points(x=x.shadow[2],y=x.shadow[1], col="red",pch=20)
    #points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
    #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    if(multiplier > 2){
      x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
      points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20 ,lwd=0)}
    #x.plot<-x.current+(multiplier*c(0.15,0))
    #x.ring<-x.plot
    #points(x=x.ring[2], y=x.ring[1], col = "black")
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    
    
  } else {
    x.shadow<-x.current+c(-0.1, 0.1)
    points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    #points(x=x.current[2], y=x.current[1], col="blue", pch=20)
    multiplier=0
    #points(x=x.current[2], y=x.current[1], col="black")
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
    #x.shadow<-x.store+c(-0.1, 0.1)
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    points(x=x.store[2], y=x.store[1], col = df[x,]$cols,pch=20, cex=1)
    points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    #if(multiplier > 2){
    #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)}
    x.plot<-x.current+(multiplier*c(0.15,0))
    x.ring<-x.plot
    points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
    points(x=x.plot[2], y=x.plot[1], col = df[x,]$cols, pch=20)
    
    
  } else {
    #x.shadow<-x.current+c(-0.1, 0.1)
    #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    #points(x=x.plot[2], y=x.plot[1],col = df[x,]$cols, pch=20)
    multiplier=0
    points(x=x.current[2], y=x.current[1], col = df[x,]$cols, pch=20)
    points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
  }
  
  x.store<-x.current
  
}


df_blu<-df[df$cols=="blue",]

legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c("red","green","blue"), pch=20, cex=1.5)
legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c("black"), pch=c(1), bg = 'transparent', pt.cex = 1.05, cex=1.5)



df_NEO<-df[df$POP=="NEO_NEW",]
df_<-df[df$col=="blue",]



dev.off()



tiff(filename="v50_rebut2_redo_comparison_longlat_newshad_stack_hues.tiff", res=250, compression = "lzw", height=2000, width=2000)
#par(mfrow=c(1,4))
#plot.new()
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

#points(x=df_split$NEO_NEW$LONG, y=df_split$NEO_NEW$LAT, col=add.alpha("blue", 0.3), pch=20, cex=4)
rm(x.store)
multiplier=0
for(x in 1:nrow(df)){
  x.current<-as.numeric(df[x,4:5])
  x.col<-df[x,]$cols
  
  
  if(!exists("x.store")){
    x.store<-as.numeric(df[x,4:5])
    x.shadow<-x.store+c(-0.1, 0.1)
    points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    #points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
    #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    if(multiplier > 2){
      x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
      points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20 ,lwd=0)}
    #x.plot<-x.current+(multiplier*c(0.15,0))
    #x.ring<-x.plot
    #points(x=x.ring[2], y=x.ring[1], col = "black")
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    
    
  } else {
    x.shadow<-x.current+c(-0.1, 0.1)
    points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    #points(x=x.current[2], y=x.current[1], col="blue", pch=20)
    multiplier=0
    #points(x=x.current[2], y=x.current[1], col="black")
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
    #x.shadow<-x.store+c(-0.1, 0.1)
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    points(x=x.store[2], y=x.store[1], col = df[x,]$cols,pch=20, cex=1)
    points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    #if(multiplier > 2){
    #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)}
    x.plot<-x.current+(multiplier*c(0.15,0))
    x.ring<-x.plot
    points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
    points(x=x.plot[2], y=x.plot[1], col = df[x,]$cols, pch=20)
    
    
  } else {
    #x.shadow<-x.current+c(-0.1, 0.1)
    #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    #points(x=x.plot[2], y=x.plot[1],col = df[x,]$cols, pch=20)
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


## Projection test . . .
map(xlim=c(-12,40), ylim=c(35,70), bg = "white", fill = T, col="grey90", lforce = "e")

map(xlim=c(-12,40), ylim=c(35,70), bg = "white", fill = T, col="grey90", orientation = c(120,-180,180), lforce = "e", projection = "globular")

map(xlim=c(-12,40), ylim=c(35,70), bg = "white", fill = T, col="grey90", orientation = c(-50,0,0), lforce = "e", projection = "azequalarea", type="l")
library(mapproj)
zz<-mapproject(x = df$LONG, y=df$LAT, projection = "azequalarea", orientation = c(-50,0,0))



#### Rotation for 3d effect..


### Newshadow + Proper stacking between pops

#HG minmax
#tiff(filename="v50_v2_comparison_longlat_newshad_stack.tiff", res=250, compression = "lzw", height=2000, width=2000)
#par(mfrow=c(1,4))
cols<-as.vector(NULL)
cols[which(df$POP=="HG_NEW")]<-"green"
cols[which(df$POP=="NEO_NEW")]<-"blue"
cols[which(df$POP=="MNEO_NEW")]<-"red"
df$cols<-cols

#plot.new()
par(mar=c(0,0,0,0))
map(xlim=c(-12,40), ylim=c(24,70), bg = "white", fill = T, col="grey90", orientation = c(-80,0,-10), lforce = "e", type="l")

grey2<-add.alpha(grey2,0.3)
df<-do.call("rbind", df_split)
df<-df[order(df$LAT,df$POP,decreasing = c(T,T)),]
df$uid<-seq_along(df[,1])


multi=0
multicount=as.vector(NULL)
for(x in 1:nrow(df)){
  #print(paste("x is", x))
  prev<-df[x-1,]
  
  if(length(unlist(prev))==0){
    # print("First Line - no comparator")
    multicount[x]<-multi
  }
  
  else if(df[x,]$LAT==prev$LAT & df[x,]$LONG==prev$LONG){
    # print("Mulitplier Combo +1")
    multi<-multi+1
    multicount[x]<-multi
    
    
  } else {
    # print("Combo streak broken")
    multi<-0
    multicount[x]<-multi
  }
} ##


df$multicount<-multicount


df$shadowx<-df$LONG+(df$multicount*0.1)
df$shadowy<-df$LAT+(df$multicount*-0.1)



shadowproj<-mapproject(x = df$shadowx, y=df$shadowy, orientation = c(-80,0,-10))

df$shadowx_proj<-shadowproj$x
df$shadowy_proj<-shadowproj$y


df$stacky<-df$LAT+(df$multicount*0.15)

zz<-mapproject(x = df$LONG, y=df$stacky, orientation = c(-80,0,-10))
df$projx<-zz$x
df$projy<-zz$y



points(x=df$shadowx_proj, y=df$shadowy_proj, col=grey2, pch=20, lwd=0)
points(x=df$projx, y=df$projy, col=df$cols, pch=20)
points(x=df$projx, y=df$projy, col="black", pch=1, cex=0.95,lwd=0.5)


points(x=df$projx, y=df$stacky)


### SHADOW

#points(x=df_split$NEO_NEW$LONG, y=df_split$NEO_NEW$LAT, col=add.alpha("blue", 0.3), pch=20, cex=4)
rm(x.store)
multiplier=0
for(x in 1:nrow(df)){
  x.current<-df[x,7:8]
  
  
  
  if(!exists("x.store")){
    x.store<-df[x,7:8]
    x.shadow<-df[x,9:10]
    points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    #points(x=x.store[2], y=x.store[1], col = "blue",pch=20, cex=1)
    #points(x=x.store[2], y=x.store[1], col = "black", cex=1)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    if(multiplier > 2){
      x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
      points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20 ,lwd=0)}
    #x.plot<-x.current+(multiplier*c(0.15,0))
    #x.ring<-x.plot
    #points(x=x.ring[2], y=x.ring[1], col = "black")
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    
    
  } else {
    x.shadow<-x.current+c(-0.1, 0.1)
    points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    #points(x=x.plot[2], y=x.plot[1], col="blue", pch=20)
    #points(x=x.current[2], y=x.current[1], col="blue", pch=20)
    multiplier=0
    #points(x=x.current[2], y=x.current[1], col="black")
  }
  
  x.store<-x.current
  
}


##POLE

rm(x.store)
multiplier=0


for(x in 1:nrow(df)){
  x.current<-df[x,4:5]
  
  
  
  if(!exists("x.store")){
    x.store<-NEO_NEW[x,4:5]
    #x.shadow<-x.store+c(-0.1, 0.1)
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2,pch=20)
    points(x=x.store[2], y=x.store[1], col = cols[x],pch=20, cex=1)
    points(x=x.store[2], y=x.store[1], col = "black", cex=0.95)
    
  } else if (isTRUE(x.current[1]==x.store[1]) & isTRUE(x.current[2]==x.store[2])) {
    multiplier=multiplier+1
    #if(multiplier > 2){
    #x.shadow<-x.store+(multiplier*c(-0.1, 0.1))
    #points(x=x.shadow[2],y=x.shadow[1], col=grey2, pch=20)}
    x.plot<-x.current+(multiplier*c(0.15,0))
    x.ring<-x.plot
    points(x=x.ring[2], y=x.ring[1], col = "black", cex=0.95)
    points(x=x.plot[2], y=x.plot[1], col=cols[x], pch=20)
    
    
  } else {
    #x.shadow<-x.current+c(-0.1, 0.1)
    #points(x=x.shadow[2], y=x.shadow[1], col=grey2, pch=20)
    points(x=x.plot[2], y=x.plot[1], col=cols[x], pch=20)
    multiplier=0
    points(x=x.current[2], y=x.current[1], col=cols[x], pch=20)
    points(x=x.current[2], y=x.current[1], col="black", cex=0.95)
  }
  
  x.store<-x.current
  
}

legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c("blue","red","green"), pch=20, cex=1.5)
legend(x="topright", legend=c("Neolithic","Admixed Neolithic","Mesolithic"), col=c("black"), pch=c(1), bg = 'transparent', pt.cex = 0.9, cex=1.5)

dev.off()





#HG LONGLAT
rect(xright = max(as.numeric(df_split$MNEO_NEW$LONG)), xleft=min(as.numeric(df_split$MNEO_NEW$LONG)), ytop=max(as.numeric(df_split$MNEO_NEW$LAT)), ybot=min(as.numeric(df_split$MNEO_NEW$LAT)), lwd=2, col = new, border=NA)
rect(xright = max((as.numeric(df_split$HG_OLD$LONG)), na.rm = T), xleft=min((as.numeric(df_split$HG_OLD$LONG)), na.rm = T), ytop=max((as.numeric(df_split$HG_OLD$LAT)), na.rm = T), ybot=min(as.numeric(df_split$HG_OLD$LAT), na.rm = T), lwd=2, col = old, border=NA)



MNEO_NEW_MEAN_LAT<-mean(as.numeric(df_split$MNEO_NEW$LAT))
MNEO_NEW_MEAN_LONG<-mean(as.numeric(df_split$MNEO_NEW$LONG))

MNEO_NEW_SD_LAT<-as.numeric(meta_sd$MNEO_NEW[4,2])
MNEO_NEW_SD_LONG<-as.numeric(meta_sd$MNEO_NEW[5,2])

rect(xright=MNEO_NEW_MEAN_LONG, xleft=MNEO_NEW_MEAN_LONG, ytop=(MNEO_NEW_MEAN_LAT + MNEO_NEW_SD_LAT), ybot=(MNEO_NEW_MEAN_LAT - MNEO_NEW_SD_LAT), col = "red", border="deepskyblue", lwd=3)
rect(ytop=MNEO_NEW_MEAN_LAT, ybot=MNEO_NEW_MEAN_LAT, xright=(MNEO_NEW_MEAN_LONG + MNEO_NEW_SD_LONG), xleft=(MNEO_NEW_MEAN_LONG - MNEO_NEW_SD_LONG), col = "red", border="deepskyblue", lwd=3)



HG_OLD_MEAN_LAT<-mean(as.numeric(df_split$HG_OLD$LAT), na.rm = T)
HG_OLD_MEAN_LONG<-mean(as.numeric(df_split$HG_OLD$LONG), na.rm = T)

HG_OLD_SD_LAT<-as.numeric(meta_sd$HG_OLD[4,2])
HG_OLD_SD_LONG<-as.numeric(meta_sd$HG_OLD[5,2])


rect(xright=HG_OLD_MEAN_LONG, xleft=HG_OLD_MEAN_LONG, ytop=(HG_OLD_MEAN_LAT + HG_OLD_SD_LAT), ybot=(HG_OLD_MEAN_LAT - HG_OLD_SD_LAT), col = "red", border="darkolivegreen1", lwd=3)
rect(ytop=HG_OLD_MEAN_LAT, ybot=HG_OLD_MEAN_LAT, xright=(HG_OLD_MEAN_LONG + HG_OLD_SD_LONG), xleft=(HG_OLD_MEAN_LONG - HG_OLD_SD_LONG), col = "red", border="darkolivegreen1", lwd=3)


title(main = "HG")


#map(fill = T, xlim=c(-10,50), ylim=c(20,70), col = "grey50")

# 
# #NEO LONGLAT
# rect(xright = max(as.numeric(df_split$NEO_NEW$LONG)), xleft=min(as.numeric(df_split$NEO_NEW$LONG)), ytop=max(as.numeric(df_split$NEO_NEW$LAT)), ybot=min(as.numeric(df_split$NEO_NEW$LAT)), lwd=2, col = new, border=NA)
# rect(xright = max((as.numeric(df_split$NEO_OLD$LONG)), na.rm = T), xleft=min((as.numeric(df_split$NEO_OLD$LONG)), na.rm = T), ytop=max((as.numeric(df_split$NEO_OLD$LAT)), na.rm = T), ybot=min(as.numeric(df_split$NEO_OLD$LAT), na.rm = T), lwd=2, col = old, border=NA)
# 
# 
# 
# NEO_NEW_MEAN_LAT<-mean(as.numeric(df_split$NEO_NEW$LAT))
# NEO_NEW_MEAN_LONG<-mean(as.numeric(df_split$NEO_NEW$LONG))
# 
# NEO_NEW_SD_LAT<-as.numeric(meta_sd$NEO_NEW[4,2])
# NEO_NEW_SD_LONG<-as.numeric(meta_sd$NEO_NEW[5,2])
# 
# rect(xright=NEO_NEW_MEAN_LONG, xleft=NEO_NEW_MEAN_LONG, ytop=(NEO_NEW_MEAN_LAT + NEO_NEW_SD_LAT), ybot=(NEO_NEW_MEAN_LAT - NEO_NEW_SD_LAT), col = "red", border="deepskyblue", lwd=3)
# rect(ytop=NEO_NEW_MEAN_LAT, ybot=NEO_NEW_MEAN_LAT, xright=(NEO_NEW_MEAN_LONG + NEO_NEW_SD_LONG), xleft=(NEO_NEW_MEAN_LONG - NEO_NEW_SD_LONG), col = "red", border="deepskyblue", lwd=3)
# 
# 
# 
# NEO_OLD_MEAN_LAT<-mean(as.numeric(df_split$NEO_OLD$LAT), na.rm = T)
# NEO_OLD_MEAN_LONG<-mean(as.numeric(df_split$NEO_OLD$LONG), na.rm = T)
# 
# NEO_OLD_SD_LAT<-as.numeric(meta_sd$NEO_OLD[4,2])
# NEO_OLD_SD_LONG<-as.numeric(meta_sd$NEO_OLD[5,2])
# 
# 
# rect(xright=NEO_OLD_MEAN_LONG, xleft=NEO_OLD_MEAN_LONG, ytop=(NEO_OLD_MEAN_LAT + NEO_OLD_SD_LAT), ybot=(NEO_OLD_MEAN_LAT - NEO_OLD_SD_LAT), col = "red", border="darkolivegreen1", lwd=3)
# rect(ytop=NEO_OLD_MEAN_LAT, ybot=NEO_OLD_MEAN_LAT, xright=(NEO_OLD_MEAN_LONG + NEO_OLD_SD_LONG), xleft=(NEO_OLD_MEAN_LONG - NEO_OLD_SD_LONG), col = "red", border="darkolivegreen1", lwd=3)
# 

# 
# map(fill = T, xlim=c(-10,50), ylim=c(20,70), col = "grey50")
# title(main = "MNEO")
# 
# 
# #MNEO LONGLAT
# rect(xright = max(as.numeric(df_split$MNEO_NEW$LONG), na.rm=T), xleft=min(as.numeric(df_split$MNEO_NEW$LONG), na.rm=T), ytop=max(as.numeric(df_split$MNEO_NEW$LAT), na.rm=T), ybot=min(as.numeric(df_split$MNEO_NEW$LAT), na.rm=T), lwd=2, col = new, border=NA)
# rect(xright = max((as.numeric(df_split$MNEO_OLD$LONG)), na.rm = T), xleft=min((as.numeric(df_split$MNEO_OLD$LONG)), na.rm = T), ytop=max((as.numeric(df_split$MNEO_OLD$LAT)), na.rm = T), ybot=min(as.numeric(df_split$MNEO_OLD$LAT), na.rm = T), lwd=2, col = old, border=NA)
# 
# 
# 
# MNEO_NEW_MEAN_LAT<-mean(as.numeric(df_split$MNEO_NEW$LAT), na.rm=T)
# MNEO_NEW_MEAN_LONG<-mean(as.numeric(df_split$MNEO_NEW$LONG), na.rm=T)
# 
# MNEO_NEW_SD_LAT<-as.numeric(meta_sd$MNEO_NEW[4,2])
# MNEO_NEW_SD_LONG<-as.numeric(meta_sd$MNEO_NEW[5,2])
# 
# rect(xright=MNEO_NEW_MEAN_LONG, xleft=MNEO_NEW_MEAN_LONG, ytop=(MNEO_NEW_MEAN_LAT + MNEO_NEW_SD_LAT), ybot=(MNEO_NEW_MEAN_LAT - MNEO_NEW_SD_LAT), col = "red", border="deepskyblue", lwd=3)
# rect(ytop=MNEO_NEW_MEAN_LAT, ybot=MNEO_NEW_MEAN_LAT, xright=(MNEO_NEW_MEAN_LONG + MNEO_NEW_SD_LONG), xleft=(MNEO_NEW_MEAN_LONG - MNEO_NEW_SD_LONG), col = "red", border="deepskyblue", lwd=3)
# 
# 
# 
# MNEO_OLD_MEAN_LAT<-mean(as.numeric(df_split$MNEO_OLD$LAT), na.rm = T)
# MNEO_OLD_MEAN_LONG<-mean(as.numeric(df_split$MNEO_OLD$LONG), na.rm = T)
# 
# MNEO_OLD_SD_LAT<-as.numeric(meta_sd$MNEO_OLD[4,2])
# MNEO_OLD_SD_LONG<-as.numeric(meta_sd$MNEO_OLD[5,2])
# 
# 
# rect(xright=MNEO_OLD_MEAN_LONG, xleft=MNEO_OLD_MEAN_LONG, ytop=(MNEO_OLD_MEAN_LAT + MNEO_OLD_SD_LAT), ybot=(MNEO_OLD_MEAN_LAT - MNEO_OLD_SD_LAT), col = "red", border="darkolivegreen1", lwd=3)
# rect(ytop=MNEO_OLD_MEAN_LAT, ybot=MNEO_OLD_MEAN_LAT, xright=(MNEO_OLD_MEAN_LONG + MNEO_OLD_SD_LONG), xleft=(MNEO_OLD_MEAN_LONG - MNEO_OLD_SD_LONG), col = "red", border="darkolivegreen1", lwd=3)
# 
# #title(main = "NEO")
# 
# 


### Date Ranges

tiff(filename="v50_rebut_fix_metadata_time.tiff",width = 1000, height=1000, res=150, compression = "lzw")
plot.new()
box()
plot.window(xlim = c(0, 4), ylim = c(min(df$AGE)-1000,max(df$AGE)+1000))
grid(nx=F,ny=NULL)
par(mar=c(5, 4, 4, 2) + 0.1)
boxplot(df_split$HG$AGE, df_split$NEO$AGE, df_split$MNEO$AGE, col=c(add.alpha(HG.col, 1),add.alpha(NEO.col,1),add.alpha(MNEO.col,1)), yaxt=F,add=T, lwd=1.5, pch=21)
axis(1, labels = c("Mesolithic","Neolithic","Admixed Neolithic"), at =c(1,2,3))
axis(2)
title(ylab = "Time (kya)")
dev.off()

tiff(filename="v50_rebut_redo_metadata_coveragelog.tiff",width = 1000, height=1000, res=150, compression = "lzw")

par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
plot.window(xlim = c(0, 4), ylim = c(min(df$COV),  1240000))
grid(nx=F,ny=NULL)
boxplot(df_split$HG$COV, df_split$NEO$COV, df_split$MNEO$COV, col=c(add.alpha(HG.col, 1),add.alpha(NEO.col,1),add.alpha(MNEO.col,1)), add=T, lwd=1.5, pch=21)
axis(1, labels = c("Mesolithic","Neolithic","Admixed Neolithic"), at =c(1,2,3))
axis(2)
title(ylab = "Coverage (Number of SNPs)")
dev.off()



tiff(filename="v50_metadata_coverage_discrete.tiff",width = 1000, height=1000, res=150, compression = "lzw")


HGcovs<-read.table("~/slurmp/HG/HG.counts.txt")
NEOcovs<-read.table("~/slurmp/NEO/NEO.counts.txt")
MNEOcovs<-read.table("~/slurmp/MNEO/MNEO.counts.txt")

par(mar=c(5, 4, 4, 2) + 0.1)
plot.new()
plot.window(xlim = c(0, 4), ylim = c(min(-log10(df$COV)-1),  max(-log10(df$COV)))+1)
grid(nx=F,ny=NULL)

list=as.list(NULL)
list[[1]]<-HGcovs
list[[2]]<-NEOcovs
list[[3]]<-MNEOcovs

plot.new()
plot.window(xlim = c(0,4) ,ylim=c(0,1250000))
grid(nx=F,ny=NULL)
boxplot(c(HGcovs, NEOcovs, MNEOcovs), col=c(add.alpha("green", 0.5),add.alpha("blue",0.5),add.alpha("red",0.5)), lwd=1.5, pch=21, add = T,y=1, yaxt=F)
#boxplot(x=NEOcovs, col=c(add.alpha("blue", 0.5)), lwd=1.5, pch=21, add = T, y=2)
#boxplot(x=HGcovs, col=c(add.alpha("red", 0.5)), lwd=1.5, pch=21, add = T,y=3)

#boxplot(list$A, list$B, list$C, col=c(add.alpha("red", 0.5)), lwd=1.5, pch=2, yaxt=F, add=T)




axis(1, labels = c("Mesolithic","Neolithic","Admixed Neolithic"), at =c(1,2,3))
axis(2, cex.axis=0.8)
title(ylab = "Coverage (No. 1240k SNPs)")
dev.off()





dev.off()

#abline(v=MNEO_NEW_MEAN_LONG)
#abline(h=MNEO_NEW_MEAN_LAT)

rect(xright=)


abline(v=MNEO_NEW_MEAN_LONG+as.numeric(meta_sd$HG_OLD[4,2]))

lines(x=c(MNEO_NEW_MEAN_LONG+as.numeric(meta_sd$HG_OLD[4,2])),(MNEO_NEW_MEAN_LONG-as.numeric(meta_sd$HG_OLD[4,2])), y=c(MNEO_NEW_MEAN_LAT+as.numeric(meta_sd$HG_OLD[5,2])))
# NEO minmax

map(fill = T, xlim=c(-10,50), ylim=c(20,70), col = "grey50")


#HG LONGLAT
rect(xright = max(as.numeric(df_split$NEO_NEW$LONG)), xleft=min(as.numeric(df_split$NEO_NEW$LONG)), ytop=max(as.numeric(df_split$NEO_NEW$LAT)), ybot=min(as.numeric(df_split$NEO_NEW$LAT)), lwd=2, col = new)
rect(xright = max((as.numeric(df_split$NEO_OLD$LONG)), na.rm = T), xleft=min((as.numeric(df_split$NEO_OLD$LONG)), na.rm = T), ytop=max((as.numeric(df_split$NEO_OLD$LAT)), na.rm = T), ybot=min(as.numeric(df_split$NEO_OLD$LAT), na.rm = T), lwd=2, col = old)




### MNEO LONGLAT

map(fill = T, xlim=c(-10,50), ylim=c(20,70), col = "grey50")


#MNEO LONGLAT
rect(xright = max(as.numeric(df_split$MNEO_NEW$LONG), na.rm=T), xleft=min(as.numeric(df_split$MNEO_NEW$LONG), na.rm=T), ytop=max(as.numeric(df_split$MNEO_NEW$LAT), na.rm = T), ybot=min(as.numeric(df_split$MNEO_NEW$LAT), na.rm=T), lwd=2, col = new)
rect(xright = max((as.numeric(df_split$MNEO_OLD$LONG)), na.rm = T), xleft=min((as.numeric(df_split$MNEO_OLD$LONG)), na.rm = T), ytop=max((as.numeric(df_split$MNEO_OLD$LAT)), na.rm = T), ybot=min(as.numeric(df_split$MNEO_OLD$LAT), na.rm = T), lwd=2, col = old)





##### old

mapGriddedData(mapRegion = "europe", dataset = "", colourPalette = "white2Black")
title(main="HG")

rect(ybot=48.762753536875-5.74011572082874, ytop = 48.762753536875+5.74011572082874, xright = 17.622089693625+9.96352308496157, xleft = 17.622089693625 - 9.96352308496157, col=old)
rect(ybot=48.8569132246195-5.6688814499509, ytop = 48.8569132246195+5.6688814499509, xright = 13.9614563895978+11.3028992112263, xleft = 13.9614563895978 - 11.3028992112263, col=new)



mapGriddedData(mapRegion = "europe", dataset = "", colourPalette = "white2Black")

title(main="MNEO")

rect(ybot=46.2726637522727-6.81168384963594, ytop = 46.2726637522727+6.81168384963594, xright = 4.5196194122201+12.2356810680617, xleft = 4.5196194122201 - 12.2356810680617, col=old)
rect(ybot=47.0663306146882-5.23095067279302, ytop = 47.0663306146882+5.23095067279302, xright = 7.32524399455234+10.8575950664404, xleft = 7.32524399455234 - 10.8575950664404, col=new)



mapGriddedData(mapRegion = "europe", dataset = "", colourPalette = "white2Black")

title(main="NEO")

rect(ybot=43.5492353927907-3.17288872325644, ytop = 43.5492353927907+3.17288872325644, xright = 23.376492025814+5.84653732934908, xleft = 23.376492025814 - 5.84653732934908, col=old)
rect(ybot=46.418000382186-4.4567753330452, ytop = 46.418000382186+4.4567753330452, xright = 17.5723644865582+7.6547248293931, xleft = 17.5723644865582 - 7.6547248293931, col=new)

NEO_old=c(7809.82978723404-867.284078693523,7809.82978723404,7809.82978723404+867.284078693523)
NEO_new=c(7324.10697674419-540.690344856548,7324.10697674419,7324.10697674419+540.690344856548)


HG_old=c(8641.0-2088.24389709883,8641.0,8641.0+2088.24389709883)
MNEO_NEW=c(9045.22282608696-2874.53849969314,9045.22282608696,9045.22282608696+2874.53849969314)

MNEO_old=c(5175.40486725664-1111.08931930022,5175.40486725664,5175.40486725664+1111.08931930022)
MNEO_new=c(5867.59513274336-943.853298172911,5867.59513274336,5867.59513274336+943.853298172911)



boxplot.df<-cbind(NEO_old,NEO_new,HG_old,MNEO_NEW,MNEO_old,MNEO_new)
plot.window(xlim=c(0,5), ylim=c(0,1))
boxplot(x = boxplot.df, axes=F, col=c(old, new))
axis(2)
axis(1, at = c(1.5,3.5,5.5), labels=c("NEO","HG","MNEO"))
title(main=expression("Differences in Ages of Populations"), xlab = "Population",ylab="Years (bp)")
legend(x="topright", pch=15, col=c(old,new), legend = c("Old","New"))




NEO_old=c(1.5988465106383-1.81267019994794,1.5988465106383,1.5988465106383+1.81267019994794)
NEO_new=c(2.04415564651163-4.13083336456383,2.04415564651163,2.04415564651163+4.13083336456383)


HG_old=c(2.35931909016393-2.34963405823734,2.35931909016393,2.35931909016393+2.34963405823734)
MNEO_NEW=c(2.51476001086957-5.74221698027765,2.51476001086957,2.51476001086957+5.74221698027765)

MNEO_old=c(1.60437907743363-3.43924962311748,1.60437907743363,1.60437907743363+3.43924962311748)
MNEO_new=c(1.38306888314607-1.38306888314607,1.38306888314607,1.38306888314607+1.38306888314607)

boxplot.cov.df<-cbind(NEO_old,NEO_new,HG_old,MNEO_NEW,MNEO_old,MNEO_new)
plot.window(xlim=c(0,5), ylim=c(0,1))
boxplot(x = boxplot.cov.df, axes=F, col=c(old, new))
axis(2)
axis(1, at = c(1.5,3.5,5.5), labels=c("NEO","HG","MNEO"))
title(main=expression("Differences in Coverage of Populations"), xlab = "Population",ylab="Coverage")
legend(x="topright", pch=15, col=c(old,new), legend = c("Old","New"))



# numbers

MNEO_NEW=184
HG_old=122


NEO_old=47
NEO_new=215

MNEO_old=
  MNEO_new=904

# 
# metadata.hg_old.txt.txt
# 8641.0 2088.24389709883 48.762753536875 5.74011572082874 17.622089693625 9.96352308496157 2.35931909016393 2.34963405823734
# metadata.mneo_old.txt.txt
# 5175.40486725664 1111.08931930022 46.2726637522727 6.81168384963594 4.5196194122201 12.2356810680617 1.60437907743363 3.43924962311748
# metadata.neo_old.txt.txt
# 7809.82978723404 867.284078693523 43.5492353927907 3.17288872325644 23.376492025814 5.84653732934908 1.5988465106383 1.81267019994794
# metadata.MNEO_NEW.txt.txt
# 9045.22282608696 2874.53849969314 48.8569132246195 5.6688814499509 13.9614563895978 11.3028992112263 2.51476001086957 5.74221698027765
# metadata.mneo_new.txt.txt
# 5867.59513274336 943.853298172911 47.0663306146882 5.23095067279302 7.32524399455234 10.8575950664404 1.38306888314607 1.38306888314607
# metadata.neo_new.txt.txt
# 7324.10697674419 540.690344856548 46.418000382186 4.4567753330452 17.5723644865582 7.6547248293931 2.04415564651163 4.13083336456383

