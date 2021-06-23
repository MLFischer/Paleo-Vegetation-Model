# MLFischer 2020-08
# Fischer et al., 2021
# A Phytolith supported Biosphere-Hydrosphere Predictive Model for southern Ethiopia: 
# insights into paleoenvironmental changes and human landscape preferences since the Last Glacial Maximum
###
# LC catchment for TS 1 to 3 
###

remove(list=ls()) #clean up#
removeTmpFiles(0)

#get some packages
library(raster)
library(rgdal)  
library(gbm)
library(dismo)
library(rasterVis)
library(dplyr)

#load shapefiles and data
setwd("../01_data/00_Shapes/")
#area.shp <- readOGR(".",layer = "20190902_masc_studyarea_stage1")
#a.shp <- readOGR(".",layer = "20191024_abaya")
#b.shp <- readOGR(".",layer = "20191024_bahir")
#c.shp <- readOGR(".",layer = "20191024_chamo")

get(load("20201119_masc_studyarea_stage1.RData"))
get(load("20201119__abaya.RData"))
get(load("20201119__bahir.RData"))
get(load("20201119__chamo.RData"))

modistsp_file         <- "../01_data/02_MODIS/LandCover_Type_Yearly_500m_v6/LC2/MCD12Q1_LC2_2018_001.tif"
MCD12Q1_2001          <-raster(modistsp_file)

path_2ndlevel<- "../02_data_2nd_level/"
path_3ndlevel<- "../02_data_3nd_brt/"
path_4ndlevel<- "../02_data_4nd_scenarios/"
path_img     <- "../00_Results/20200705_GIF"

setwd(path_3ndlevel)
get(load("20200821_gbm_LC_final.RData"))

setwd(path_2ndlevel)
get(load("data_all_pr.RData"))
get(load("data_all.RData"))
get(load("mask_all.RData"))
#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################
#TIME SLICE 1: MODERN DAY 
d_pr <- as.data.frame(data_all_pr,na.rm=F,xy=T)

par(mfrow=c(1,1))

itree <-    gbm.perf(gbm_LC, plot.it = F, oobag.curve = FALSE, overlay = TRUE)[1]
gbm_pr <-predict.gbm(gbm_LC , d_pr[c(3:15)]     ,type = "response",  n.trees=itree) #prediction 

d_pr[,16]<- apply(gbm_pr, 1, function(row) which.max(row))

d_pr[which(d_pr[,16]==4),16]<-10 
d_pr[which(d_pr[,16]==3),16]<-9 
d_pr[which(d_pr[,16]==2),16]<-7 
d_pr[which(d_pr[,16]==1),16]<-2 
d_pr[which(is.na(d_pr[,4])==T),16]<-NA

r3<- rasterFromXYZ(d_pr[,c(1,2,16)])
r3@crs@projargs<- MCD12Q1_2001@crs@projargs
r3          <- projectRaster(r3,crs= area.shp@proj4string,method="ngb" )
masc_all    <- projectRaster(masc_all,crs= area.shp@proj4string,method="ngb" )
data_all_lc <- projectRaster(data_all$lc,crs= area.shp@proj4string,method="ngb" )
r_pr_MD<- r3
r3<-r3*masc_all
r_dif<- r3-data_all_lc
r_dif[which(r_dif[]!=0)]<-1
############################################################################################################
##probability raster for each class:
# ..$ : chr [1:4] "2 Evergreen Broadleaf Forest" "7 Open Shrubland" "9 Savanna" "10 Grassland"
d_pr[,17]<- gbm_pr[,1,]
d_pr[,18]<- gbm_pr[,2,]
d_pr[,19]<- gbm_pr[,3,]
d_pr[,20]<- gbm_pr[,4,]
par(mfrow=c(2,2))
for(i in 17:20){
  r4              <- rasterFromXYZ(d_pr[,c(1,2,i)])
  r4@crs@projargs <- MCD12Q1_2001@crs@projargs
  r4              <- projectRaster(r4,crs= area.shp@proj4string,method="ngb" )
  plot(r4)
  oldwd<-getwd()
  setwd("../00_Results")
  writeRaster(r4,paste("04_20200831_MD_pr_", colnames(gbm_pr)[i-16], ".tif",sep=""),overwrite=TRUE)
  setwd(oldwd)
  
}
d_pr[which(d_pr[,17]>0.5),17] <- d_pr[which(d_pr[,17]>0.5),17]+2
d_pr[which(d_pr[,18]>0.5),18] <- d_pr[which(d_pr[,18]>0.5),18]+7
d_pr[which(d_pr[,19]>0.5),19] <- d_pr[which(d_pr[,19]>0.5),19]+9
d_pr[which(d_pr[,20]>0.5),20] <- d_pr[which(d_pr[,20]>0.5),20]+10
d_pr[,21]<- d_pr[,17]+d_pr[,18]+d_pr[,19]+d_pr[,20]
r4              <- rasterFromXYZ(d_pr[,c(1,2,21)])
r4@crs@projargs <- MCD12Q1_2001@crs@projargs
r4              <- projectRaster(r4,crs= area.shp@proj4string,method="ngb" )
par(mfrow=c(1,1))
plot(r4)
hist(d_pr[,21])
############################################################################################################


###
oldwd<-getwd()
setwd(path_4ndlevel)
writeRaster(r3,"04_20201117_PI_pr.tif",overwrite=TRUE)
writeRaster(r_pr_MD,"04_20201117_PI_pr_all.tif",overwrite=TRUE)
writeRaster(r_dif,"04_20201117_PI_pr_dif.tif",overwrite=TRUE)
setwd(oldwd)
#################################################################################################
#################################################################################################

###elevation gradient PI 
par(mfrow=c(1,1))


d_pr_l <- d_pr[,c(15,16)]
d_pr_l$srtm[d_pr_l$srtm>2000]<- 3
d_pr_l$srtm[d_pr_l$srtm>1000]<- 2
d_pr_l$srtm[d_pr_l$srtm>5]   <- 1


require(ggplot2)

require(dplyr)
#
df <- na.omit(d_pr_l)
df_d <- df %>% 
  count(srtm, tpi)%>%
  group_by(srtm) %>%
  mutate(freq = n / sum(n))

df_d$tpi<- as.factor(df_d$tpi)
eg_1 <- ggplot(df_d, aes(x="", y=freq,fill=tpi)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void() + facet_wrap(~srtm,ncol = 4)
plot(eg_1) 

####elevation gradient MD
d_pr_l <- as.data.frame(stack(data_all$lc, data_all$srtm))
d_pr_l$srtm[d_pr_l$srtm>2000]<- 3
d_pr_l$srtm[d_pr_l$srtm>1000]<- 2
d_pr_l$srtm[d_pr_l$srtm>5]   <- 1


require(ggplot2)

require(dplyr)
#
df <- na.omit(d_pr_l)
df_d <- df %>% 
  count(srtm, lc)%>%
  group_by(srtm) %>%
  mutate(freq = n / sum(n))

df_d$lc<- as.factor(df_d$lc)
eg_1md <- ggplot(df_d, aes(x="", y=freq,fill=lc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void() + facet_wrap(~srtm,ncol = 4)
plot(eg_1md) 
########################################################
####MD new with unmasked dataset
setwd(path_4ndlevel)

r_md_lc <- raster("01_20201116_LC_MD_all.tif")
plot(r_md_lc)

d_pr_l <- as.data.frame(stack(r_md_lc, data_all_pr$srtm))
colnames(d_pr_l)<- c("lc","srtm")
d_pr_l$srtm[d_pr_l$srtm>2000]<- 3
d_pr_l$srtm[d_pr_l$srtm>1000]<- 2
d_pr_l$srtm[d_pr_l$srtm>5]   <- 1
df <- na.omit(d_pr_l)

df_d <- df %>% 
  count(srtm, lc)%>%
  group_by(srtm) %>%
  mutate(freq = n / sum(n))

df_d$lc<- as.factor(df_d$lc)
eg_1md <- ggplot(df_d, aes(x="", y=freq,fill=lc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void() + facet_wrap(~srtm,ncol = 4)
plot(eg_1md) 


#################################################################################################
#################################################################################################
# Time Slice 2a: AHP 
# NO CAB
#
df_s <- data.frame(matrix(ncol = 15, nrow = 101))
colnames(df_s) <- c("pre", "a_2","a_7","a_9","a_10", "b_2","b_7","b_9","b_10","c_2","c_7","c_9","c_10")
d_pr_n <- as.data.frame(data_all_pr,na.rm=F,xy=T)
  
for (i in 1:101){

  d_pr <- d_pr_n
  j<- i-1
  f_pr <- 0.5+   j/100


d_pr[c(3:14)]<- d_pr[c(3:14)]*f_pr

itree <-    gbm.perf(gbm_LC, plot.it = F, oobag.curve = FALSE, overlay = TRUE)[1]

gbm_pr <-predict.gbm(gbm_LC , d_pr[c(3:15)]     ,type = "response",  n.trees=itree) #prediction 

#unique(apply(gbm_pr, 1, function(row) which.max(row)))
#unique(colnames(gbm_pr)[apply(gbm_pr, 1, function(row) which.max(row))])

d_pr[,16]<- apply(gbm_pr, 1, function(row) which.max(row))


d_pr[which(d_pr[,16]==4),16]<-10 
d_pr[which(d_pr[,16]==3),16]<-9 
d_pr[which(d_pr[,16]==2),16]<-7 
d_pr[which(d_pr[,16]==1),16]<-2 
d_pr[which(is.na(d_pr[,4])==T),16]<-NA


r3<- rasterFromXYZ(d_pr[,c(1,2,16)])
r3@crs@projargs<- MCD12Q1_2001@crs@projargs
r3<-projectRaster(r3,crs= area.shp@proj4string,method="ngb" )

r3_a <- mask(r3,a.shp)
#plot(r3_a)

r3_b <- mask(r3,b.shp)
#plot(r3_b)

r3_c <- mask(r3,c.shp)
#plot(r3_c)

df_s$pre[i]<-f_pr

a_all <-freq(r3_a,value=2)+freq(r3_a,value=7)+freq(r3_a,value=8)+freq(r3_a,value=9)+freq(r3_a,value=10)
b_all <-freq(r3_b,value=2)+freq(r3_b,value=7)+freq(r3_b,value=8)+freq(r3_b,value=9)+freq(r3_b,value=10)
c_all <-freq(r3_c,value=2)+freq(r3_c,value=7)+freq(r3_c,value=8)+freq(r3_c,value=9)+freq(r3_c,value=10)

df_s$a_2[i]<- freq(r3_a,value=2)/a_all
df_s$b_2[i]<- freq(r3_b,value=2)/b_all
df_s$c_2[i]<- freq(r3_c,value=2)/c_all

df_s$a_7[i]<- freq(r3_a,value=7)/a_all
df_s$b_7[i]<- freq(r3_b,value=7)/b_all
df_s$c_7[i]<- freq(r3_c,value=7)/c_all

df_s$a_9[i]<- freq(r3_a,value=9)/a_all
df_s$b_9[i]<- freq(r3_b,value=9)/b_all
df_s$c_9[i]<- freq(r3_c,value=9)/c_all

df_s$a_10[i]<- freq(r3_a,value=10)/a_all
df_s$b_10[i]<- freq(r3_b,value=10)/b_all
df_s$c_10[i]<- freq(r3_c,value=10)/c_all

print(i)
}
olwd<- getwd()
setwd(path_4ndlevel)
save(df_s,file= "20201120_AHP01_noCAB_LC.RData")
#setwd(olwd)


require(ggplot2)
require(reshape2)

d <- melt(df_s, id.vars="pre")

ggplot(d[1:404,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_AHP01_noCAB_LC_Abaya",".png"))


  ggplot(d[405:808,], aes(pre,value, col=variable)) + 
    geom_line() + 
    geom_point()
  ggsave(paste("20201104_AHP01_noCAB_LC_Bahir",".png"))
  
  ggplot(d[809:1212,], aes(pre,value, col=variable)) + 
    geom_line() + 
    geom_point()
  ggsave(paste("20201104_AHP01_noCAB_LC_Chamo",".png"))
  
  
  
  
  
  
  


  
#################################################################################################
#################################################################################################
  # Time Slice 2b: AHP 
  # CAB Szenario
  # Temperature: 
  #
  df_s <- data.frame(matrix(ncol = 13, nrow = 101))
  colnames(df_s) <- c("pre", "a_2","a_7","a_9","a_10", "b_2","b_7","b_9","b_10","c_2","c_7","c_9","c_10")
  
  
  d_pr_n <- as.data.frame(data_all_pr,na.rm=F,xy=T)
  d_pr_n$gpm_a <- d_pr_n[c(3)]+d_pr_n[c(4)]+d_pr_n[c(5)]+d_pr_n[c(6)]+d_pr_n[c(7)]+d_pr_n[c(8)]+d_pr_n[c(9)]+d_pr_n[c(10)]+d_pr_n[c(11)]+d_pr_n[c(12)]+d_pr_n[c(13)]+d_pr_n[c(14)]


  d_pr_n2<- d_pr_n
  
  d_pr_n2<- d_pr_n
  
  for (i in 1:101){
    d_pr_n <- d_pr_n2
    
    j<- i-1
    f_pr <- 0.5+   j/100
    
    
    fac <- 1+  ( d_pr_n[,17]*(f_pr-1)/  (d_pr_n[,8]+d_pr_n[,9]+d_pr_n[,10]+d_pr_n[,11])   )
    
    d_pr_n[,c(8)] <- d_pr_n[,c(8)]* fac[,1]
    d_pr_n[,c(9)] <- d_pr_n[,c(9)]* fac[,1]
    d_pr_n[,c(10)] <- d_pr_n[,c(10)]* fac[,1]
    d_pr_n[,c(11)] <- d_pr_n[,c(11)]* fac[,1]
    

    d_pr <- d_pr_n
    itree <-    gbm.perf(gbm_LC, plot.it = F, oobag.curve = FALSE, overlay = TRUE)[1]
    
    gbm_pr <-predict.gbm(gbm_LC , d_pr[c(3:15)]     ,type = "response",  n.trees=itree) #prediction 
    
    #unique(apply(gbm_pr, 1, function(row) which.max(row)))
    #unique(colnames(gbm_pr)[apply(gbm_pr, 1, function(row) which.max(row))])
    
    d_pr[,16]<- apply(gbm_pr, 1, function(row) which.max(row))
    
    
    d_pr[which(d_pr[,16]==4),16]<-10 
    d_pr[which(d_pr[,16]==3),16]<-9 
    d_pr[which(d_pr[,16]==2),16]<-7 
    d_pr[which(d_pr[,16]==1),16]<-2 
    d_pr[which(is.na(d_pr[,4])==T),16]<-NA
    
    
    r3<- rasterFromXYZ(d_pr[,c(1,2,16)])
    r3@crs@projargs<- MCD12Q1_2001@crs@projargs
    r3<-projectRaster(r3,crs= area.shp@proj4string,method="ngb" )
    plot(r3)
    
    r3_a <- mask(r3,a.shp)
    #plot(r3_a)
    
    r3_b <- mask(r3,b.shp)
    #plot(r3_b)
    
    r3_c <- mask(r3,c.shp)
    #plot(r3_c)
    
    df_s$pre[i]<-f_pr
    
    a_all <-freq(r3_a,value=2)+freq(r3_a,value=7)+freq(r3_a,value=8)+freq(r3_a,value=9)+freq(r3_a,value=10)
    b_all <-freq(r3_b,value=2)+freq(r3_b,value=7)+freq(r3_b,value=8)+freq(r3_b,value=9)+freq(r3_b,value=10)
    c_all <-freq(r3_c,value=2)+freq(r3_c,value=7)+freq(r3_c,value=8)+freq(r3_c,value=9)+freq(r3_c,value=10)
    
    df_s$a_2[i]<- freq(r3_a,value=2)/a_all
    df_s$b_2[i]<- freq(r3_b,value=2)/b_all
    df_s$c_2[i]<- freq(r3_c,value=2)/c_all
    
    df_s$a_7[i]<- freq(r3_a,value=7)/a_all
    df_s$b_7[i]<- freq(r3_b,value=7)/b_all
    df_s$c_7[i]<- freq(r3_c,value=7)/c_all
    
    df_s$a_9[i]<- freq(r3_a,value=9)/a_all
    df_s$b_9[i]<- freq(r3_b,value=9)/b_all
    df_s$c_9[i]<- freq(r3_c,value=9)/c_all
    
    df_s$a_10[i]<- freq(r3_a,value=10)/a_all
    df_s$b_10[i]<- freq(r3_b,value=10)/b_all
    df_s$c_10[i]<- freq(r3_c,value=10)/c_all
    
    print(i)
  }
  olwd<- getwd()
  setwd(path_4ndlevel)
  save(df_s,file= "20201120_AHP01_CAB_LC.RData")
  #setwd(olwd)
  
  
  require(ggplot2)
  require(reshape2)
  
  d <- melt(df_s, id.vars="pre")
  
  ggplot(d[1:404,], aes(pre,value, col=variable)) + 
    geom_line() + 
    geom_point()
  ggsave(paste("20201120_AHP01_CAB_LC_Abaya",".png"))
  
  
  ggplot(d[405:808,], aes(pre,value, col=variable)) + 
    geom_line() + 
    geom_point()
  ggsave(paste("20201120_AHP01_CAB_LC_Bahir",".png"))
  
  ggplot(d[809:1212,], aes(pre,value, col=variable)) + 
    geom_line() + 
    geom_point()
  ggsave(paste("20201120_AHP01_CAB_LC_Chamo",".png"))
  
#################################################################################################
#################################################################################################
# Time Slice 3: LGM
# Temperature change based on Loomis et al. 2012, 2017
  
d_pr <- as.data.frame(data_all_pr,na.rm=F,xy=T)

# Temperature lowering for LGM lower liit based on: 
# Loomis 2017 laps rate differenfce of 0.9 over whole Eastern Africa (6.7 LGM, 5.9 PI)
# Loomis 2012 general cooling of 3 to 4°C

d_pr_tl <- d_pr
d_pr_tl$srtm <- 1000*(0.9*(d_pr_tl$srtm/1000) +3)/5.8 + d_pr_tl$srtm 
  
  
d_pr_tu <- d_pr
d_pr_tu$srtm <- 1000*(0.9*(d_pr_tu$srtm/1000) +4)/5.8 + d_pr_tu$srtm 

###kernel densitiy plot
e_md_d <- density(d_pr$srtm ,na.rm=T)
e_tl_d <- density(d_pr_tl$srtm ,na.rm=T)
e_tu_d <- density(d_pr_tu$srtm ,na.rm=T)
par(mfrow=c(1,1))
plot(e_md_d)
lines(e_tl_d,col="blue")
lines(e_tu_d,col="red")

##########
################################################
df_s <- data.frame(matrix(ncol = 16, nrow = 101))
colnames(df_s) <- c("pre", "a_2","a_7","a_9","a_10","a_16", "b_2","b_7","b_9","b_10","b_16","c_2","c_7","c_9","c_10","c_16")
d_pr_n <- d_pr_tl

for (i in 1:101){
  
  d_pr <- d_pr_n
  j<- i-1
  f_pr <- 0.5+   j/100
  
  
  d_pr[c(3:14)]<- d_pr[c(3:14)]*f_pr
  
  itree <-    gbm.perf(gbm_LC, plot.it = F, oobag.curve = FALSE, overlay = TRUE)[1]
  
  gbm_pr <-predict.gbm(gbm_LC , d_pr[c(3:15)]     ,type = "response",  n.trees=itree) #prediction 
  
  #unique(apply(gbm_pr, 1, function(row) which.max(row)))
  #unique(colnames(gbm_pr)[apply(gbm_pr, 1, function(row) which.max(row))])
  
  d_pr[,16]<- apply(gbm_pr, 1, function(row) which.max(row))
  
  #elevation threshold for trees
  d_pr[which(d_pr[,15]>3000),16]<-16
  
  d_pr[which(d_pr[,16]==4),16]<-10 
  d_pr[which(d_pr[,16]==3),16]<-9 
  d_pr[which(d_pr[,16]==2),16]<-7 
  d_pr[which(d_pr[,16]==1),16]<-2 
  d_pr[which(is.na(d_pr[,4])==T),16]<-NA
  
  
  r3<- rasterFromXYZ(d_pr[,c(1,2,16)])
  r3@crs@projargs<- MCD12Q1_2001@crs@projargs
  r3<-projectRaster(r3,crs= area.shp@proj4string,method="ngb" )
  
  r3_a <- mask(r3,a.shp)
  #plot(r3_a)
  
  r3_b <- mask(r3,b.shp)
  #plot(r3_b)
  
  r3_c <- mask(r3,c.shp)
  #plot(r3_c)
  
  df_s$pre[i]<-f_pr
  
  a_all <-freq(r3_a,value=2)+freq(r3_a,value=7)+freq(r3_a,value=8)+freq(r3_a,value=9)+freq(r3_a,value=10)+freq(r3_a,value=16)
  b_all <-freq(r3_b,value=2)+freq(r3_b,value=7)+freq(r3_b,value=8)+freq(r3_b,value=9)+freq(r3_b,value=10)+freq(r3_b,value=16)
  c_all <-freq(r3_c,value=2)+freq(r3_c,value=7)+freq(r3_c,value=8)+freq(r3_c,value=9)+freq(r3_c,value=10)+freq(r3_c,value=16)
  
  df_s$a_2[i]<- freq(r3_a,value=2)/a_all
  df_s$b_2[i]<- freq(r3_b,value=2)/b_all
  df_s$c_2[i]<- freq(r3_c,value=2)/c_all
  
  df_s$a_7[i]<- freq(r3_a,value=7)/a_all
  df_s$b_7[i]<- freq(r3_b,value=7)/b_all
  df_s$c_7[i]<- freq(r3_c,value=7)/c_all
  
  df_s$a_9[i]<- freq(r3_a,value=9)/a_all
  df_s$b_9[i]<- freq(r3_b,value=9)/b_all
  df_s$c_9[i]<- freq(r3_c,value=9)/c_all
  
  df_s$a_10[i]<- freq(r3_a,value=10)/a_all
  df_s$b_10[i]<- freq(r3_b,value=10)/b_all
  df_s$c_10[i]<- freq(r3_c,value=10)/c_all
  
  df_s$a_16[i]<- freq(r3_a,value=16)/a_all
  df_s$b_16[i]<- freq(r3_b,value=16)/b_all
  df_s$c_16[i]<- freq(r3_c,value=16)/c_all
  
  #print(paste(i,c_all,df_s$c_10[i],freq(r3_c,value=10)))
  print(i)
}
olwd<- getwd()
setwd(path_4ndlevel)
save(df_s,file= "20201104_LGM01_TL_LC.RData")
#setwd(olwd)


require(ggplot2)
require(reshape2)

d <- melt(df_s, id.vars="pre")

ggplot(d[1:505,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Abaya",".png"))


ggplot(d[506:1010,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Bahir",".png"))

ggplot(d[1011:1515,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Chamo",".png"))

#################################################################################################
################################################
################################################
#UPPER TEMPERATURE BOUNDARY
df_s <- data.frame(matrix(ncol = 16, nrow = 101))
colnames(df_s) <- c("pre", "a_2","a_7","a_9","a_10","a_16", "b_2","b_7","b_9","b_10","b_16","c_2","c_7","c_9","c_10","c_16")
d_pr_n <- d_pr_tu

for (i in 1:101){
  
  d_pr <- d_pr_n
  j<- i-1
  f_pr <- 0.5+   j/100
  
  
  d_pr[c(3:14)]<- d_pr[c(3:14)]*f_pr
  
  itree <-    gbm.perf(gbm_LC, plot.it = F, oobag.curve = FALSE, overlay = TRUE)[1]
  gbm_pr <-predict.gbm(gbm_LC , d_pr[c(3:15)]     ,type = "response",  n.trees=itree) #prediction 
  
  #unique(apply(gbm_pr, 1, function(row) which.max(row)))
  #unique(colnames(gbm_pr)[apply(gbm_pr, 1, function(row) which.max(row))])
  
  d_pr[,16]<- apply(gbm_pr, 1, function(row) which.max(row))
  
  d_pr[which(d_pr[,15]>3000),16]<-16
  
  d_pr[which(d_pr[,16]==4),16]<-10 
  d_pr[which(d_pr[,16]==3),16]<-9 
  d_pr[which(d_pr[,16]==2),16]<-7 
  d_pr[which(d_pr[,16]==1),16]<-2 
  d_pr[which(is.na(d_pr[,4])==T),16]<-NA
  
  
  r3<- rasterFromXYZ(d_pr[,c(1,2,16)])
  r3@crs@projargs<- MCD12Q1_2001@crs@projargs
  r3<-projectRaster(r3,crs= area.shp@proj4string,method="ngb" )
  
  r3_a <- mask(r3,a.shp)
  #plot(r3_a)
  
  r3_b <- mask(r3,b.shp)
  #plot(r3_b)
  
  r3_c <- mask(r3,c.shp)
  #plot(r3_c)
  
  df_s$pre[i]<-f_pr
  
  a_all <-freq(r3_a,value=2)+freq(r3_a,value=7)+freq(r3_a,value=8)+freq(r3_a,value=9)+freq(r3_a,value=10)+freq(r3_a,value=16)
  b_all <-freq(r3_b,value=2)+freq(r3_b,value=7)+freq(r3_b,value=8)+freq(r3_b,value=9)+freq(r3_b,value=10)+freq(r3_b,value=16)
  c_all <-freq(r3_c,value=2)+freq(r3_c,value=7)+freq(r3_c,value=8)+freq(r3_c,value=9)+freq(r3_c,value=10)+freq(r3_c,value=16)
  
  df_s$a_2[i]<- freq(r3_a,value=2)/a_all
  df_s$b_2[i]<- freq(r3_b,value=2)/b_all
  df_s$c_2[i]<- freq(r3_c,value=2)/c_all
  
  df_s$a_7[i]<- freq(r3_a,value=7)/a_all
  df_s$b_7[i]<- freq(r3_b,value=7)/b_all
  df_s$c_7[i]<- freq(r3_c,value=7)/c_all
  
  df_s$a_9[i]<- freq(r3_a,value=9)/a_all
  df_s$b_9[i]<- freq(r3_b,value=9)/b_all
  df_s$c_9[i]<- freq(r3_c,value=9)/c_all
  
  df_s$a_10[i]<- freq(r3_a,value=10)/a_all
  df_s$b_10[i]<- freq(r3_b,value=10)/b_all
  df_s$c_10[i]<- freq(r3_c,value=10)/c_all
  
  df_s$a_16[i]<- freq(r3_a,value=16)/a_all
  df_s$b_16[i]<- freq(r3_b,value=16)/b_all
  df_s$c_16[i]<- freq(r3_c,value=16)/c_all
  
  print(i)
}
olwd<- getwd()
setwd(path_4ndlevel)
save(df_s,file= "20201104_LGM01_TU_LC.RData")
#setwd(olwd)


require(ggplot2)
require(reshape2)

d <- melt(df_s, id.vars="pre")

ggplot(d[1:505,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Abaya",".png"))


ggplot(d[506:1010,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Bahir",".png"))

ggplot(d[1011:1515,], aes(pre,value, col=variable)) + 
  geom_line() + 
  geom_point()
ggsave(paste("20201104_LGM01_TL_LC_Chamo",".png"))

 #################################################################################################