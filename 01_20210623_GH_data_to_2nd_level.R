# MLFischer 2020-08
# Fischer et al., 2021
# A Phytolith supported Biosphere-Hydrosphere Predictive Model for southern Ethiopia: 
# insights into paleoenvironmental changes and human landscape preferences since the Last Glacial Maximum
###
# data preprocessing
###
remove(list=ls()) #clean up#
removeTmpFiles(0)
#

library(raster)
library(rgdal)  

setwd("../01_data/00_Shapes/")
area.shp <- readOGR(".",layer = "20190902_masc_studyarea_stage1")
osm.shp  <- readOGR(".",layer = "20190908_masc_OSM_streets_500m")
osm.masc <- raster("masc_OSM_1km.tif")

##################################################################################################################
##################################################################################################################
##MCD12Q1
#load a single file
modistsp_file         <- "../01_data/02_MODIS/VI_Monthly_1Km_v6/EVI/MOD13A3_EVI_2001_001.tif"
MOD13A3_EVI_2001_001  <-raster(modistsp_file) #this is the resample standard

modistsp_file         <- "../01_data/02_MODIS/LandCover_Type_Yearly_500m_v6/LC2/MCD12Q1_LC2_2018_001.tif"
MCD12Q1_2001          <-raster(modistsp_file)

area.shp@proj4string #check projections of MODIS and shapefile of study area and reproject the shapefile
MCD12Q1_2001@crs@projargs
area.shp<- spTransform(area.shp,MCD12Q1_2001@crs@projargs) #transform shapefile 
osm.shp<- spTransform(osm.shp,MCD12Q1_2001@crs@projargs) #transform shapefile 

plot(MCD12Q1_2001);plot(area.shp,add=T) #show both, if it works, good ;-)


#load the raster stat time series
in_virtual_file       <- "../01_data/02_MODIS/LandCover_Type_Yearly_500m_v6/Time_Series/RData/Terra/LC2/MCD12Q1_LC2_1_2001_1_2018_RData.RData"
MCD12Q1_TS_2001_2018  <-get(load(in_virtual_file))
plot(MCD12Q1_TS_2001_2018)


par(mfrow=c(1,2))
#processing loop for all years
#1. masc 
#2. set NA for all unwanted targets (agriculture, water, urban, ...)
MCD_masc <- MCD12Q1_TS_2001_2018[[1]]
for (i in 1:18){
  
  plot(MCD12Q1_TS_2001_2018[[i]]);plot(area.shp,add=T) #show both, if it works, good ;-)
  
  MCD12Q1_TS_2001_2018[[i]]  <- mask(MCD12Q1_TS_2001_2018[[i]],area.shp)
  #Class Value Class Description MODIS
  #0 Water Bodies
  #1 Evergreen Needleleaf Forest
  #2 Evergreen Broadleaf Forest
  #3 Deciduous Needleleaf Forest
  #4 Deciduous Broadleaf Forest
  #5 Mixed Forest
  #6 Closed Shrubland
  #7 Open Shrubland
  #8 Woody Savanna
  #9 Savanna
  #10 Grassland
  #11 Permanent Wetland
  #12 Cropland
  #13 Urban or Built-Up
  #14 Croptland/ Natural Vegetation Mosaics
  #15 Barren or Sparsely Vegetated/ Non vegetated
  #254 Unclassified
  #255 Missing Data
  #d_m<-d[c(which(d$MCD12Q1=="6 Closed Shrubland"),which(d$MCD12Q1=="4 Deciduous Broadleaf Forest"),which(d$MCD12Q1=="5 Mixed Forest")),]
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==6)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==4)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==5)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==8)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==0)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==11)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==12)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==13)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==14)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==15)]<- NA

  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==254)]<- NA
  MCD12Q1_TS_2001_2018[[i]][which(  MCD12Q1_TS_2001_2018[[i]][]==255)]<- NA

  MCD_masc <-  MCD_masc * MCD12Q1_TS_2001_2018[[i]]
  plot(MCD12Q1_TS_2001_2018[[i]]);plot(area.shp,add=T) #show both, if it works, good ;-)
  MCD_masc[which(  MCD_masc[]>1)]<- 1
}
MCD_masc[which(  MCD_masc[]>1)]<- 1


par(mfrow=c(1,1))
plot(MCD_masc,add=F);plot(area.shp,add=T) #show both, if it works, good ;-)


##show the data for an example year (2018)
plot(MCD12Q1_TS_2001_2018[[1]],col = terrain.colors(10));plot(area.shp,add=T)#;plot(osm.shp,add=T) #show both, if it works, good ;-)


###take the OSM street data and a 1km buffer and clip the extent from the masc:
osm.masc[which(  osm.masc[]==1)]<- 0;
osm.masc[which(  is.na(osm.masc[]) == T  )]<- 1
osm.masc[which(  osm.masc[]==0)]<- NA;

osm.masc<- mask(osm.masc,area.shp)
plot(osm.masc)

osm.masc <- resample(osm.masc, MCD_masc, method="ngb")
par(mfrow=c(1,3))
plot(MCD_masc)
plot(osm.masc)
plot(MCD_masc*osm.masc)



cellStats(MCD_masc, 'sum') #how many data points we have?
cellStats(MCD_masc*osm.masc, 'sum') #how many data points we have?

#create complete masc and safe it
masc_all<- MCD_masc*osm.masc
plot(masc_all)

oldwd<-getwd()
setwd("../Model/00_Results")
writeRaster(masc_all,"01_20200823_masc_all.tif",overwrite=T)

setwd(oldwd)

###majority decision for land use classification:

#mode of a vector:
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

MCD12Q1_TS <- calc(MCD12Q1_TS_2001_2018, getmode)

par(mfrow= c(2,2))
plot(MCD12Q1_TS)
barplot(MCD12Q1_TS)

MCD12Q1_TS <- MCD12Q1_TS * masc_all
plot(MCD12Q1_TS)
barplot(MCD12Q1_TS)
par(mfrow= c(1,1))
plot(MCD12Q1_TS,col=terrain.colors(6))


#example of a single datapoint:
MCD12Q1_TS_2001_2018[351,600]
getmode(MCD12Q1_TS_2001_2018[351,600])

###
MCD12Q1_TS<- resample(MCD12Q1_TS,MOD13A3_EVI_2001_001,method= "ngb")
#MCD12Q1_TS<- MCD12Q1_TS * masc_all

names(MCD12Q1_TS)<- "MCD12Q1"
plot(MCD12Q1_TS*masc_all*MCD12Q1_TS)
##################################################################################################################
##################################################################################################################
##MOD13A3

modistsp_file         <- "../01_data/02_MODIS/VI_Monthly_1Km_v6/EVI/MOD13A3_EVI_2001_001.tif"
MOD13A3_EVI_2001_001  <-raster(modistsp_file)
plot(MOD13A3_EVI_2001_001 );plot(area.shp,add=T) #show both

in_virtual_file       <- "../01_data/02_MODIS/VI_Monthly_1Km_v6/Time_Series/RData/Terra/EVI/MOD13A3_EVI_1_2001_335_2018_RData.RData"
MOD13A3_EVI_TS  <-get(load(in_virtual_file))
#plot(MOD13A3_EVI_TS)#;plot(area.shp,add=T) #show both

#SDS Name 	Description 	Units 	Data Type 	Fill Value 	No Data Value 	Valid Range 	Scale Factor
#1 km monthly EVI 	1 km monthly EVI 	EVI 	16-bit signed integer 	-3000 	N/A 	-2000 to 10000 	0.0001

#create a vektor that aims at a month for each year
targets <- seq(1,216,12)

#create sublayer dummy for th later means:
MOD13A3_EVI_TS_mean<- MOD13A3_EVI_TS[[1:12]]



# resample the masc with natural neighbors
masc_all    <- resample(masc_all, MOD13A3_EVI_2001_001, method="ngb")
plot(masc_all)

#process every image individually:
for(i in 1:216){
  MOD13A3_EVI_TS[[i]]<- MOD13A3_EVI_TS[[i]] * 	masc_all
  MOD13A3_EVI_TS[[i]][which(  MOD13A3_EVI_TS[[i]][]==-3000)]<- NA
  MOD13A3_EVI_TS[[i]][which(  MOD13A3_EVI_TS[[i]][]==0)]<- NA
  MOD13A3_EVI_TS[[i]]<- MOD13A3_EVI_TS[[i]] * 	0.0001
  print(i)
}
#show it
plot(MOD13A3_EVI_TS$layer.1)


# show a month for each year
plot(  MOD13A3_EVI_TS[[targets+1]])
##

#take the average for every year
for(i in 0:11){
  MOD13A3_EVI_TS_mean[[i+1]]<- stackApply(MOD13A3_EVI_TS[[targets+i]],indices=1,fun=mean)
  print(i)
}



MOD13A3_EVI_TS_mean[[1]][320,200]==mean(MOD13A3_EVI_TS[[targets]][320,200])
pdf(file= "00_evi.pdf",width=20 ,height=20,bg="transparent",family="Helvetica")#,paper="a4")
plot(MOD13A3_EVI_TS_mean,main="",xlab="",legend=F,ylab="",axes=F,box=F)

dev.off()


yr <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
names(MOD13A3_EVI_TS_mean)<-paste("evi",yr,sep="_")

##################################################################################################################
##################################################################################################################
##TRMM

trmmpath<- "../01_data/03_TRMM/"
setwd(trmmpath)
get(load("TRMM_annual.RData"))
plot(TRMM_annual)
names(TRMM_annual)
TRMM_annual    <- TRMM_annual[[-4]]
TRMM_annual    <- TRMM_annual[[c(6,1,4,7,10,3,5,11,12,8,9,2)]]
TRMM_annual    <-  projectRaster(TRMM_annual, MOD13A3_EVI_2001_001)

TRMM_annual    <- resample(TRMM_annual, MOD13A3_EVI_2001_001, method="bilinear")
TRMM_annual_pr <- TRMM_annual
TRMM_annual    <- TRMM_annual * masc_all
names(TRMM_annual_pr)<-paste("trmm",yr,sep="_")
names(TRMM_annual)<-paste("trmm",yr,sep="_")
par(mfrow=c(1,2))
plot(TRMM_annual)
plot(TRMM_annual_pr)

pdf(file= "00_evi.pdf",width=10 ,height=10,bg="transparent",family="Helvetica")#,paper="a4")
plot(TRMM_annual_pr,main="",xlab="",legend=F,ylab="",axes=F,box=F,col= colorRampPalette(c("white","blue"))(5000))
dev.off()
##################################################################################################################

##################################################################################################################
##GPM

gpmpath<- "../01_data/03_GPM_TS/"
setwd(gpmpath)
gpm_f<- list.files(gpmpath,pattern=".tif")
d <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)
#gpm<-brick()
#gpm[[1]]    <-  projectRaster(  gpm[[1]], MOD13A3_EVI_2001_001)
par(mfrow=c(1,2))
for (i in 1:12){
  ra   <-  raster(gpm_f[i])
  ra   <-  ra* 24 * d[i]
  plot(ra)
  ra   <-  projectRaster(  ra, MOD13A3_EVI_2001_001)
  ra   <-  mask(ra,area.shp)
  ra   <-  crop(ra,area.shp)
  ra   <-  resample(  ra, MOD13A3_EVI_2001_001, method="bilinear")
  plot(ra)
  ra2<- ra* masc_all
  if (i==1){gpm<-brick(ra2)} else{  gpm[[i]]<-ra2 }
  if (i==1){gpm_pr<-brick(ra)} else{  gpm_pr[[i]]<-ra }
}
plot(gpm)
plot(gpm_pr)

##################################################################################################################
##################################################################################################################
#GPM vs. TRMM

i=6
  plot(TRMM_annual_pr[[i]])
  plot(gpm[[i]])
  dif <- TRMM_annual_pr[[i]]-gpm[[i]]
  plot(dif)
cellStats(TRMM_annual_pr[[i]],stat="mean")
cellStats(gpm[[i]],stat="mean")

##################################################################################################################
##################################################################################################################
##SRTM

srtmpath<- "../01_data/01_SRTM/"
setwd(srtmpath)
SRTM <- raster("SRTM_aggf10.tif")
plot(SRTM,xlab="",legend=F,ylab="",axes=F,box=F)
#SRTM
#SRTM <- aggregate(SRTM,fact=10,fun=mean)

SRTM    <- projectRaster(SRTM, MOD13A3_EVI_2001_001)
SRTM    <- resample(SRTM, MOD13A3_EVI_2001_001, method="bilinear")
SRTM_pr <- SRTM
SRTM    <- SRTM * masc_all
par(mfrow=c(1,2))
plot(SRTM)

plot(SRTM_pr)
names(SRTM)<-"srtm"
names(SRTM_pr)<-"srtm"

##################################################################################################################
##################################################################################################################
##Potected Areas

PApath<- "C:/Users/Marcus/Desktop/PaleoVegetationModel/Model/01_data/04_NationalParcs/"
setwd(PApath)
# 1 = Wildlife Reserve
# 2 = National Parc
# 3 = Controlled Hunting Area
# 6 = Nothing

PA    <- raster("ETH_ProtectedAreas.tif")
PA    <-  projectRaster(PA, MOD13A3_EVI_2001_001)
PA    <- resample(PA, MOD13A3_EVI_2001_001, method="ngb")
PA[which(  is.na(PA[]) == T  )]<- 6

PA    <- mask(PA,area.shp)
PA    <- PA * masc_all
plot(PA);plot(area.shp,add=T)
names(PA)<-"pa"

##################################################################################################################
##################################################################################################################
##TREE HEIGHT

PApath<- "../01_data/06_TreeHeight/"
setwd(PApath)

TH <- raster("Simard_Pinto_3DGlobalVeg_L3C.tif")
TH    <- projectRaster(TH, MOD13A3_EVI_2001_001)

plot(TH);plot(area.shp,add=T)

TH    <- mask(TH,area.shp)
TH    <- TH * masc_all
plot(TH)
##################################################################################################################
##################################################################################################################
##VEG COVER

#load the files:
in_virtual_file       <- "../01_data/07_VegCoverPerc/Veg_Cont_Fields_Yearly_250m_v6/Time_Series/RData/Terra/Perc_NonTree/MOD44B_Perc_NonTree_65_2001_65_2018_RData.RData"
MOD44B_NTC  <-get(load(in_virtual_file))
plot(MOD44B_NTC)

in_virtual_file       <- "../01_data/07_VegCoverPerc/Veg_Cont_Fields_Yearly_250m_v6/Time_Series/RData/Terra/Perc_NonVeg/MOD44B_Perc_NonVeg_65_2001_65_2018_RData.RData"
MOD44B_NVC  <-get(load(in_virtual_file))
plot(MOD44B_NVC)

in_virtual_file       <- "../01_data/07_VegCoverPerc/Veg_Cont_Fields_Yearly_250m_v6/Time_Series/RData/Terra/Perc_TreeCov/MOD44B_Perc_TreeCov_65_2001_65_2018_RData.RData"
MOD44B_TC  <-get(load(in_virtual_file))
plot(MOD44B_TC)

#mean along the time series:
MOD44B_NTC_m <- stackApply(MOD44B_NTC,indices=1,fun=mean)
plot(MOD44B_NTC_m);plot(area.shp,add=T)

MOD44B_NVC_m <- stackApply(MOD44B_NVC,indices=1,fun=mean)
plot(MOD44B_NVC_m) ;plot(area.shp,add=T)

MOD44B_TC_m <- stackApply(MOD44B_TC,indices=1,fun=mean)
plot(MOD44B_TC_m) ;plot(area.shp,add=T)

#resampe and masc:
MOD44B_NTC_m    <- resample(MOD44B_NTC_m, MOD13A3_EVI_2001_001, method="bilinear")
MOD44B_NTC_m    <- mask(MOD44B_NTC_m,area.shp)
MOD44B_NTC_m    <- MOD44B_NTC_m * masc_all
plot(MOD44B_NTC_m)
MOD44B_NVC_m    <- resample(MOD44B_NVC_m, MOD13A3_EVI_2001_001, method="bilinear")
MOD44B_NVC_m    <- mask(MOD44B_NVC_m,area.shp)
MOD44B_NVC_m    <- MOD44B_NVC_m * masc_all
plot(MOD44B_NVC_m)

MOD44B_TC_m    <- resample(MOD44B_TC_m, MOD13A3_EVI_2001_001, method="bilinear")
MOD44B_TC_m    <- mask(MOD44B_TC_m,area.shp)
MOD44B_TC_m    <- MOD44B_TC_m * masc_all
plot(MOD44B_TC_m)

##################################################################################################################
##################################################################################################################
##TPI

srtmpath<- "../01_data/01_SRTM/"
setwd(srtmpath)
TPI_srtm <- raster("SRTM_aggf10.tif")
#TPI <- terrain(TPI_srtm, opt="TPI", unit="radians")
#plot(TPI,xlab="",legend=F,ylab="",axes=F,box=F)
require(spatialEco)
TPI <- tpi(TPI_srtm, scale = 0.01, win = "circle", normalize = FALSE, zero.correct = T)
par(mfrow=c(1,1))
plot(TPI)
plot(TPI_srtm)

#plot(TPI_s3)
#writeRaster(TPI_s3,"TPIs3.tif")
#SRTM
#SRTM <- aggregate(SRTM,fact=10,fun=mean)

TPI    <- projectRaster(TPI, MOD13A3_EVI_2001_001)
TPI    <- resample(TPI, MOD13A3_EVI_2001_001, method="bilinear")
TPI_pr <- TPI
TPI    <- TPI * masc_all
par(mfrow=c(1,2))
plot(TPI)

plot(TPI_pr)
names(TPI)<-"tpi"
names(TPI_pr)<-"tpi"


##################################################################################################################
##################################################################################################################
##stack all together and save them to 2nd level
data_all <- stack(MCD12Q1_TS,MOD13A3_EVI_TS_mean, gpm,SRTM,PA ,TH,MOD44B_NTC_m,MOD44B_NVC_m,MOD44B_TC_m,TPI )


####only targets with no NA in all layers:
masc_new <- data_all[[1]]
masc_new[which(  masc_new[]>1)]<- 1
cellStats(masc_new, 'sum') #how many data points we have?

for (i in 1:length(data_all@layers))
{
  masc_new<- masc_new * data_all[[i]]
}
masc_new[which(  is.na(masc_new[])==F)]<- 1
masc_new[which(  is.na(masc_new[])==T)]<- NA
cellStats(masc_new, 'sum') #how many data points we have?
plot(masc_new)
for (i in 1:length(data_all@layers))
{
  data_all[[i]]<- masc_new * data_all[[i]]
}
names(data_all)<- c( "lc",  "evi_jan","evi_feb" , "evi_mar","evi_apr","evi_may" ,"evi_jun","evi_jul",     
                     "evi_aug","evi_sep","evi_oct","evi_nov" ,"evi_dec" ,
                     "gpm_jan","gpm_feb","gpm_mar","gpm_apr","gpm_may","gpm_jun","gpm_jul","gpm_aug","gpm_sep",
                     "gpm_oct","gpm_nov","gpm_dec","srtm","pa","th","ntc","nvc","tc","tpi")
####
data_all_pr <- stack(gpm_pr,SRTM_pr,TPI_pr)
names(data_all_pr)<-c( "gpm_jan","gpm_feb","gpm_mar","gpm_apr","gpm_may","gpm_jun","gpm_jul","gpm_aug","gpm_sep",
                       "gpm_oct","gpm_nov","gpm_dec","srtm","tpi")
path_2ndlevel<- "../02_data_2nd_level/"
setwd(path_2ndlevel)

test<- sum(data_all)
plot(test)

save(data_all, file="data_all.RData")
save(data_all_pr, file="data_all_pr.RData")
save(masc_all, file="mask_all.RData")

plot(data_all)


oldwd<-getwd()
setwd("../00_Results")
writeRaster(data_all$lc,"01_20200823_data_tr_lc.tif",overwrite=T)

setwd(oldwd)
