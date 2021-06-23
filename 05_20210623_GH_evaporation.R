# MLFischer 2020-09
# evaporation parametrization
remove(list=ls()) #clean up#
removeTmpFiles(0)

#get some packages
library(raster)
library(rgdal)  
library(gbm)
library(dismo)
library(rasterVis)
library(dplyr)
library(shapefiles)
require(ggplot2)

path_source<- "../00_Script/"
setwd(path_source)
source("05_20201109_ET_model_function.R")

path_4ndlevel<- "../02_data_4nd_scenarios/"
path_5ndlevel<- "../02_data_5nd_evaporation/"
path_6ndlevel<- "../02_data_6nd_model_link/"

setwd(path_5ndlevel)

para <- read.csv2("20201029_LC_Parametrization.csv")

#################################################################################################
#################################################################################################
#TIME SLICE 1a: MODERN DAY all LC


#################################################################################################
#First a plot of the modern day distribution
f_a <- read.csv2("01_20201105_A_all_LC_F.csv")[c(-12),]
f_b <- read.csv2("01_20201105_B_all_LC_F.csv")[c(-10),]
f_c <- read.csv2("01_20201105_C_all_LC_F.csv")[c(-9),]

f_a$perc <-f_a[,3]/ sum(f_a[,3])
f_b$perc <-f_b[,3]/ sum(f_b[,3])
f_c$perc <-f_c[,3]/ sum(f_c[,3])

f_a$area <- "a"
f_b$area <- "b"
f_c$area <- "c"

f_a$value<- as.factor(f_a$value)
f_b$value<- as.factor(f_b$value)
f_c$value<- as.factor(f_c$value)

f<-rbind(f_a,f_b,f_c)
f_MD <- f
levels(f$value)

MD_dist <- ggplot(f, aes(x="", y=perc, fill=value)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void() + facet_wrap(~area)
plot(MD_dist)
#################################################################################################



####
#ETa and para tuning:
f_a <- read.csv2("01_20201105_A_all_LC_F.csv")[c(-1,-12),]
f_b <- read.csv2("01_20201105_B_all_LC_F.csv")[c(-1,-10),]
f_c <- read.csv2("01_20201105_C_all_LC_F.csv")[c(-1,-9),]

#percentage coverage calculation without lake areas:
f_a$perc <-f_a[,3]/ sum(f_a[,3])
f_b$perc <-f_b[,3]/ sum(f_b[,3])
f_c$perc <-f_c[,3]/ sum(f_c[,3])


#albedo:
f_a$alb <- para[(f_a$value+1),2] * f_a$perc
f_b$alb <- para[(f_b$value+1),2] * f_b$perc
f_c$alb <- para[(f_c$value+1),2] * f_c$perc

#albedo catchment average above land with no tuning
sum(f_a$alb) # 0.136
sum(f_b$alb) # 0.14 
sum(f_c$alb) # 0.126

f_a$alb <- para[(f_a$value+1),10] * f_a$perc
f_b$alb <- para[(f_b$value+1),10] * f_b$perc
f_c$alb <- para[(f_c$value+1),10] * f_c$perc

#albedo catchment average above land: after tuning:
sum(f_a$alb) # 0.136
sum(f_b$alb) # 0.14 
sum(f_c$alb) # 0.126


#emissivity:
f_a$em <- para[(f_a$value+1),3] * f_a$perc
f_b$em <- para[(f_b$value+1),3] * f_b$perc
f_c$em <- para[(f_c$value+1),3] * f_c$perc

sum(f_a$em) # 0.96
sum(f_b$em) # 0.98
sum(f_c$em) # 0.96


#soil moisture availability not tuned:
f_a$sma <- para[(f_a$value+1),4] * f_a$perc
f_b$sma <- para[(f_b$value+1),4] * f_b$perc
f_c$sma <- para[(f_c$value+1),4] * f_c$perc

sum(f_a$sma) # 0.25
sum(f_b$sma) # 0.137
sum(f_c$sma) # 0.198

#soil moisture availability tuned:
f_a$sma <- para[(f_a$value+1),11] * f_a$perc
f_b$sma <- para[(f_b$value+1),11] * f_b$perc
f_c$sma <- para[(f_c$value+1),11] * f_c$perc

sum(f_a$sma) # 0.25
sum(f_b$sma) # 0.137
sum(f_c$sma) # 0.198


###Surface Drag coefficient:
f_a$z0 <- para[(f_a$value+1),5] * f_a$perc
f_b$z0 <- para[(f_b$value+1),5] * f_b$perc
f_c$z0 <- para[(f_c$value+1),5] * f_c$perc

###Surface Drag coefficient tuned:
f_a$z0 <- para[(f_a$value+1),12] * f_a$perc
f_b$z0 <- para[(f_b$value+1),12] * f_b$perc
f_c$z0 <- para[(f_c$value+1),12] * f_c$perc

#roughness length: 
z0_a <- sum(f_a$z0)/100 #(m)
z0_b <- sum(f_b$z0)/100 
z0_c <- sum(f_c$z0)/100 

#geostrophic wind: 
#u_z_a = 1.42 * (log(1:500/z0_a))/(log(10/z0_a))
#u_z_b = 2.13 * (log(1:500/z0_b))/(log(10/z0_b))
#u_z_c = 1.72 * (log(1:500/z0_c))/(log(10/z0_c))
#logarithmic wind profile:
#plot(y=1:500,x= u_z_a,type="l",xlab="Windsped",ylab="height (m)",xlim=c(0,5))
#lines(y=1:500,x= u_z_b,col="red")
#lines(y=1:500,x= u_z_c,col="blue")

#Rowntree In; Schmugge und Andre p.9f, Bookhagen 2001, Bergner 2003

#z0_a<- 0.5
u_a      <- (0.4*1.42)/ log(10/z0_a) #friction velocity
r_a_a    <- log(10/z0_a)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0_a))^-2) #surface drag coefficient
cd_a*1.4 # # 0.0076

u_b      <- (0.4*2.13)/ log(10/z0_b)
r_a_b    <- log(10/z0_b)/(u_b * 0.4)
cd_b     <- 0.4*0.4 * (log(r_a_b /(z0_b))^-2)
cd_b # 0.0059

u_c      <- (0.4*1.72)/ log(10/z0_c)
r_a_c    <- log(10/z0_c)/(u_c * 0.4)
cd_c     <- 0.4*0.4 * (log(r_a_c /(z0_c))^-2)
cd_c*1.2# 0.0066


# check all single class surface drag coefficients: 
#z0_c <- (para$roughness.length.tuned/100)
#u_c      <- (0.4*1.72)/ log(10/z0_c)
#r_a_c    <- log(10/z0_c)/(u_c * 0.4)
#cd     <- 0.4*0.4 * (log(r_a_c /(z0_c))^-2)
#cd# 0.0066
#para$cd_c <- cd
#para[,c(7,5,12,14)]


#Abaya ETa 1,123 mm/a
ETa_a_MD <- ETa( t_a_ld    = 291.5, 
  emis_ld   = sum(f_a$em), 
  albedo_ld = sum(f_a$alb),
  rh        = 0.58,
  f_ld      = sum(f_a$sma),
  cc        = 0.54,
  ws        = 1.42,
  a         = 0.39,
  b         = 0.38,
  a_2       = 0.22,
  b_2       = 2.00,
  cds       = cd_a*1.4,
  p         = 81735,
  r_swc     = 415)
ETa_a_MD


#Bahir ETa 892 mm/a
ETa_b_MD <- ETa( t_a_ld    = 295.9, 
     emis_ld   = sum(f_b$em), 
     albedo_ld = sum(f_b$alb),
     rh        = 0.57,
     f_ld      = sum(f_b$sma),
     cc        = 0.61,
     ws        = 2.13,
     a         = 0.39,
     b         = 0.38,
     a_2       = 0.22,
     b_2       = 2.00,
     cds       = cd_b,
     p         = 88466.63,
     r_swc     = 415)
ETa_b_MD


#Chamo ETa 1060 mm/a 
ETa_c_MD <- ETa( t_a_ld    = 293.3, 
     emis_ld   = sum(f_c$em), 
     albedo_ld = sum(f_c$alb),
     rh        = 0.57,
     f_ld      = sum(f_c$sma),
     cc        = 0.55,
     ws        = 1.72,
     a         = 0.39,
     b         = 0.38,
     a_2       = 0.22,
     b_2       = 2.00,
     cds       = cd_c*1.2,
     p         = 85128,
     r_swc     = 415)#
ETa_c_MD
  
#comparisons actual parametrization and results from Fischer et al., 2020
  ETa_a_MD -  1123
  ETa_b_MD -   892
  ETa_c_MD -  1060
  
  
#################################################################################################
#################################################################################################  
#################################################################################################
#################################################################################################
#TIME SLICE 1b: MODERN DAY potential vegetation or Pre-Industrial: 
  setwd("../01_data/00_Shapes/")
  #area.shp <- readOGR(".",layer = "20190902_masc_studyarea_stage1")
  #a.shp <- readOGR(".",layer = "20191024_abaya")
  #b.shp <- readOGR(".",layer = "20191024_bahir")
  #c.shp <- readOGR(".",layer = "20191024_chamo")
  
  #save(area.shp,file= "20201119_masc_studyarea_stage1.RData")
  #save(a.shp,file= "20201119__abaya.RData")
  #save(b.shp,file= "20201119__bahir.RData")
  #save(c.shp,file= "20201119__chamo.RData")
  
  get(load("20201119_masc_studyarea_stage1.RData"))
  get(load("20201119__abaya.RData"))
  get(load("20201119__bahir.RData"))
  get(load("20201119__chamo.RData"))
  
  setwd(path_4ndlevel)
  r_pi <- raster("04_20201117_PI_pr_all.tif")
  r_md <- raster("01_20201116_LC_MD_all.tif")
  r_md <- projectRaster(r_md,crs= a.shp@proj4string,method="ngb" )
  
  r_pi[which(r_md[]==0)]<-0
  r_pi[which(r_md[]==15)]<-15
  
  par(mfrow=c(2,2))
  plot(r_pi)
  plot(area.shp,add=T)

  r3_a <- mask(r_pi,a.shp)
  plot(r3_a)
  f_a <- as.data.frame(freq(r3_a))[c(-6),]
  
  r3_b <- mask(r_pi,b.shp)
  plot(r3_b)
  f_b <- as.data.frame(freq(r3_b))[c(-6),]
  
  r3_c <- mask(r_pi,c.shp)
  plot(r3_c)
  f_c <- as.data.frame(freq(r3_c))[c(-6),]
  
  f_a$perc <-f_a$count/ sum(f_a$count)
  f_b$perc <-f_b$count/ sum(f_b$count)
  f_c$perc <-f_c$count/ sum(f_c$count)
  
  f_a$area <- "a"
  f_b$area <- "b"
  f_c$area <- "c"
  
  f_a$value<- as.factor(f_a$value)
  f_b$value<- as.factor(f_b$value)
  f_c$value<- as.factor(f_c$value)
  
  f   <-rbind(f_a,f_b,f_c)
  f_PI<- f
  
  levels(f$value)
  
  PI_dist <- ggplot(f, aes(x="", y=perc, fill=value)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0)+
    theme_void() + facet_wrap(~area)
  
  plot(PI_dist)
  
  
  #############
  f_a <- as.data.frame(freq(r3_a))[c(-1,-6),]
  f_b <- as.data.frame(freq(r3_b))[c(-6),]
  f_c <- as.data.frame(freq(r3_c))[c(-1,-6),]
  
  f_a$perc <-f_a$count/ sum(f_a$count)
  f_b$perc <-f_b$count/ sum(f_b$count)
  f_c$perc <-f_c$count/ sum(f_c$count)
  
  f_a$alb <- para[(f_a$value+1),10] * f_a$perc
  f_b$alb <- para[(f_b$value+1),10] * f_b$perc
  f_c$alb <- para[(f_c$value+1),10] * f_c$perc
  
  #albedo catchment average above land: after tuning:
  sum(f_a$alb) # 0.136
  sum(f_b$alb) # 0.14 
  sum(f_c$alb) # 0.126
  
  
  #emissivity:
  f_a$em <- para[(f_a$value+1),3] * f_a$perc
  f_b$em <- para[(f_b$value+1),3] * f_b$perc
  f_c$em <- para[(f_c$value+1),3] * f_c$perc
  
  sum(f_a$em) # 0.96
  sum(f_b$em) # 0.98
  sum(f_c$em) # 0.96
  
  #soil moisture availability tuned:
  f_a$sma <- para[(f_a$value+1),11] * f_a$perc
  f_b$sma <- para[(f_b$value+1),11] * f_b$perc
  f_c$sma <- para[(f_c$value+1),11] * f_c$perc
  
  sum(f_a$sma) # 0.25
  sum(f_b$sma) # 0.137
  sum(f_c$sma) # 0.198
  
  ###Surface Drag coefficient tuned:
  f_a$z0 <- para[(f_a$value+1),12] * f_a$perc
  f_b$z0 <- para[(f_b$value+1),12] * f_b$perc
  f_c$z0 <- para[(f_c$value+1),12] * f_c$perc
  
  #roughness length: 
  z0_a <- sum(f_a$z0)/100 #(m)
  z0_b <- sum(f_b$z0)/100 
  z0_c <- sum(f_c$z0)/100 
  
  #geostrophic wind: 
  #u_z_a = 1.42 * (log(1:500/z0_a))/(log(10/z0_a))
  #u_z_b = 2.13 * (log(1:500/z0_b))/(log(10/z0_b))
  #u_z_c = 1.72 * (log(1:500/z0_c))/(log(10/z0_c))
  #logarithmic wind profile:
  #plot(y=1:500,x= u_z_a,type="l",xlab="Windsped",ylab="height (m)",xlim=c(0,5))
  #lines(y=1:500,x= u_z_b,col="red")
  #lines(y=1:500,x= u_z_c,col="blue")
  
  #Rowntree In; Schmugge und Andre p.9f, Bookhagen 2001, Bergner 2003
  
  #z0_a<- 0.5
  u_a      <- (0.4*1.42)/ log(10/z0_a) #friction velocity
  r_a_a    <- log(10/z0_a)/(u_a * 0.4) #atmospheric resitance
  cd_a     <- 0.4*0.4 * (log(r_a_a /(z0_a))^-2) #surface drag coefficient
  cd_a # # 0.0076
  
  u_b      <- (0.4*2.13)/ log(10/z0_b)
  r_a_b    <- log(10/z0_b)/(u_b * 0.4)
  cd_b     <- 0.4*0.4 * (log(r_a_b /(z0_b))^-2)
  cd_b # 0.0059
  
  u_c      <- (0.4*1.72)/ log(10/z0_c)
  r_a_c    <- log(10/z0_c)/(u_c * 0.4)
  cd_c     <- 0.4*0.4 * (log(r_a_c /(z0_c))^-2)
  cd_c# 0.0066
  
  #Abaya ETa 1,123 mm/a
  ETa_a_PI <- ETa( t_a_ld    = 291.5, 
                  emis_ld   = sum(f_a$em), 
                  albedo_ld = sum(f_a$alb),
                  rh        = 0.58,
                  f_ld      = sum(f_a$sma),
                  cc        = 0.54,
                  ws        = 1.42,
                  a         = 0.39,
                  b         = 0.38,
                  a_2       = 0.22,
                  b_2       = 2.00,
                  cds       = cd_a*1.4,
                  p         = 81735,
                  r_swc     = 415)
  ETa_a_PI
  
  
  #Bahir ETa 892 mm/a
  ETa_b_PI <- ETa( t_a_ld    = 295.9, 
                  emis_ld   = sum(f_b$em), 
                  albedo_ld = sum(f_b$alb),
                  rh        = 0.57,
                  f_ld      = sum(f_b$sma),
                  cc        = 0.61,
                  ws        = 2.13,
                  a         = 0.39,
                  b         = 0.38,
                  a_2       = 0.22,
                  b_2       = 2.00,
                  cds       = cd_b,
                  p         = 88466.63,
                  r_swc     = 415)
  ETa_b_PI
  
  
  #Chamo ETa 1060 mm/a 
  ETa_c_PI <- ETa( t_a_ld    = 293.3, 
                  emis_ld   = sum(f_c$em), 
                  albedo_ld = sum(f_c$alb),
                  rh        = 0.57,
                  f_ld      = sum(f_c$sma),
                  cc        = 0.55,
                  ws        = 1.72,
                  a         = 0.39,
                  b         = 0.38,
                  a_2       = 0.22,
                  b_2       = 2.00,
                  cds       = cd_c*1.2,
                  p         = 85128,
                  r_swc     = 415)#
  ETa_c_PI
  
  #comparisons Modern Day and Pre Industrial Evaporation rates
  ETa_a_MD -  1123
  ETa_b_MD -   892
  ETa_c_MD -  1060
  
  ETa_a_PI -  1123
  ETa_b_PI -   892
  ETa_c_PI -  1060
  
  (ETa_a_MD -  ETa_a_PI)/ETa_a_MD * 100
  (ETa_b_MD -  ETa_b_PI)/ETa_b_MD * 100
  (ETa_c_MD -  ETa_c_PI)/ETa_c_MD * 100
  
  
  
  
  
  
  #################################################################################################
  #################################################################################################
  # Time Slice 2a: AHP 
  # NO CAB
  #
  xlim_nr<- 140
  
  setwd(path_4ndlevel)
  get(load("20201120_AHP01_noCAB_LC.RData"))  
  
  d2a <- df_s[,1:13]
  d2a$pre <- d2a$pre*100
  
  #albedo for each catchment
  p_n <- 10
  d2a$alb_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
  d2a$alb_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
  d2a$alb_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]

  
  #emissivity for each catchment
  p_n <- 3
  d2a$em_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
  d2a$em_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
  d2a$em_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]
  
  #soil moisture availabiltiy for each catchment
  p_n <- 11
  d2a$sma_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
  d2a$sma_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
  d2a$sma_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]
  
  
  ###Surface drag coefficient for each catchment:
  
  #roughness length:
  p_n <- 12
  d2a$z0_a <- (d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n])/100
  d2a$z0_b <- (d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n])/100
  d2a$z0_c <- (d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n])/100
  
  #friction velocity:
  d2a$u_a <- (0.4*1.42)/ log(10/d2a$z0_a)
  d2a$u_b <- (0.4*2.13)/ log(10/d2a$z0_b)
  d2a$u_c <- (0.4*1.72)/ log(10/d2a$z0_c)
  
  #atmospheric resitance:
  d2a$r_a_a    <- log(10/d2a$z0_a)/(d2a$u_a * 0.4) #atmospheric resitance
  d2a$r_a_b    <- log(10/d2a$z0_b)/(d2a$u_b * 0.4) #atmospheric resitance
  d2a$r_a_c    <- log(10/d2a$z0_c)/(d2a$u_c * 0.4) #atmospheric resitance
  
  #surface drag coefficient
  d2a$cd_a     <- 0.4*0.4 * (log(d2a$r_a_a /(d2a$z0_a))^-2) #surface drag coefficient
  d2a$cd_b     <- 0.4*0.4 * (log(d2a$r_a_b /(d2a$z0_b))^-2) #surface drag coefficient
  d2a$cd_c     <- 0.4*0.4 * (log(d2a$r_a_c /(d2a$z0_c))^-2) #surface drag coefficient

  d2a$ETa_a<-NULL
  d2a$ETa_b<-NULL
  d2a$ETa_c<-NULL
  
  d2a$ETa_a_u<-NULL
  d2a$ETa_b_u<-NULL
  d2a$ETa_c_u<-NULL
  
  d2a$ETa_a_l<-NULL
  d2a$ETa_b_l<-NULL
  d2a$ETa_c_l<-NULL
  
  for(i in 1:101){
  d2a$ETa_a[i] <- ETa(    t_a_ld    = 291.5, 
                   emis_ld   = d2a$em_a[i], 
                   albedo_ld = d2a$alb_a[i],
                   rh        = 0.58,
                   f_ld      = d2a$sma_a[i],
                   cc        = 0.54,
                   ws        = 1.42,
                   a         = 0.39,
                   b         = 0.38,
                   a_2       = 0.22,
                   b_2       = 2.00,
                   cds       = d2a$cd_a[i]*1.4,
                   p         = 81735,
                   r_swc     = 415)
  
  d2a$ETa_b[i] <- ETa( t_a_ld    = 295.9, 
                       emis_ld   = d2a$em_b[i],
                       albedo_ld = d2a$alb_b[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_b[i],
                       cc        = 0.61,
                       ws        = 2.13,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = d2a$cd_b[i],
                       p         = 88466.63,
                       r_swc     = 415)
  
  d2a$ETa_c[i] <- ETa( t_a_ld    = 293.3, 
                       emis_ld   = d2a$em_c[i],
                       albedo_ld = d2a$alb_c[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_c[i],
                       cc        = 0.55,
                       ws        = 1.72,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = d2a$cd_c[i]*1.2,
                       p         = 85128,
                       r_swc     = 415)#

   }
  
  #comparisons actual parametrization and results from Fischer et al., 2020
  ETa_a_MD -  1123
  ETa_b_MD -   892
  ETa_c_MD -  1060  
  
  d2a$pre[51]
  
  ETa_a_PI -   d2a$ETa_a[51]
  ETa_b_PI -   d2a$ETa_b[51]
  ETa_c_PI -   d2a$ETa_c[51]
  
  
  
 print( paste(sum(f_a$alb) , d2a$alb_a[51], "0.136"))
 print( paste(sum(f_b$alb) , d2a$alb_b[51], "0.14")) 
 print( paste(sum(f_c$alb) , d2a$alb_c[51], "0.126"))
  
 print( paste(sum(f_a$sma) , d2a$sma_a[51], "0.25"))
 print( paste(sum(f_b$sma) , d2a$sma_b[51], "0.137")) 
 print( paste(sum(f_c$sma) , d2a$sma_c[51], "0.198")) 
  
 print( paste(cd_a , d2a$cd_a[51], "0.0076"))
 print( paste(cd_b , d2a$cd_b[51], "0.0059")) 
 print( paste(cd_c , d2a$cd_c[51], "0.0066")) 
  
  d2a$ETa_a[51]
  d2a$ETa_b[51]
  d2a$ETa_c[51]
  
  d2a$ETa_a_p <- (d2a$ETa_a/ETa_a_MD )*100
  d2a$ETa_b_p <- (d2a$ETa_b/ETa_b_MD )*100
  d2a$ETa_c_p <- (d2a$ETa_c/ETa_c_MD )*100
  
  
  #albedo catchment average above land: after tuning:
  sum(f_a$alb) # 0.136
  sum(f_b$alb) # 0.14 
  sum(f_c$alb) # 0.126
  
  sum(f_b$alb) # 0.14 
  d2a$alb_b[51]
  
  d2a$cd_b[51]
  
  sum(f_b$sma) # 0.14 
  d2a$sma_b[51]
  
  d2a$sma_a[51]
  d2a$alb_b[51]
  d2a$ETa_c[51]
  
  
 f_PI[1:5,c(1,3)]
 d2a[51,2:5]
  
 f_PI[6:10,c(1,3)]
 d2a[51,6:9]

 
 d2a_real <- d2a
 d2a<- d2a[51:91,]
 
 
par(mfrow=c(2,3))  

plot(d2a$pre, d2a$ETa_a,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Abaya",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_a_MD,col="darkgreen",pch=16)
points(100,1123,col="red",pch=16)

abline(h=1123,col="red" )
lines(d2a$pre,d2a$pre/100*1147,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1123),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_a_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_a_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)



plot(d2a$pre, d2a$ETa_c,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Chamo",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_c_MD,col="darkgreen",pch=16)
points(100,1060,col="red",pch=16)


abline(h=1060,col="red" )
lines(d2a$pre,d2a$pre/100*1131,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1060),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_c_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_c_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)


plot(d2a$pre, d2a$ETa_b,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Bahir",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_b_MD,col="darkgreen",pch=16)
points(100,892,col="red",pch=16)

abline(h=892,col="red" )
lines(d2a$pre,d2a$pre/100*893,col="blue")

mtext(text = paste("ETa_MD_orig=",round(892),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_b_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_b_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)

b_c <- c("#267300","#ff9500","#4ecc00","#ffff00")
p_r <- rep(100,4)

plot (d2a$pre, d2a$a_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
#lines(d2a$pre, d2a$a_2*100+d2a$a_7*100,col="#ff9500")

#polygon(c(d2a$pre, d2a$pre[101:1]), c(d2a$a_2*100, rep(0,101)),
#        border = NA,   col = "#267300")
#polygon(c(d2a$pre, d2a$pre[101:1]), c( (d2a$a_2*100+d2a$a_9*100) , 100*d2a$a_2[101:1]),
#        border = NA,   col = "#4ecc00")
lines(d2a$pre, d2a$a_2*100+d2a$a_9*100,col="#4ecc00")
lines(d2a$pre, d2a$a_2*100+d2a$a_9*100+d2a$a_10*100,col="#ffff00")
points(p_r ,100*f_MD[c(2,3,5,6),4],pch=16,col=b_c)

plot (d2a$pre, d2a$c_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
#lines(d2a$pre, d2a$c_7*100,col="#ff9500")
lines(d2a$pre, d2a$c_2*100+d2a$c_9*100,col="#4ecc00")
lines(d2a$pre, d2a$c_2*100+d2a$c_9*100+d2a$c_10*100,col="#ffff00")
points(p_r[1:3] ,100*f_MD[c(22,23,24),4],pch=16,col=b_c[c(1,3,4)])



plot (d2a$pre, d2a$b_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100,col="#4ecc00")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100+d2a$b_10*100,col="#ffff00")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100+d2a$b_10*100+d2a$b_7*100,col="#ff9500")

points(p_r ,100*f_MD[c(12,14,16,17),4],pch=16,col=b_c)

d2a_real -> d2a
setwd(path_6ndlevel)
save(d2a,file= "20201118_AHP01_noCAB_LC_ET.RData")
#################################################################################################
#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################
#################################################################################################

# Time Slice 2b: AHP 
# CAB
#

setwd(path_4ndlevel)
get(load("20201120_AHP01_CAB_LC.RData"))  

d2a <- df_s[,1:13]
d2a$pre <- d2a$pre*100

#albedo for each catchment
p_n <- 10
d2a$alb_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
d2a$alb_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
d2a$alb_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]


#emissivity for each catchment
p_n <- 3
d2a$em_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
d2a$em_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
d2a$em_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]

#soil moisture availabiltiy for each catchment
p_n <- 11
d2a$sma_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]
d2a$sma_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]
d2a$sma_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]


###Surface drag coefficient for each catchment:

#roughness length:
p_n <- 12
d2a$z0_a <- (d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n])/100
d2a$z0_b <- (d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n])/100
d2a$z0_c <- (d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n])/100

#friction velocity:
d2a$u_a <- (0.4*1.42)/ log(10/d2a$z0_a)
d2a$u_b <- (0.4*2.13)/ log(10/d2a$z0_b)
d2a$u_c <- (0.4*1.72)/ log(10/d2a$z0_c)

#atmospheric resitance:
d2a$r_a_a    <- log(10/d2a$z0_a)/(d2a$u_a * 0.4) #atmospheric resitance
d2a$r_a_b    <- log(10/d2a$z0_b)/(d2a$u_b * 0.4) #atmospheric resitance
d2a$r_a_c    <- log(10/d2a$z0_c)/(d2a$u_c * 0.4) #atmospheric resitance

#surface drag coefficient
d2a$cd_a     <- 0.4*0.4 * (log(d2a$r_a_a /(d2a$z0_a))^-2) #surface drag coefficient
d2a$cd_b     <- 0.4*0.4 * (log(d2a$r_a_b /(d2a$z0_b))^-2) #surface drag coefficient
d2a$cd_c     <- 0.4*0.4 * (log(d2a$r_a_c /(d2a$z0_c))^-2) #surface drag coefficient

d2a$ETa_a<-NULL
d2a$ETa_b<-NULL
d2a$ETa_c<-NULL

d2a$ETa_a_u<-NULL
d2a$ETa_b_u<-NULL
d2a$ETa_c_u<-NULL

d2a$ETa_a_l<-NULL
d2a$ETa_b_l<-NULL
d2a$ETa_c_l<-NULL

for(i in 1:101){
  d2a$ETa_a[i] <- ETa(    t_a_ld    = 291.5, 
                          emis_ld   = d2a$em_a[i], 
                          albedo_ld = d2a$alb_a[i],
                          rh        = 0.58,
                          f_ld      = d2a$sma_a[i],
                          cc        = 0.54,
                          ws        = 1.42,
                          a         = 0.39,
                          b         = 0.38,
                          a_2       = 0.22,
                          b_2       = 2.00,
                          cds       = d2a$cd_a[i]*1.4,
                          p         = 81735,
                          r_swc     = 415)
  
  d2a$ETa_b[i] <- ETa( t_a_ld    = 295.9, 
                       emis_ld   = d2a$em_b[i],
                       albedo_ld = d2a$alb_b[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_b[i],
                       cc        = 0.61,
                       ws        = 2.13,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = d2a$cd_b[i],
                       p         = 88466.63,
                       r_swc     = 415)
  
  d2a$ETa_c[i] <- ETa( t_a_ld    = 293.3, 
                       emis_ld   = d2a$em_c[i],
                       albedo_ld = d2a$alb_c[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_c[i],
                       cc        = 0.55,
                       ws        = 1.72,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00, 
                       cds       = d2a$cd_c[i]*1.2,
                       p         = 85128,
                       r_swc     = 415)#
  
}

#comparisons actual parametrization and results from Fischer et al., 2020
ETa_a_MD -  1123
ETa_b_MD -   892
ETa_c_MD -  1060  

d2a$pre[51]

ETa_a_PI -   d2a$ETa_a[51]
ETa_b_PI -   d2a$ETa_b[51]
ETa_c_PI -   d2a$ETa_c[51]

d2a_real_ahp <- d2a
d2a<- d2a[51:91,]

par(mfrow=c(2,3))  

plot(d2a$pre, d2a$ETa_a,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Abaya",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_a_MD,col="darkgreen",pch=16)
points(100,1123,col="red",pch=16)

abline(h=1123,col="red" )
lines(d2a$pre,d2a$pre/100*1147,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1123),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_a_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_a_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)



plot(d2a$pre, d2a$ETa_c,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Chamo",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_c_MD,col="darkgreen",pch=16)
points(100,1060,col="red",pch=16)


abline(h=1060,col="red" )
lines(d2a$pre,d2a$pre/100*1131,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1060),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_c_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_c_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)


plot(d2a$pre, d2a$ETa_b,type="l",xlim=c(100,xlim_nr),ylim=c(750,2000),col="black",
     main="Bahir",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_b_MD,col="darkgreen",pch=16)
points(100,892,col="red",pch=16)

abline(h=892,col="red" )
lines(d2a$pre,d2a$pre/100*893,col="blue")

mtext(text = paste("ETa_MD_orig=",round(892),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_b_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_b_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)

b_c <- c("#267300","#ff9500","#4ecc00","#ffff00")
p_r <- rep(100,4)

plot (d2a$pre, d2a$a_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
#lines(d2a$pre, d2a$a_7*100,col="#ff9500")

#polygon(c(d2a$pre, d2a$pre[101:1]), c(d2a$a_2*100, rep(0,101)),
#        border = NA,   col = "#267300")
#polygon(c(d2a$pre, d2a$pre[101:1]), c( (d2a$a_2*100+d2a$a_9*100) , 100*d2a$a_2[101:1]),
#        border = NA,   col = "#4ecc00")
lines(d2a$pre, d2a$a_2*100+d2a$a_9*100,col="#4ecc00")
lines(d2a$pre, d2a$a_2*100+d2a$a_9*100+d2a$a_10*100,col="#ffff00")
points(p_r ,100*f_MD[c(2,3,5,6),4],pch=16,col=b_c)

plot (d2a$pre, d2a$c_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
#lines(d2a$pre, d2a$c_7*100,col="#ff9500")
lines(d2a$pre, d2a$c_2*100+d2a$c_9*100,col="#4ecc00")
lines(d2a$pre, d2a$c_2*100+d2a$c_9*100+d2a$c_10*100,col="#ffff00")
points(p_r[1:3] ,100*f_MD[c(22,23,24),4],pch=16,col=b_c[c(1,3,4)])



plot (d2a$pre, d2a$b_2*100,type="l",xlim=c(100,xlim_nr),ylim=c(0,100),col="#267300")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100,col="#4ecc00")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100+d2a$b_10*100,col="#ffff00")
lines(d2a$pre, d2a$b_2*100+d2a$b_9*100+d2a$b_10*100+d2a$b_7*100,col="#ff9500")

points(p_r ,100*f_MD[c(12,14,16,17),4],pch=16,col=b_c)



d2a_real_ahp -> d2a

setwd(path_6ndlevel)
save(d2a,file= "20201125_AHP01_CAB_LC_ET.RData")

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
# Time Slice 3: LGM 
#LGM UPPER

setwd(path_4ndlevel)
d3 <- get(load("20201104_LGM01_TU_LC.RData"))  
d2a<-d3


d2a$pre <- d2a$pre*100

#albedo for each catchment
p_n <- 10
d2a$alb_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$alb_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$alb_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]


#emissivity for each catchment
p_n <- 3
d2a$em_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$em_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$em_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]

#soil moisture availabiltiy for each catchment
p_n <- 11
d2a$sma_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$sma_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$sma_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]


###Surface drag coefficient for each catchment:

#roughness length:
p_n <- 12
d2a$z0_a <- (d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n])/100
d2a$z0_b <- (d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n])/100
d2a$z0_c <- (d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n])/100

#friction velocity:
d2a$u_a <- (0.4*1.42)/ log(10/d2a$z0_a)
d2a$u_b <- (0.4*2.13)/ log(10/d2a$z0_b)
d2a$u_c <- (0.4*1.72)/ log(10/d2a$z0_c)

#atmospheric resitance:
d2a$r_a_a    <- log(10/d2a$z0_a)/(d2a$u_a * 0.4) #atmospheric resitance
d2a$r_a_b    <- log(10/d2a$z0_b)/(d2a$u_b * 0.4) #atmospheric resitance
d2a$r_a_c    <- log(10/d2a$z0_c)/(d2a$u_c * 0.4) #atmospheric resitance

#surface drag coefficient
d2a$cd_a     <- 0.4*0.4 * (log(d2a$r_a_a /(d2a$z0_a))^-2) #surface drag coefficient
d2a$cd_b     <- 0.4*0.4 * (log(d2a$r_a_b /(d2a$z0_b))^-2) #surface drag coefficient
d2a$cd_c     <- 0.4*0.4 * (log(d2a$r_a_c /(d2a$z0_c))^-2) #surface drag coefficient

d2a$ETa_a_u<-NULL
d2a$ETa_b_u<-NULL
d2a$ETa_c_u<-NULL

d2a$ETa_a_l<-NULL
d2a$ETa_b_l<-NULL
d2a$ETa_c_l<-NULL

#########
#Temperatur recalculation for a,b and c
# dT = 0.9 *h [km] + 3 #Loomis et al. 2012, 2017

dTa_u <- 0.9* 1.78 + 3
dTb_u <- 0.9* 1.16 + 3
dTc_u <- 0.9* 1.45 + 3

dTa_l <- 0.9* 1.78 + 4
dTb_l <- 0.9* 1.16 + 4
dTc_l <- 0.9* 1.45 + 4

for(i in 1:101){
  d2a$ETa_a_u[i] <- ETa(    t_a_ld    = 291.5-dTa_u, 
                          emis_ld   = d2a$em_a[i], 
                          albedo_ld = d2a$alb_a[i],
                          rh        = 0.58,
                          f_ld      = d2a$sma_a[i],
                          cc        = 0.54,
                          ws        = 1.42,
                          a         = 0.39,
                          b         = 0.38,
                          a_2       = 0.22,
                          b_2       = 2.00,
                          cds       = d2a$cd_a[i]*1.4,
                          p         = 81735,
                          r_swc     = 415)
  
  d2a$ETa_b_u[i] <- ETa( t_a_ld    = 295.9-dTb_u, 
                       emis_ld   = d2a$em_b[i],
                       albedo_ld = d2a$alb_b[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_b[i],
                       cc        = 0.61,
                       ws        = 2.13,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = d2a$cd_b[i],
                       p         = 88466.63,
                       r_swc     = 415)
  
  d2a$ETa_c_u[i] <- ETa( t_a_ld    = 293.3-dTc_u, 
                       emis_ld   = d2a$em_c[i],
                       albedo_ld = d2a$alb_c[i],
                       rh        = 0.57,
                       f_ld      = d2a$sma_c[i],
                       cc        = 0.55,
                       ws        = 1.72,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = d2a$cd_c[i]*1.2,
                       p         = 85128,
                       r_swc     = 415)#
  
  
  
}
d3u <- d2a
d3u_r <- d3u
##########################################################################################################
#########################################LGM LOWER
setwd(path_4ndlevel)
d3 <- get(load("20201104_LGM01_TL_LC.RData"))  
d2a<-d3


d2a$pre <- d2a$pre*100

#albedo for each catchment
p_n <- 10
d2a$alb_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$alb_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$alb_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]


#emissivity for each catchment
p_n <- 3
d2a$em_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$em_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$em_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]

#soil moisture availabiltiy for each catchment
p_n <- 11
d2a$sma_a <- d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n]
d2a$sma_b <- d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n]
d2a$sma_c <- d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n]


###Surface drag coefficient for each catchment:

#roughness length:
p_n <- 12
d2a$z0_a <- (d2a$a_2*para[3,p_n]+d2a$a_7*para[8,p_n]+d2a$a_9*para[10,p_n]+d2a$a_10*para[11,p_n]+d2a$a_16*para[17,p_n])/100
d2a$z0_b <- (d2a$b_2*para[3,p_n]+d2a$b_7*para[8,p_n]+d2a$b_9*para[10,p_n]+d2a$b_10*para[11,p_n]+d2a$b_16*para[17,p_n])/100
d2a$z0_c <- (d2a$c_2*para[3,p_n]+d2a$c_7*para[8,p_n]+d2a$c_9*para[10,p_n]+d2a$c_10*para[11,p_n]+d2a$c_16*para[17,p_n])/100

#friction velocity:
d2a$u_a <- (0.4*1.42)/ log(10/d2a$z0_a)
d2a$u_b <- (0.4*2.13)/ log(10/d2a$z0_b)
d2a$u_c <- (0.4*1.72)/ log(10/d2a$z0_c)

#atmospheric resitance:
d2a$r_a_a    <- log(10/d2a$z0_a)/(d2a$u_a * 0.4) #atmospheric resitance
d2a$r_a_b    <- log(10/d2a$z0_b)/(d2a$u_b * 0.4) #atmospheric resitance
d2a$r_a_c    <- log(10/d2a$z0_c)/(d2a$u_c * 0.4) #atmospheric resitance

#surface drag coefficient
d2a$cd_a     <- 0.4*0.4 * (log(d2a$r_a_a /(d2a$z0_a))^-2) #surface drag coefficient
d2a$cd_b     <- 0.4*0.4 * (log(d2a$r_a_b /(d2a$z0_b))^-2) #surface drag coefficient
d2a$cd_c     <- 0.4*0.4 * (log(d2a$r_a_c /(d2a$z0_c))^-2) #surface drag coefficient

d2a$ETa_a_u<-NULL
d2a$ETa_b_u<-NULL
d2a$ETa_c_u<-NULL

d2a$ETa_a_l<-NULL
d2a$ETa_b_l<-NULL
d2a$ETa_c_l<-NULL

#########
#Temperatur recalculation for a,b and c
# dT = 0.9 *h [km] + 3 #Loomis et al. 2012, 2017

dTa_u <- 0.9* 1.78 + 3
dTb_u <- 0.9* 1.16 + 3
dTc_u <- 0.9* 1.45 + 3

dTa_l <- 0.9* 1.78 + 4
dTb_l <- 0.9* 1.16 + 4
dTc_l <- 0.9* 1.45 + 4

for(i in 1:101){
  
  #############loweer
  d2a$ETa_a_l[i] <- ETa(    t_a_ld    = 291.5-dTa_l, 
                            emis_ld   = d2a$em_a[i], 
                            albedo_ld = d2a$alb_a[i],
                            rh        = 0.58,
                            f_ld      = d2a$sma_a[i],
                            cc        = 0.54,
                            ws        = 1.42,
                            a         = 0.39,
                            b         = 0.38,
                            a_2       = 0.22,
                            b_2       = 2.00,
                            cds       = d2a$cd_a[i]*1.4,
                            p         = 81735,
                            r_swc     = 415)
  
  d2a$ETa_b_l[i] <- ETa( t_a_ld    = 295.9-dTb_l, 
                         emis_ld   = d2a$em_b[i],
                         albedo_ld = d2a$alb_b[i],
                         rh        = 0.57,
                         f_ld      = d2a$sma_b[i],
                         cc        = 0.61,
                         ws        = 2.13,
                         a         = 0.39,
                         b         = 0.38,
                         a_2       = 0.22,
                         b_2       = 2.00,
                         cds       = d2a$cd_b[i],
                         p         = 88466.63,
                         r_swc     = 415)
  
  d2a$ETa_c_l[i] <- ETa( t_a_ld    = 293.3-dTc_l, 
                         emis_ld   = d2a$em_c[i],
                         albedo_ld = d2a$alb_c[i],
                         rh        = 0.57,
                         f_ld      = d2a$sma_c[i],
                         cc        = 0.55,
                         ws        = 1.72,
                         a         = 0.39,
                         b         = 0.38,
                         a_2       = 0.22,
                         b_2       = 2.00,
                         cds       = d2a$cd_c[i]*1.2,
                         p         = 85128,
                         r_swc     = 415)#
  
  
}
d3l <- d2a
d3l_r <- d3l













#########################################
d2a_lgm_l <- d2a

d2a<- d2a[1:51,]
d3l<- d3l[1:51,]
d3u<- d3u[1:51,]


par(mfrow=c(2,3))  
ra <- c(50,100)


plot(d2a$pre, d3u$ETa_a_u,type="l",xlim=ra,ylim=c(200,1400),col="black",
     main="Abaya",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_a_MD,col="darkgreen",pch=16)
points(100,1123,col="red",pch=16)

lines(d2a$pre, d3l$ETa_a_l)

abline(h=1123,col="red" )
lines(d2a$pre,d2a$pre/100*1147,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1123),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_a_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_a_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)



plot(d2a$pre, d3u$ETa_c_u,type="l",xlim=ra,ylim=c(200,1400),col="black",
     main="Chamo",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_c_MD,col="darkgreen",pch=16)
points(100,1060,col="red",pch=16)

lines(d2a$pre, d3l$ETa_c_l)


abline(h=1060,col="red" )
lines(d2a$pre,d2a$pre/100*1131,col="blue")

mtext(text = paste("ETa_MD_orig=",round(1060),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_c_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_c_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="black",adj=0,cex=0.8)


plot(d2a$pre, d3u$ETa_b_u,type="l",xlim=ra,ylim=c(200,1400),col="black",
     main="Bahir",xlab="Precipitation compared to Modern Day (%)",ylab="Water Flux (mm/a)")
points(100,ETa_b_MD,col="darkgreen",pch=16)
points(100,892,col="red",pch=16)

lines(d2a$pre, d3l$ETa_b_l)

abline(h=892,col="red" )
lines(d2a$pre,d2a$pre/100*893,col="blue")

mtext(text = paste("ETa_MD_orig=",round(892),"mm"),
      side = 3,#side 1 = bottom
      line = 0,col="red",adj=0,cex=0.8)
mtext(text = paste("ETa_MD=",round(ETa_b_MD),"mm"),
      side = 3,#side 1 = bottom
      line = -2,col="darkgreen",adj=0,cex=0.8)
mtext(text = paste("ETa_PI=",round(ETa_b_PI),"mm"),
      side = 3,#side 1 = bottom
      line = -4,col="darkgreen",adj=0,cex=0.8)

b_c <- c("#267300","#ff9500","#4ecc00","#ffff00")
p_r <- rep(100,4)

plot (d2a$pre, d2a$a_16*100 ,type="l",xlim=ra,ylim=c(0,100),col="#267300")
lines(d2a$pre, d2a$a_2*100+ d2a$a_16*100 ,col="darkgreen")


#polygon(c(d2a$pre, d2a$pre[101:1]), c(d2a$a_2*100, rep(0,101)),
#        border = NA,   col = "#267300")
#polygon(c(d2a$pre, d2a$pre[101:1]), c( (d2a$a_2*100+d2a$a_9*100) , 100*d2a$a_2[101:1]),
#        border = NA,   col = "#4ecc00")
lines(d2a$pre, d2a$a_9*100+d2a$a_2*100+ d2a$a_16*100 ,col="#4ecc00")
lines(d2a$pre, d2a$a_10*100+d2a$a_9*100+d2a$a_2*100+ d2a$a_16*100,col="#ffff00")
lines(d2a$pre, d2a$a_7*100+d2a$a_10*100+d2a$a_9*100+d2a$a_2*100+ d2a$a_16*100 ,col="#ff9500")


#points(p_r ,100*f_MD[c(2,3,5,6),4],pch=16,col=b_c)

plot (d2a$pre, d2a$c_16*100,type="l",xlim=ra,ylim=c(0,100),col="#267300")
lines(d2a$pre, d2a$c_2*100+d2a$c_16*100,col="darkgreen")
lines(d2a$pre, d2a$c_9*100+d2a$c_2*100+d2a$c_16*100,col="#4ecc00")
lines(d2a$pre, d2a$c_10*100+d2a$c_9*100+d2a$c_2*100+d2a$c_16*100,col="#ffff00")
lines(d2a$pre, d2a$c_7*100+ d2a$c_10*100+d2a$c_9*100+d2a$c_2*100+d2a$c_16*100,col="#ff9500")

points(p_r[1:3] ,100*f_MD[c(22,23,24),4],pch=16,col=b_c[c(1,3,4)])



plot (d2a$pre, d2a$b_16*100,type="l",xlim=ra,ylim=c(0,100),col="#267300")
lines(d2a$pre, d2a$b_2*100+d2a$b_16*100,col="black")

lines(d2a$pre, d2a$b_9*100+d2a$b_2*100+d2a$b_16*100,col="#4ecc00")
lines(d2a$pre, d2a$b_10*100+d2a$b_9*100+d2a$b_2*100+d2a$b_16*100,col="#ffff00")
lines(d2a$pre, d2a$b_7*100+d2a$b_10*100+d2a$b_9*100+d2a$b_2*100+d2a$b_16*100,col="#ff9500")

points(p_r ,100*f_MD[c(12,14,16,17),4],pch=16,col=b_c)

##########LAKE EVAPORATION:



dTa_l <- 0.9* 1.18  + 3.5
dTb_l <- 0.9* 0.5 + 3.5
dTc_l <- 0.9* 1.11 + 3.5

a_lake <- ETa(    t_a_ld    = 294.7-dTa_l, 
                          emis_ld   = 0.98, 
                          albedo_ld = 0.06,
                          rh        = 0.57,
                          f_ld      = 1,
                          cc        = 0.31,
                          ws        = 1.46,
                          a         = 0.39,
                          b         = 0.38,
                          a_2       = 0.22,
                          b_2       = 2.00,
                          cds       = 7.3E-4,
                          p         = 87974,
                          r_swc     = 415)

b_lake <- ETa( t_a_ld    = 299.5-dTb_l, 
                       emis_ld   =  0.98,
                       albedo_ld = 0.06,
                       rh        = 0.55,
                       f_ld      = 1,
                       cc        = 0.44,
                       ws        = 2.56,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = 7.3E-4,
                       p         = 95671.49,
                       r_swc     = 415)

c_lake <- ETa( t_a_ld    = 295.1-dTc_l, 
                       emis_ld   = 0.98,
                       albedo_ld = 0.06,
                       rh        = 0.57,
                       f_ld      = 1,
                       cc        = 0.37,
                       ws        = 1.71,
                       a         = 0.39,
                       b         = 0.38,
                       a_2       = 0.22,
                       b_2       = 2.00,
                       cds       = 7.3E-4,
                       p         = 88695,
                       r_swc     = 415)#
a_lake - 1513
b_lake - 1908
c_lake - 1550

d3u <- d3u_r
d3l <- d3l_r
setwd(path_6ndlevel)
d3 <- d3u
d3$ETa_a_l <- d3l$ETa_a_l
d3$ETa_b_l <- d3l$ETa_b_l
d3$ETa_c_l <- d3l$ETa_c_l


save(d3,file= "20201120_LGM_LC_ET.RData")








