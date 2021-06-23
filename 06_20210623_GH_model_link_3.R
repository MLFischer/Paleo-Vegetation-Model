
###
#preparing
remove(list=ls()) #clean up#
library(tiff) #just load some packages
library(raster)
library(rgdal)
library(rgeos)
library(sp)

par(mfrow=c(1,1))


path_6ndlevel<- "../02_data_6nd_model_link/"
setwd(path_6ndlevel)

freq_a<-read.csv2("lake_altitude_a.csv")#abaya
freq_b<-read.csv2("lake_altitude_b.csv")#bahir
freq_c<-read.csv2("lake_altitude_c.csv")#chamo

get(load("20201120_LGM_LC_ET.RData"))


b <- list()
b$e_ld  <- 0
b$e_lk  <- 1908-222.664
b$s_bas <- 68
b$p_bas <- 916.5
b$size <- 20649.74

a <- list()
a$e_ld  <- 0
a$e_lk  <- 1513-193.8741
a$s_bas <- 260
a$p_bas <- 1407
a$size <-  16203.78
freq_a$lake.area.ration <- freq_a$lake.area.ration/ 100

c <- list()
c$e_ld  <- 0
c$e_lk  <- 1550-200.7402
c$s_bas <- 80
c$p_bas <- 1210.5
c$size <-  1793.164
freq_c$lake.area.ration <- freq_c$lake.area.ration/ 100


#########################################
a$vol_delta <- 0 
b$vol_delta <- 0 
c$vol_delta <- 0 

a$vol_total <- 0
b$vol_total <- 0
c$vol_total <- 0

a$aw_run<- freq_a[1,9]
b$aw_run<- freq_b[1,9]
c$aw_run<- freq_c[1,9]

a$altitude_year <- 0
b$altitude_year <- 0
c$altitude_year <- 0

##################
anu <- 500 #runtime for the model
##################

time <- vector() 
PRE  <-vector()
PRE_2<-vector()

alti_a<-vector()
alti_b<-vector()
alti_c<-vector()

i1<-1

dva<-vector()
dvc<-vector()
dvb<-vector()


b_vol_e_lk <- vector()
b_vol_e_ld <- vector()
b_vol_p    <- vector()


for (k in seq(-50,0,1)){
  
j <- k
k_i <- k+51

a$vol_delta <- 0 
b$vol_delta <- 0 
c$vol_delta <- 0 

a$vol_total <- 0
b$vol_total <- 0
c$vol_total <- 0

a$aw_run<- freq_a[1,9]
b$aw_run<- freq_b[1,9]
c$aw_run<- freq_c[1,9]

a$altitude_year <- 0
b$altitude_year <- 0
c$altitude_year <- 0

p_incr <- 1 + (j/100)  #percentage increased precipitation


comparison <- round(p_incr*100,0) 
d3$pre    <- round(d3$pre,0)
nr <- which(d3$pre == comparison)


a$e_ld  <- d3$ETa_a_u[nr]
b$e_ld  <- d3$ETa_b_u[nr]
c$e_ld  <- d3$ETa_c_u[nr]

a_to_c <- 0
c_to_b <- 0
b_out <- 0

for (i in 1 : anu){
  ################
  #ABAYA DYNAMICS:
  a$vol_p<- (a$p_bas*p_incr) * a$size / 1000000
  a$vol_e_lk<- a$e_lk * a$aw_run[i] * a$size / 1000000
  a$vol_e_ld<- a$e_ld  * (1-a$aw_run[i] ) * a$size / 1000000
  a$vol_delta[i]<-  a$vol_p -  a$vol_e_lk - a$vol_e_ld - (a$s_bas*a$size / 1000000)
  
  if (i == 1){a$vol_total[i] <- a$vol_delta[i]}
  if (i > 1){a$vol_total[i] <- a$vol_total[i-1] + a$vol_delta[i]}
  a$vol_total[a$vol_total<0]<-0
  
  if (a$vol_total[i]<freq_a[1,7]){
    a$altitude_year[i]=freq_a[1,2]
    a$aw_run[i+1]=  freq_a[1,9]
    a_to_c[i] <- 0
    
  }else{
    a$altitude_year[i] = approx(freq_a[,7], freq_a[,2], xout=a$vol_total[i])$y
    a$aw_run[i+1] = approx(freq_a[,7], freq_a[,9], xout=a$vol_total[i])$y
    a_to_c[i] <- 0
  }
  if (a$altitude_year[i] > 1194){
    a$altitude_year[i] = 1194
    a$vol_total[i] <- 25.244876
    a_to_c[i] <- a$vol_delta[i]
  } 
  # a_to_c[i] <- 0 #turn the drainage on or off <3
  ################
  #CHAMO DYNAMICS:
  c$vol_p<- (c$p_bas*p_incr) * c$size / 1000000
  c$vol_e_lk<- c$e_lk * c$aw_run[i] * c$size / 1000000
  c$vol_e_ld<- c$e_ld  * (1-c$aw_run[i] ) * c$size / 1000000
  c$vol_delta[i]<-  c$vol_p -  c$vol_e_lk - c$vol_e_ld + a_to_c[i] - (c$s_bas*c$size / 1000000)
  if (i == 1){c$vol_total[i] <- c$vol_delta[i]}
  if (i > 1){c$vol_total[i] <- c$vol_total[i-1] + c$vol_delta[i]}
  c$vol_total[c$vol_total<0]<-0
  
  if (c$vol_total[i]<freq_c[1,7]){
    c$altitude_year[i]=freq_c[1,2]
    c$aw_run[i+1]=  freq_c[1,9]
    c_to_b[i] <- 0
  }else{
    c$altitude_year[i] = approx(freq_c[,7], freq_c[,2], xout=c$vol_total[i])$y
    c$aw_run[i+1] = approx(freq_c[,7], freq_c[,9], xout=c$vol_total[i])$y
    c_to_b[i] <- 0
  }
  
  if (c$altitude_year[i] > 1123){
    c$altitude_year[i] = 1123
    c$vol_total[i] <- 5.1600304
    c_to_b[i] <- c$vol_delta[i]
  }  
  # c_to_b[i] <- 0 #turn the drainage on or off <3
  
  ################
  #BAHIR DYNAMICS:
  b$vol_p<- (b$p_bas*p_incr) * b$size / 1000000
  b$vol_e_lk<- b$e_lk * b$aw_run[i] * b$size / 1000000
  b$vol_e_ld<- b$e_ld  * (1-b$aw_run[i] ) * b$size / 1000000
  b$vol_delta[i]<-  b$vol_p -  b$vol_e_lk - b$vol_e_ld + c_to_b[i] - (b$s_bas* b$size / 1000000)
  if (i == 1){b$vol_total[i] <- b$vol_delta[i]}
  if (i > 1){b$vol_total[i] <- b$vol_total[i-1] + b$vol_delta[i]} 
  b$vol_total[b$vol_total<0]<-0
  
  if (b$vol_total[i]<freq_b[1,7]){
    b$altitude_year[i]=freq_b[1,2]
    b$aw_run[i+1]=  freq_b[1,9]
    b_out[i]<- 0
  }else{
    b$altitude_year[i] = approx(freq_b[,7], freq_b[,2], xout=b$vol_total[i])$y
    b$aw_run[i+1] = approx(freq_b[,7], freq_b[,9], xout=b$vol_total[i])$y
    b_out[i]<- 0
  }
  
  if (b$altitude_year[i] > 543){
    b$altitude_year[i] = 543
    b$vol_total[i] <- 83.1698643
    b_out[i] <- b$vol_delta[i]
  }  
  
}


if (b$altitude_year[anu]==543){
print(paste(p_incr,"-",b$altitude_year[500],"-",min(which(b$altitude_year==543))))
  time[i1]<-min(which(b$altitude_year==543))
  PRE_2[i1]<-p_incr
  
  b_vol_e_ld[i1] <- b$vol_e_ld
  b_vol_e_lk[i1] <-  b$vol_e_lk
  i1 <- i1 +1
}else{print(paste(p_incr,"-",b$altitude_year[anu]))}




PRE[k_i]<-p_incr
dva[k_i]<-max(a_to_c)
dvb[k_i]<-max(c_to_b)
dvc[k_i]<-max(b_out)
b_vol_p[k_i]    <-   b$vol_p
b$altitude_year<- c(498,b$altitude_year)
a$altitude_year<- c(1176,a$altitude_year)
c$altitude_year<- c(1109,c$altitude_year)

a$altitude_year <- a$altitude_year -1176
c$altitude_year <- c$altitude_year - 1109
b$altitude_year <- b$altitude_year - 498

alti_a[k_i]<-a$altitude_year[anu]
alti_b[k_i]<-b$altitude_year[anu]
alti_c[k_i]<-c$altitude_year[anu]
}

#########################################

par(mfrow = c(1,3))
ra <- c(0.50,1)

plot(PRE,alti_b, type="l",ylim=c(0,50),xlab="Precipiation",ylab="Chew Bahir Lake Level",xlim=ra,pch=16)
lines(PRE,alti_a,col="blue",pch=16)
lines(PRE,alti_c,col="red",pch=16)
##
plot(PRE,dvc,pch=16 ,ylim=c(0,10),xlim=ra,type="l")
lines(PRE,dva,col="blue",pch=16)
lines(PRE,dvb, col="red",pch=16)
mtext("#Bahir= black
#Abya=blue
#Chamo=red")


rel2<- dvb/(b_vol_p+dvb)
plot(PRE,rel2,type="l",xlim=ra,ylim=c(0,0.25))

alti_b[which(alti_b[-1]-alti_b[-51]!=0)][-1]
PRE[which(alti_b[-1]-alti_b[-51]!=0)][-1]


