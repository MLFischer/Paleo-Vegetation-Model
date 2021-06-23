# MLFischer 2020-08
# Fischer et al., 2021
# A Phytolith supported Biosphere-Hydrosphere Predictive Model for southern Ethiopia: 
# insights into paleoenvironmental changes and human landscape preferences since the Last Glacial Maximum
###
# data insights
###
remove(list=ls()) #clean up#
removeTmpFiles(0)

require(ggplot2)
require(raster)
require(dplyr)

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

path_2ndlevel<- "../02_data_2nd_level/"

setwd(path_2ndlevel)
get(load("data_all.RData"))
get(load("data_all_pr.RData"))
####
names(data_all)
d   <- as.data.frame(data_all,na.rm=F,xy=T)

d$evi_anu<- rowMeans(d[,4:15])
d$lc <- as.factor(d$lc)

levels(d$lc)
levels(d$lc)<- c("2 Evergreen Broadleaf Forest",  "7 Open Shrubland", "9 Savanna","10 Grassland" )
as.vector(levels(d$lc))

d$count <- seq(1,length(d$x),1)

#d1 ist trainingsmenge
#d2 ist auffüller und dann noch order
#d3  subsamle wo alle gruppen gleich groß sind
d1<- na.omit(d) 
d2 <- d[!is.element(d$count,d1$count),]


set.seed(42)
d3 <- d1 %>% group_by(lc) %>% sample_n(size = 2000)
d1$gpm_ann <- d1$gpm_jan+d1$gpm_feb +d1$gpm_mar+d1$gpm_apr+d1$gpm_may+d1$gpm_jun+d1$gpm_jul+d1$gpm_aug+d1$gpm_sep+d1$gpm_oct+d1$gpm_nov+d1$gpm_dec

################################
##################################################

#size of the study area:
area(area.shp)/1000000

nr    <- data_all$lc
nr_pr <- data_all_pr$gpm_jan

#pixelnumber of tr and pr
nr[which(is.na(nr[])==F)]<-1
sum(nr[],na.rm=T)

nr_pr[which(is.na(nr_pr[])==F)]<-1
sum(nr_pr[],na.rm=T)

#elevation summary
cellStats(data_all_pr$srtm,stat=mean)
cellStats(data_all_pr$srtm,stat=max)
cellStats(data_all_pr$srtm,stat=min)

cellStats(data_all$srtm,stat=mean)
cellStats(data_all$srtm,stat=max)
cellStats(data_all$srtm,stat=min)

#precipitation summary annual:
data_all$gpm_anu <- (data_all$gpm_jan+data_all$gpm_feb
+data_all$gpm_mar
+data_all$gpm_apr
+data_all$gpm_jun
+data_all$gpm_jul
+data_all$gpm_aug
+data_all$gpm_sep
+data_all$gpm_oct
+data_all$gpm_nov
+data_all$gpm_dec)

data_all_pr$gpm_anu <- (data_all_pr$gpm_jan+data_all_pr$gpm_feb
+data_all_pr$gpm_mar
+data_all_pr$gpm_apr
+data_all_pr$gpm_jun
+data_all_pr$gpm_jul
+data_all_pr$gpm_aug
+data_all_pr$gpm_sep
+data_all_pr$gpm_oct
+data_all_pr$gpm_nov
+data_all_pr$gpm_dec)

cellStats(data_all_pr$gpm_anu,stat=mean)
cellStats(data_all_pr$gpm_anu,stat=max)
cellStats(data_all_pr$gpm_anu,stat=min)

cellStats(data_all$gpm_anu,stat=mean)
cellStats(data_all$gpm_anu,stat=max)
cellStats(data_all$gpm_anu,stat=min)

##LC distribution:
LC_d <- as.data.frame(freq(data_all$lc)[1:4,])
LC_d$perc <- 100* LC_d$count/sum(LC_d[,2])
LC_d


#####################
#EVI annual summary
data_all$evi_anu <- (data_all$evi_jan+data_all$evi_feb
                     +data_all$evi_mar
                     +data_all$evi_apr
                     +data_all$evi_jun
                     +data_all$evi_jul
                     +data_all$evi_aug
                     +data_all$evi_sep
                     +data_all$evi_oct
                     +data_all$evi_nov
                     +data_all$evi_dec)/12

cellStats(data_all$evi_anu ,stat=mean)
cellStats(data_all$evi_anu ,stat=max)
cellStats(data_all$evi_anu ,stat=min)

#TC, NTC NVC:
cellStats(data_all$tc ,stat=mean)
cellStats(data_all$tc ,stat=max)
cellStats(data_all$tc ,stat=min)

cellStats(data_all$ntc ,stat=mean)
cellStats(data_all$ntc ,stat=max)
cellStats(data_all$ntc ,stat=min)

cellStats(data_all$nvc ,stat=mean)
cellStats(data_all$nvc ,stat=max)
cellStats(data_all$nvc ,stat=min)


#TC NTC and NVC for the classes
d1_2<- d1[which(d1$lc=="2 Evergreen Broadleaf Forest"),]
median(d1_2$tc)
median(d1_2$ntc)
median(d1_2$nvc)


d1_9<- d1[which(d1$lc=="9 Savanna"),]
median(d1_9$tc)
median(d1_9$ntc)
median(d1_9$nvc)

d1_10<- d1[which(d1$lc=="10 Grassland"),]
median(d1_10$tc)
median(d1_10$ntc)
median(d1_10$nvc)

d1_7<- d1[which(d1$lc=="7 Open Shrubland"),]
median(d1_7$tc)
median(d1_7$ntc)
median(d1_7$nvc)

# Mittelwert und Mediane von den Klassen plotten (EVI und GPM)
# Mittelwert der TC NTC NVC SRTM

#plot EVI and quantiles in one plot:
setwd("../00_Results")
pdf(file= "02_plots_EVI_GPM.pdf",width=6.6 ,height=2.5,bg="transparent",family="sans")#,paper="a4")

par(mfrow=c(1,2)
    ,mai = c(0.5, 0.25, 0.1, 0.25)
    ,oma = c(0.1, 1, 0.1, 0.1))


d6 <- d1  %>% group_by(lc)  %>%  dplyr::select(contains("evi"))%>% summarise_each( funs(quantile))
d6<- (t(d6))
a<- rep(c("30","00","50","00","30"),4)
b<- c(rep("#267300",5),rep("#ff9500",5),rep("#4ecc00",5),rep("#ffff00",5))
d6<- rbind(d6, paste(b,a,sep=""))

plot(d6[2:13,3],type="n",ylim=c(0,1),col=d6[15,1],xlim=c(1,12),ylab="EVI",xlab="Month of the Year")

for(i in c(3,8,13,18)){
  lines(d6[2:13,i],col="black",pch=16)
}

polygon(c(1:12, 12:1), c(d6[2:13,4], rev(d6[2:13,2])),
         border = NA,   col = d6[15,3])
polygon(c(1:12, 12:1), c(d6[2:13,9], rev(d6[2:13,7])),
        border = NA,   col = d6[15,8])
polygon(c(1:12, 12:1), c(d6[2:13,14], rev(d6[2:13,12])),
        border = NA,   col = d6[15,13])
polygon(c(1:12, 12:1), c(d6[2:13,19], rev(d6[2:13,17])),
        border = NA,   col = d6[15,18])


#plot GPM quantiles in one plot per class:
d6 <- d1  %>% group_by(lc)  %>%  dplyr::select(contains("gpm"))%>% summarise_each( funs(quantile))
d6<- (t(d6))
a<- rep(c("30","00","50","00","30"),4)
b<- c(rep("#267300",5),rep("#ff9500",5),rep("#4ecc00",5),rep("#ffff00",5))
d6<- rbind(d6, paste(b,a,sep=""))

plot(d6[2:13,3],type="n",ylim=c(0,300),col=d6[15,1],xlim=c(1,12),ylab="Precipitation (mm per year)",xlab="Month of the Year")

for(i in c(3,8,13,18)){
  lines(d6[2:13,i],col="black",pch=16)
}


polygon(c(1:12, 12:1), c(d6[2:13,4], rev(d6[2:13,2])),
        border = NA,   col = d6[15,3])
polygon(c(1:12, 12:1), c(d6[2:13,9], rev(d6[2:13,7])),
        border = NA,   col = d6[15,8])
polygon(c(1:12, 12:1), c(d6[2:13,14], rev(d6[2:13,12])),
        border = NA,   col = d6[15,13])
polygon(c(1:12, 12:1), c(d6[2:13,19], rev(d6[2:13,17])),
        border = NA,   col = d6[15,18])
dev.off()


###Boxplots for each class:
pdf(file= "02_plots_VCF_SRTM.pdf",width=6.6 ,height=2.5,bg="transparent",family="sans")#,paper="a4")
par(mfrow=c(1,4)
    ,mai = c(0.5, 0.25, 0.1, 0.25)
    ,oma = c(0.1, 1, 0.1, 0.1))


#plot TC NTC, NVC, SRTM boxplots per group:
b_c <- c("#267300","#ff9500","#4ecc00","#ffff00")

#pdf(file= "02_VCF_SRTM_boxplots.pdf",width=10 ,height=5,bg="transparent",family="sans")#,paper="a4")
#par(mfrow=c(1,4))

boxplot(d1$tc~d1$lc,col=b_c,ylim=c(0,100),ylab="Tree Cover (%)",xlab="Landcover",xaxt="n")
boxplot(d1$ntc~d1$lc,col=b_c,ylim=c(0,100),ylab="Non Tree Cover (%)",xlab="",xaxt="n")
boxplot(d1$nvc~d1$lc,col=b_c,ylim=c(0,100),ylab="Non Vegetation Cover (%)",xlab="",xaxt="n")
boxplot(d1$srtm~d1$lc,col=b_c,ylim=c(0,3500),ylab="Elevation (m a.s.l.)",xlab="",xaxt="n")

dev.off()

setwd(path_2ndlevel)
save(d1,file= "20201211_d1_test_dataset.RData")
###############################

xf <- freq(data_all$srtm, digits=0, useNA="no", merge=T)
xf<-as.data.frame(xf)
plot(xf)
for (i in 1:length(xf[,1])){
  xf[i,3]<- sum(xf[1:i,2])
  
}
plot(xf[,3],xf[,1])
