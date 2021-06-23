# MLFischer 2020-08
# Fischer et al., 2021
# A Phytolith supported Biosphere-Hydrosphere Predictive Model for southern Ethiopia: 
# insights into paleoenvironmental changes and human landscape preferences since the Last Glacial Maximum
###
# BRT training and validation
###
remove(list=ls()) #clean up#
removeTmpFiles(0)

library(raster)
library(rgdal)  
library(gbm)
library(dismo)
library(pROC)

path_2ndlevel<- "../02_data_2nd_level/"
path_3ndlevel<- "../02_data_3nd_brt/"

setwd(path_2ndlevel)
get(load("data_all_test.RData"))
data_all_test<-data_all

get(load("data_all.RData"))

# transmute to data frame... with empty and filled rows
d   <- as.data.frame(data_all,na.rm=F,xy=T)
d_f <- as.data.frame(data_all,na.rm=F,xy=T)
d_n <- as.data.frame(data_all,na.rm=T,xy=T)

# calculate annual average EVI
d$evi_anu<- rowMeans(d[,4:15])

# lc to factor
d$lc <- as.factor(d$lc)

# check levels and set them manually
levels(d$lc)
levels(d$lc)<- c("2 Evergreen Broadleaf Forest",  "7 Open Shrubland", "9 Savanna","10 Grassland" )
as.vector(levels(d$lc))

d$count <- seq(1,length(d$x),1)

# create subgroups d1 for train, d2 fill, d3 with same groud size 
d1<- na.omit(d) 
d2 <- d[!is.element(d$count,d1$count),]
require(dplyr)
set.seed(420)
d3 <- d1 %>% group_by(lc) %>% sample_n(size = 1000)

# calculate gpm annual
d1$gpm_ann <- d1$gpm_jan+d1$gpm_feb +d1$gpm_mar+d1$gpm_apr+d1$gpm_may+d1$gpm_jun+d1$gpm_jul+d1$gpm_aug+d1$gpm_sep+d1$gpm_oct+d1$gpm_nov+d1$gpm_dec

# create weight factor
lct <- 1/table(d1$lc)
d1$wf <- lct[d1$lc]
sum(d1$wf)

#check some stuff
length(d1[,1])+length(d2[,1])
summary(is.element(d$count,d1$count))

###test area:
d_t   <- as.data.frame(data_all_test,na.rm=T,xy=T)
d_t$lc <- as.factor(d_t$lc)

# check levels and set them manually
levels(d_t$lc)
levels(d_t$lc)<- c("2 Evergreen Broadleaf Forest",  "7 Open Shrubland", "9 Savanna","10 Grassland" )
as.vector(levels(d_t$lc))
#####################################################################################
#####################################################################################
# LAND COVER MODEL
# PALEO VEGETATION MODEL
# PREDICTIVE VEGETATION MODEL



# exemplary range of hyperparameters
bafr <- rep(c(0.5,0.75),each=12)
trco <- rep(c(2,3,4,5),each=3,times=2)
lera <- rep(c(0.05,0.01,0.005),times=8)
ntree <- rep(c(50,100,250,300,400,500,1000,10000),1)

# some dummys
rocresult <- NULL

# set seed for reproduceability and ...GO!
set.seed(420)
gbm_LC<- gbm(lc ~.,
                       data = d1[c(3,16:28)],
                       distribution = "multinomial",
                       interaction.depth=1, 
                       shrinkage = 0.05,
                       bag.fraction = 0.5,
                       n.trees= 1000,
                       verbose=T,
                       cv.folds=5,
                       weights=d1$wf , 
                       train.fraction = 0.75,  
                       n.minobsinnode = 5,
                       class.stratify.cv=T)
  # check some model results
  gbm_LC
  
  gbm.perf(  gbm_LC , plot.it = T, oobag.curve = T, overlay =T, method="cv")
  gbm.perf(  gbm_LC , plot.it = T, oobag.curve = T, overlay =T)
  
  print(gbm_LC)
  summary( gbm_LC)
  
  # use the model to predict and take the most probable class as result
  gbm_pr <-predict.gbm(gbm_LC , d1[c(16:28,34)],type="response",  n.trees=298) #prediction
  gbm_pr_cl <-as.factor( apply(gbm_pr, 1, function(row) which.max(row)))
  levels(gbm_pr_cl)<- c("2 Evergreen Broadleaf Forest",  "7 Open Shrubland", "9 Savanna","10 Grassland" )
  
  #calculate the confusion matrix
  require(caret)
  con_mat <- confusionMatrix(gbm_pr_cl, d1$lc, positive = "1")
  
  #calculate the muliclass ROC
  roc_j<- multiclass.roc(d1$lc, apply(gbm_pr, 1, function(row) which.max(row)))
  rs <- roc_j[['rocs']]
  plot.roc(rs[[1]])
  sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
  print(paste("auc", roc_j$auc[[1]]))

  
  ######################## conMat ROC for other treatments
  d1[,39]<- gbm_pr_cl
  # 1 = Wildlife Reserve
  # 3 = Controlled Hunting Area
  # 6 = Nothing
  cr <- 1
  con_mat <- confusionMatrix(d1[which(d1$pa==cr),39], d1$lc[which(d1$pa==cr)], positive = "1")
  con_mat
  
  roc_1<- multiclass.roc(d1[which(d1$pa==cr),39], as.numeric(d1$lc[which(d1$pa==cr)]))
  print(paste("auc", roc_1$auc[[1]]))
  
  # 2 = National Parc
  cr <- 2
  con_mat <- confusionMatrix(d1[which(d1$pa==cr),39], d1$lc[which(d1$pa==cr)], positive = "1")
  con_mat  
  
  roc_2<- multiclass.roc(d1[which(d1$pa==cr),39], as.numeric(d1$lc[which(d1$pa==cr)]))
  print(paste("auc", roc_2$auc[[1]]))
  
  # 3 = Controlled Hunting Area
  cr <- 3
  con_mat <- confusionMatrix(d1[which(d1$pa==cr),39], d1$lc[which(d1$pa==cr)], positive = "1")
  con_mat    
      
  roc_3<- multiclass.roc(d1[which(d1$pa==cr),39], as.numeric(d1$lc[which(d1$pa==cr)]))
  print(paste("auc", roc_3$auc[[1]]))       
  
  
################################################################################################
################################################################################################
##Test area:

  
  gbm_pr <-predict.gbm(gbm_LC , d_t[c(16:28)],type="response",  n.trees=298) #prediction
  
  gbm_pr_cl <-as.factor( apply(gbm_pr, 1, function(row) which.max(row)))
  levels(gbm_pr_cl)<- c("2 Evergreen Broadleaf Forest",  "7 Open Shrubland", "9 Savanna","10 Grassland" )
  
  #calculate the confusion matrix
  require(caret)
  con_mat <- confusionMatrix(gbm_pr_cl, d_t$lc, positive = "1")
  con_mat 
  
  #calculate the muliclass ROC
  roc_j<- multiclass.roc(d_t$lc, apply(gbm_pr, 1, function(row) which.max(row)))
  rs <- roc_j[['rocs']]
  plot.roc(rs[[1]])
  sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
  print(paste("auc", roc_j$auc[[1]]))
  
  
###  
  

setwd(path_3ndlevel)
save(gbm_LC,file= "20200821_gbm_LC_final.RData")
write.csv2(con_mat$table, "20201130_LC_conmat.csv")

#####################################################################################
#####################################################################################
# EVI ANNUAL MODEEL
set.seed(420)
gbm_EVI <- gbm(evi_anu ~.,
                     data = d1[c(35,16:28)],
                     distribution = "gaussian",
                     interaction.depth=1, # richtigen Parameter wählen
                     shrinkage = 0.05,
                     bag.fraction = 0.75,
                     n.trees= 5000,
                     cv.folds=5,
                     verbose=T,
                     train.fraction = 0.75)

#check some model parameter:
  gbm.perf(gbm_EVI, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)[1]
  summary(gbm_EVI)
  gbm_EVI
  plot.gbm(gbm_EVI,level.plot=F,return.grid=F,i.var=c(9,5,13))
  i <-  gbm.perf(gbm_EVI, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)[1]

#predict
gbm_pr <-predict.gbm(gbm_EVI , d1[c(16:28)],type = "response",  n.trees=i) #prediction 


#check scatterplot and r2
plot(d1$evi_anu, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
cor(d1$evi_anu, gbm_pr)^2


#####
#check results for other treatments

# 1 = Wildlife Reserve
# 2 = National Parc
# 3 = Controlled Hunting Area
# 6 = Nothing
cr <- 1 
#correlation:
cor(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
#visualization:
points(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(1,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
#absolute difference
sum((d1$evi_anu[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
sum(d1$evi_anu[which(d1$pa==cr)]) / length(gbm_pr[which(d1$pa==cr)])
#percentage difference:
100* (mean(d1$evi_anu[which(d1$pa==cr)])- mean(gbm_pr[which(d1$pa==cr)]) )/ mean(d1$evi_anu[which(d1$pa==cr)])


cr <- 2 
cor(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,1,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$evi_anu[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
sum(d1$evi_anu[which(d1$pa==cr)]) / length(gbm_pr[which(d1$pa==cr)])
100* (mean(d1$evi_anu[which(d1$pa==cr)])- mean(gbm_pr[which(d1$pa==cr)]) )/ mean(d1$evi_anu[which(d1$pa==cr)])


cr <- 3 
cor(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,0,1,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$evi_anu[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
sum(d1$evi_anu[which(d1$pa==cr)]) / length(gbm_pr[which(d1$pa==cr)])
100* (mean(d1$evi_anu[which(d1$pa==cr)])- mean(gbm_pr[which(d1$pa==cr)]) )/ mean(d1$evi_anu[which(d1$pa==cr)])


cr <- 6 
cor(d1$evi_anu[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$evi_anu[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
sum(d1$evi_anu[which(d1$pa==cr)]) / length(gbm_pr[which(d1$pa==cr)])


sum((d1$evi_anu- gbm_pr))/ length(gbm_pr )
sum(d1$evi_anu) / length(gbm_pr)

########################
#Whats happening in the test dataset?
#look here:
gbm_pr <-predict.gbm(gbm_EVI , d_t[c(16:28)],type="response",  n.trees=298) #prediction
d_t$evi_anu <- (d_t$evi_jan+d_t$evi_feb+d_t$evi_mar+d_t$evi_apr+d_t$evi_may+d_t$evi_jun+d_t$evi_jul+d_t$evi_aug+d_t$evi_sep+d_t$evi_oct+d_t$evi_nov+d_t$evi_dec)/12

plot(d_t$evi_anu, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
cor(d_t$evi_anu, gbm_pr)^2
100*(mean(d_t$evi_anu)-mean(gbm_pr))/mean(d_t$evi_anu)


###
setwd(path_3ndlevel)
save(gbm_EVI,file= "20200821_gbm_EVI_anu_final.RData")


#####################################################################################
#####################################################################################
#####################################################################################
# TC, NTC, NVC

#################################
# Tree Cover TC:
set.seed(420)
gbm_TC <- gbm(tc ~.,
                     data = d1[c(33,16:28)],
                     distribution = "gaussian",
                     interaction.depth=3, # richtigen Parameter wählen
                     shrinkage = 0.05,
                     bag.fraction = 0.5,
                     n.trees= 5000,
                     cv.folds=5,
                     verbose=T,
                     train.fraction = 0.75)


gbm.perf(  gbm_TC , plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)#, method="cv")
summary( gbm_TC )

i <- gbm.perf(  gbm_TC , plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)[1]
gbm_pr <-predict.gbm(gbm_TC , d1[c(16:28)],type = "response",  n.trees=i) #prediction 
gbm.perf(gbm_TC, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)

TC_pr <- gbm_pr
#scatterplot and r2
plot(d1$tc, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlab="Observed Tree Cover (%)",ylab="Predicted Tree Cover (%)")
cor(d1$tc, gbm_pr)^2


###
#check other treatments:

# 1 = Wildlife Reserve
# 2 = National Parc
# 3 = Controlled Hunting Area
# 6 = Nothing
cr <- 1 
cor(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(1,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$tc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 2 
cor(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,1,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$tc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 3 
cor(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,0,1,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$tc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 6 
cor(d1$tc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$tc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

##check test area:
gbm_pr <-predict.gbm(gbm_TC , d_t[c(16:28)],type="response",  n.trees=298) #prediction
plot(d_t$tc, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlim=c(0,100),ylim=c(0,100),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
cor(d_t$tc, gbm_pr)^2
(mean(d_t$tc)-mean(gbm_pr))

#################################
#################################
#################################
# Non Tree Cover NTC:
set.seed(420)
gbm_NTC <-  gbm(ntc ~.,
                data = d1[c(31,16:28)],
                distribution = "gaussian",
                interaction.depth=5, # richtigen Parameter wählen
                shrinkage = 0.05,
                bag.fraction = 0.5,
                n.trees= 5000,
                cv.folds=5,
                verbose=T,
                train.fraction = 0.75)

par(mfrow=c(1,1))
summary( gbm_NTC )
gbm.perf(gbm_NTC, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)
i <- gbm.perf(  gbm_NTC , plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)[1]
gbm_pr <-predict.gbm(gbm_NTC , d1[c(16:28)],type = "response",  n.trees=i) #prediction 
NTC_pr <- gbm_pr


#scatterplot and r2
plot(x= d1$ntc, y= gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlab="Observed Non Tree Cover (%)",ylab="Predicted Non Tree Cover (%)")
cor(d1$ntc, gbm_pr)^2

# 1 = Wildlife Reserve
# 2 = National Parc
# 3 = Controlled Hunting Area
# 6 = Nothing
cr <- 1 
cor(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(1,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$ntc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 2 
cor(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,1,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$ntc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 3 
cor(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
points(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,0,1,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
sum((d1$ntc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

cr <- 6 
cor(d1$ntc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$ntc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )

##check test area:
gbm_pr <-predict.gbm(gbm_NTC , d_t[c(16:28)],type="response",  n.trees=i) #prediction
plot(d_t$ntc, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlim=c(0,100),ylim=c(0,100),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
cor(d_t$ntc, gbm_pr)^2
(mean(d_t$ntc)-mean(gbm_pr))

#################################
#################################
#################################
# Non Vegetation Cover NVC:
set.seed(420)
gbm_NVC <- gbm(nvc ~.,
                     data = d1[c(32,16:28)],
                     distribution = "gaussian",
                     interaction.depth=3, # richtigen Parameter wählen
                     shrinkage = 0.05,
                     bag.fraction = 0.5,
                     n.trees= 5000,
                     cv.folds=5,
                     verbose=T,
                     train.fraction = 0.75)

summary( gbm_NVC )
gbm_NVC
gbm.perf(gbm_NVC, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)
i <- gbm.perf(  gbm_NVC , plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE)[1]
gbm_pr <-predict.gbm(gbm_NVC , d1[c(16:28)],type = "response",  n.trees=i) #prediction 
NVC_pr <- gbm_pr

#check plot and r2
plot(d1$nvc, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlab="Observed Non Vegetation Cover (%)",ylab="Predicted Non Vegetation Cover (%)")
cor(d1$nvc, gbm_pr)^2

# 1 = Wildlife Reserve
# 2 = National Parc
# 3 = Controlled Hunting Area
# 6 = Nothing
cr <- 1 
cor(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$nvc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
points(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(1,0,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")

cr <- 2 
cor(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$nvc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
points(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,1,0,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")

cr <- 3 
cor(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$nvc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
points(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)],pch=16,cex=0.5, col=rgb(0,0,1,alpha=0.05) ,xlim=c(0,0.6),ylim=c(0,0.6),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")

cr <- 6 
cor(d1$nvc[which(d1$pa==cr)], gbm_pr[which(d1$pa==cr)])^2
sum((d1$nvc[which(d1$pa==cr)]- gbm_pr[which(d1$pa==cr)] ))/ length(gbm_pr[which(d1$pa==cr)] )
########################

gbm_pr <-predict.gbm(gbm_NVC , d_t[c(16:28)],type="response",  n.trees=298) #prediction
plot(d_t$nvc, gbm_pr,pch=16,cex=0.5, col=rgb(0,0,0,alpha=0.05) ,xlim=c(0,100),ylim=c(0,100),xlab="Observed Annual Enhanced Vegetation Index",ylab="Predicted Annual Enhanced Vegetation Index")
cor(d_t$nvc, gbm_pr)^2
(mean(d_t$ntc)-mean(gbm_pr))



##########
contr <- NVC_pr+NTC_pr+TC_pr

mean(contr)
sd(contr)

boxplot(contr)
plot(NVC_pr+NTC_pr+TC_pr,pch=16,cex=0.05)

setwd(path_3ndlevel)
save(gbm_TC,file= "20201203_gbm_TC.RData")
save(gbm_NTC,file= "20201203_gbm_NTC.RData")
save(gbm_NVC,file= "20201203_gbm_NVC.RData")

#####################################################################################
#####################################################################################
