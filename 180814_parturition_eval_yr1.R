###
### 8/14/2018 
### Alison Ketz
### Parturition prediction evaluation of anomaly detection algorithm
###


###
### Preliminaries
###

rm(list=ls())

library(geosphere)
library(lubridate)
library(Hmisc)
library(ggplot2)
library(adehabitatLT)
library(mvtnorm)
library(beepr)
library(stringr)
library(dplyr)
library(ggmap)
library(adehabitatHR)
library(maptools)
library(changepoint)
library(sp)
library(spatstat)#for "duplicate" function
library(readr)
library(RColorBrewer)
library(wgutil)

#setwd
setwd("~/Documents/Parturition/180814_parturition_eval_yr1")

###
### Load parturition search data from database
###

d.part=mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "Parturition search")
names(d.part)=tolower(gsub('[[:punct:]]',"",names(d.part)))

names(d.part)[13:16]=c("secondproperty","secondstarttime", "secondfieldcrew","secondfinishtime")
names(d.part)[22]="num_fawnscaptured"
names(d.part)[31]="num_fawnsseennotcaptured" 
names(d.part)[34:36]=c("secondcenterlat","seconndcenterlong","secondtotalnumberofsearchers")

###
###replacing miss-labeled "alertdate" individual, where it says an alert date, although there was no alert
###

d.part[d.part$lowtag==6865,1]=""

###
### Load mortality data
###

d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "Mortalities")
names(d.mort)=tolower(gsub('[[:punct:]]',"",names(d.mort)))
head(d.mort)
morts=d.mort$lowtag

###
### clean parturition search table
###

#formatting dates
d.part$alertdate=str_sub(as.character(d.part$alertdate),1,8)
d.part$searchdate=str_sub(as.character(d.part$searchdate),1,8)
d.part$finishtime=str_sub(as.character(d.part$finishtime),9,17)
d.part$secondalertdate=str_sub(as.character(d.part$secondalertdate),1,8)
d.part$secondstarttime=str_sub(as.character(d.part$secondstarttime),9,17)
d.part$secondfinishtime=str_sub(as.character(d.part$secondfinishtime),9,17)
d.part$secondsearchdate=str_sub(as.character(d.part$secondsearchdate),1,8)

#julian day of search
d.part$juliansearch=yday(mdy(d.part$searchdate))

# names(d.part)
# unique(d.part$lowtag)
n.part=length(d.part$lowtag)

#reorder part search data by lowtag/lowtag
d.part=d.part[order(d.part$lowtag),]

###
### Load data GPS location data
###

#reset wd
setwd("~/Documents/Data/GPS_Lotek/PositionData/")

datalist_temp = list.files(pattern="*.csv")
myfiles = lapply(datalist_temp, read.csv,header=TRUE,stringsAsFactors=FALSE)
nameHeader = tolower(gsub('[[:punct:]]',"",names(myfiles[[1]])))

#reset wd
setwd("~/Documents/Parturition/180814_parturition_eval_yr1")

#clean data
datafiles=lapply(myfiles,function(x){
  names(x)<-nameHeader #fix column names, remove punctuation
  x[,1]<-trimws(x[,1]) #trimwhitespace in deviceid column
  lowtag=rep(substr(x[1,1], 1, 4),dim(x)[1])
  x=data.frame(lowtag,x,stringsAsFactors = FALSE)
  x})

#get list of all ids
# datafiles=lapply(datafiles,function(x){lowtag=rep(substr(x[1,1], 1, 4),dim(x)[1]);x=data.frame(lowtag,x,stringsAsFactors = FALSE);x})
ids=unlist(lapply(datafiles,function(x){x[1,1]}))

#remove morts from overall list of dataframes
datafiles[which(ids %in% morts)]=NULL

#get list of remaining alive ids
ids=unlist(lapply(datafiles,function(x){x[1,1]}))

#remove individuals with messed up GPS collars
badcollars=c(5884,6044,7237,5173,5107,5740,5078,6884,5900,5720,6847,6821)
datafiles[which(ids %in% badcollars)]=NULL

#get list of remaining individuals to run through detector
ids=unlist(lapply(datafiles,function(x){x[1,1]}))
devices=sapply(datafiles,function(x){x[1,3]})

miss.long=sapply(datafiles,function(x){sum(x$longitude==0)})
miss.long

#removing bad GPS fixes where DOP >10
datafiles=lapply(datafiles,function(x){
  x=x[x$dop<11,]
  x
})

x=datafiles[[1]]
miss.latest=sapply(datafiles,function(x){
  x$date_time_local=as.POSIXct(x$datetimegmt*60*60*24, tz = "CST6CDT", origin = "1899-12-30")
  x$julian=yday(x$date_time_local)
  temp=tail(unique(x$julian),2)
  y=sum(x[x$julian==temp[1] |x$julian==temp[2],7]==0)
  y
})

startdate = "2018-03-15 00:00:00 CDT"

datafiles=lapply(datafiles,function(x){
  
  #remove observations with 0 long/lat
  x=x[x$longitude!=0,]
  x=x[x$latitude!=0,]
  
  #impute missing values of lat/long using geosphere midpoint
  if(sum(is.na(x$longitude))==0){x=x}#check if any missing values
  else{
    for(i in 2:(dim(x)[1]-1)){
      if(is.na(x$longitude[i])){
        a=i-1
        while(is.na(x$longitude[a])){a=a-1}
        b=i+1
        while(is.na(x$longitude[b])){b=b+1}
        save = midPoint(cbind(x$longitude[a],x$latitude[a]),cbind(x$longitude[b],x$latitude[b]))
        x$longitude[i] = save[1]
        x$latitude[i] = save[1]
      }
    }
  }#end missing values
  
  # creating date/time columns, converting julian dates
  x$date_time_gmt=as.POSIXct(x$datetimegmt*60*60*24, tz = "GMT", origin = "1899-12-30")
  x$date_time_local=as.POSIXct(x$datetimegmt*60*60*24, tz = "CST6CDT", origin = "1899-12-30")
  
  #Create time lag between successive locations and add as column to all dataframes
  timediff= diff(x[,5])*24
  x=x[-1,]
  x=data.frame(x,timediff,stringsAsFactors = FALSE)
  
  
  #calculate bearings
  x$bear=bearing(cbind(x$longitude,x$latitude))
  
  #calculate distances and step rate
  dist.out = distHaversine(cbind(x$longitude,x$latitude))
  x=x[-1,]
  x$distance = dist.out
  x$step = x$distance/x$timediff
  
  #remove observations prior to start date of spring tracking comparison 
  ### for wtd = March 15, 2018
  x=x[x$date_time_local>startdate,]
  
  #julian day
  x$julian=yday(x$date_time_local)
  
  
  x
})#end function #endlapply

### remove individuals without observations since spring start date
sizes=lapply(datafiles,function(x){
  dim(x)[1]
})
sizes=unlist(sizes)
sizes

removed.indx=which(sizes==0)
removed=ids[removed.indx]
removed 
datafiles[which(sizes==0)]=NULL
ids=ids[-removed.indx]
# length(ids)
# length(datafiles)

### Subset GPS datafiles into those that fawned
d.fawned = datafiles[which(ids %in% d.part[d.part$num_fawnscaptured>0 | !is.na(d.part$numberoffawnscaptured),]$lowtag)]
n.fawned = length(d.fawned)

d.missed = datafiles[which(ids %in% d.part[!(d.part$num_fawnscaptured>0 | !is.na(d.part$numberoffawnscaptured)),]$lowtag)]
n.missed = length(d.missed)



#get list of remaining ids
ids.fawned=unlist(lapply(d.fawned,function(x){x[1,1]}))
ids.missed=unlist(lapply(d.missed,function(x){x[1,1]}))

d.part.fawned=d.part[d.part$num_fawnscaptured>0 | !is.na(d.part$numberoffawnscaptured),]
d.part.missed=d.part[!(d.part$num_fawnscaptured>0 | !is.na(d.part$numberoffawnscaptured)),]
dim(d.part.missed)


for(j in 1:n.fawned){
  d.fawned[[j]]$fawn = 0
  d.fawned[[j]]$fawn=ifelse(d.fawned[[j]]$julian==d.part.fawned$juliansearch[j],1,0)
}


for(j in 1:n.missed){
  d.missed[[j]]$searched = 0
  d.missed[[j]]$searched=ifelse(d.missed[[j]]$julian==d.part.missed$juliansearch[j],1,0)
}




fawnlist=lapply(d.fawned,function(x){x$fawn})
obs.fawn=c()
for(j in 1:n.fawned){
  obs.fawn=c(obs.fawn,sum(fawnlist[[j]]))
}
which(obs.fawn==0) #no data for 6033 after 2 days prior to parturition?!?

#d.fawned[[11]]$julian
#d.part.fawned[d.part.fawned$lowtag==d.fawned[[11]][1,1],]

#fudging for now
d.fawned[[11]]$fawn[dim(d.fawned[[11]])[1]]=1
tail(d.fawned[[11]])

#project data 
#calculate and compile adehabitat features

d.fawned=lapply(d.fawned,function(x){
  
  # setup coordinates
  coords = cbind(x$longitude, x$latitude)
  sp = SpatialPoints(coords)
  
  # make spatial data frame
  # spdf = SpatialPointsDataFrame(coords, x)
  spdf = SpatialPointsDataFrame(sp, x)
  
  # EPSG strings
  latlong = "+init=epsg:4326"
  proj4string(spdf) = CRS(latlong)
  
  x.sp.proj = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000+y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  x=data.frame(x.sp.proj)
  
  #load into adehabitat
  x.traj <- as.ltraj(cbind(x$coords.x1,x$coords.x2), date=x$date_time_local,id=x$lowtag)
  
  #converts traj object to data frame
  x.traj.df = ld(x.traj)
  x.traj.df$id=as.character(x.traj.df$id)
  x.traj.df=rbind(rep(NA,dim(x.traj.df)[2]),x.traj.df)
  x.traj.df=x.traj.df[-dim(x.traj.df)[1],]
  x$relangle=x.traj.df$rel.angle
  x$disttraj=x.traj.df$dist
  x$R2n=x.traj.df$R2n
  x$dx=x.traj.df$dx
  x$dy=x.traj.df$dy
  
  remove.indx=which(is.na(x$disttraj))
  x <- x[-remove.indx,]
  x
})


d.missed=lapply(d.missed,function(x){
  
  # setup coordinates
  coords = cbind(x$longitude, x$latitude)
  sp = SpatialPoints(coords)
  
  # make spatial data frame
  # spdf = SpatialPointsDataFrame(coords, x)
  spdf = SpatialPointsDataFrame(sp, x)
  
  # EPSG strings
  latlong = "+init=epsg:4326"
  proj4string(spdf) = CRS(latlong)
  
  x.sp.proj = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000+y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  x=data.frame(x.sp.proj)
  
  #load into adehabitat
  x.traj <- as.ltraj(cbind(x$coords.x1,x$coords.x2), date=x$date_time_local,id=x$lowtag)
  
  #converts traj object to data frame
  x.traj.df = ld(x.traj)
  x.traj.df$id=as.character(x.traj.df$id)
  x.traj.df=rbind(rep(NA,dim(x.traj.df)[2]),x.traj.df)
  x.traj.df=x.traj.df[-dim(x.traj.df)[1],]
  x$relangle=x.traj.df$rel.angle
  x$disttraj=x.traj.df$dist
  x$R2n=x.traj.df$R2n
  x$dx=x.traj.df$dx
  x$dy=x.traj.df$dy
  
  remove.indx=which(is.na(x$disttraj))
  x <- x[-remove.indx,]
  x
})


#two day moving window summary statistics, centered and scaled, returns list of matrixes of all the features
#first column is id
#second column is julian day


features=lapply(d.fawned,function(x){
  julian.temp = unique(x$julian)
  x.covs = matrix(NA,nr = length(julian.temp),nc = 19)
  x.covs[,1]=as.numeric(x$lowtag[1])
  x.covs[,2]=julian.temp
  x.covs[,19]=0
  x.covs[which(julian.temp==x$julian[x$fawn>0][1]),19] = 1
  for(i in 1:length(julian.temp)){
    x.covs[i,3] = mean(x$step[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.step
    x.covs[i,4] = sd(x$step[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.step
    x.covs[i,5] = mean(x$altitude[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.altitude
    x.covs[i,6] = sd(x$altitude[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.altitude
    x.covs[i,7] = mean(x$tempc[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.temp
    x.covs[i,8] = sd(x$tempc[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.temp
    x.covs[i,9] = mean(x$bear[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.bearing
    x.covs[i,10] = sd(x$bear[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.bearing
    x.covs[i,11] = mean(x$relangle[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.rel.angle
    x.covs[i,12] = sd(x$relangle[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.rel.angle
    x.covs[i,13] = mean(x$R2n[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.R2n
    x.covs[i,14] = sd(x$R2n[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.R2n
    x.covs[i,15] = mean(x$dx[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dx
    x.covs[i,16] = sd(x$dx[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dx
    x.covs[i,17] = mean(x$dy[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dy
    x.covs[i,18] = sd(x$dy[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dy
  }
  x.covs[,3:18]=apply(x.covs[,3:18],2,scale)
  x.covs
})

featureNames=c("mu.step","sig.step","mu.altitude","sig.altitude","mu.temp","sig.temp","mu.bearing","sig.bearing","mu.rel.angle","sig.rel.angle","mu.R2n","sig.R2n","mu.dx","sig.dx","mu.dy","sig.dy")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


features.missed=lapply(d.missed,function(x){
  julian.temp = unique(x$julian)
  x.covs = matrix(NA,nr = length(julian.temp),nc = 19)
  x.covs[,1]=as.numeric(x$lowtag[1])
  x.covs[,2]=julian.temp
  x.covs[,19]=0
  x.covs[which(julian.temp==x$julian[x$searched>0][1]),19] = 1
  for(i in 1:length(julian.temp)){
    x.covs[i,3] = mean(x$step[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.step
    x.covs[i,4] = sd(x$step[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.step
    x.covs[i,5] = mean(x$altitude[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.altitude
    x.covs[i,6] = sd(x$altitude[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.altitude
    x.covs[i,7] = mean(x$tempc[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.temp
    x.covs[i,8] = sd(x$tempc[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.temp
    x.covs[i,9] = mean(x$bear[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.bearing
    x.covs[i,10] = sd(x$bear[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.bearing
    x.covs[i,11] = mean(x$relangle[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.rel.angle
    x.covs[i,12] = sd(x$relangle[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.rel.angle
    x.covs[i,13] = mean(x$R2n[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.R2n
    x.covs[i,14] = sd(x$R2n[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.R2n
    x.covs[i,15] = mean(x$dx[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dx
    x.covs[i,16] = sd(x$dx[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dx
    x.covs[i,17] = mean(x$dy[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dy
    x.covs[i,18] = sd(x$dy[x$julian == julian.temp[i] | x$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dy
  }
  x.covs[,3:18]=apply(x.covs[,3:18],2,scale)
  x.covs
})


###
### printplots
###

##The following is commmented out in order to not remake these plots, when running the code for mapping

###
### for all individuals that fawned with parturition alerts and fawns captured
###

# lapply(features,function(x){
# 
#   out.df = data.frame(x[,2:19])
#   out.df = out.df[-1,]
#   names(out.df)=c("Julian",featureNames,"Fawned")
#   out.df[,18]=as.factor(out.df[,18])
# 
# 
#   p <- list()
#   for(i in 1:4){
#       temp.df=out.df[,c(1,i+1,18)]
#       names(temp.df) = c("Julian","Covariate","Fawned")
#       p[[i]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Fawned),size=2)+
#                 scale_colour_manual(values=cbPalette)+theme_bw()+
#                 geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#                 ylab(featureNames[i])+xlab("Julian")
#   }
# 
#   p2 <- list()
#   for(i in 5:8){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Fawned")
#     p2[[i-4]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Fawned),size=2)+
#                              scale_colour_manual(values=cbPalette)+theme_bw()+
#                              geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#                              ylab(featureNames[i])+xlab("Julian")
#   }
# 
#   p3 <- list()
#   for(i in 9:12){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Fawned")
#     p3[[i-8]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Fawned),size=2)+
#                              scale_colour_manual(values=cbPalette)+theme_bw()+
#                              geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#                              ylab(featureNames[i])+xlab("Julian")
#   }
# 
#   p4 <- list()
#   for(i in 13:16){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Fawned")
#     p4[[i-12]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Fawned),size=2)+
#                              scale_colour_manual(values=cbPalette)+theme_bw()+
#                              geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#                              ylab(featureNames[i])+xlab("Julian")
#   }
# 
#   pdf(paste("190814_featureplots/Features_",x[1,1],".pdf",sep=""))
#   suppressWarnings(multiplot(plotlist=p,cols=1))
#   multiplot(plotlist=p2,cols=1)
#   multiplot(plotlist=p3,cols=1)
#   multiplot(plotlist=p4,cols=1)
#   dev.off()
# })


###
### printing plots of missed/falsepostives
###
# 
# 
# lapply(features.missed,function(x){
#   
#   out.df = data.frame(x[,2:19])
#   out.df = out.df[-1,]
#   names(out.df)=c("Julian",featureNames,"Searched")
#   out.df[,18]=as.factor(out.df[,18])
#   
#   
#   p <- list()
#   for(i in 1:4){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Searched")
#     p[[i]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Searched),size=2)+
#       scale_colour_manual(values=cbPalette)+theme_bw()+
#       geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#       ylab(featureNames[i])+xlab("Julian")
#   }
#   
#   p2 <- list()
#   for(i in 5:8){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Searched")
#     p2[[i-4]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Searched),size=2)+
#       scale_colour_manual(values=cbPalette)+theme_bw()+
#       geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#       ylab(featureNames[i])+xlab("Julian")
#   }
#   
#   p3 <- list()
#   for(i in 9:12){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Searched")
#     p3[[i-8]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Searched),size=2)+
#       scale_colour_manual(values=cbPalette)+theme_bw()+
#       geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#       ylab(featureNames[i])+xlab("Julian")
#   }
#   
#   p4 <- list()
#   for(i in 13:16){
#     temp.df=out.df[,c(1,i+1,18)]
#     names(temp.df) = c("Julian","Covariate","Searched")
#     p4[[i-12]] <- ggplot(temp.df,aes(x=Julian))+geom_point(aes(y=Covariate,col=Searched),size=2)+
#       scale_colour_manual(values=cbPalette)+theme_bw()+
#       geom_hline(yintercept=0,linetype="dashed")+ggtitle(paste(x[1,1],featureNames[i]))+
#       ylab(featureNames[i])+xlab("Julian")
#   }
#   
#   pdf(paste("190814_featureplots/Features_missed_",x[1,1],".pdf",sep=""))
#   suppressWarnings(multiplot(plotlist=p,cols=1))
#   multiplot(plotlist=p2,cols=1)
#   multiplot(plotlist=p3,cols=1)
#   multiplot(plotlist=p4,cols=1)
#   dev.off()
# })

  
###
### Maps
###

d.all = do.call("rbind", d.fawned) 
d.all.missed = do.call("rbind", d.missed) 

###
### Overall map of all individuals that fawns were captured

# Make a bounding box for data
box <- make_bbox(longitude,latitude,data = d.all)
calc_zoom(box)
map=get_map(location=c(mean(d.all$longitude),mean(d.all$latitude)),source="google",zoom="auto",maptype="satellite",crop=box)
#Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=43.069755,-90.280971&zoom=14&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
newPalette=getPalette(18)

pdf("Map_overall_Parturition_2018.pdf")
ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=lowtag,alpha=.5),data = d.all, size = 2)+
  scale_colour_manual(values = newPalette,"ID")+xlab("Longitude")+ylab("Latitude")
dev.off()

###
### Overall map of all individuals that fawns were captured

# Make a bounding box for data
box <- make_bbox(longitude,latitude,data = d.all.missed)
calc_zoom(box)
map=get_map(location=c(mean(d.all.missed$longitude),mean(d.all.missed$latitude)),source="google",zoom="auto",maptype="satellite",crop=box)
#Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=43.069755,-90.280971&zoom=14&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
newPalette=getPalette(18)

pdf("Map_overall_missed_2018.pdf")
ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=lowtag,alpha=.5),data = d.all.missed, size = 2)+
  scale_colour_manual(values = newPalette,"ID")+xlab("Longitude")+ylab("Latitude")+ggtitle("Missed")
dev.off()

###
### mapping points for each individual, good, with fawns captured
###


lapply(d.fawned,function(x){

  d.part.sub=d.part.fawned[d.part.fawned$lowtag==as.integer(x[1,1]),]
  reflong=mean(x$longitude)
  reflat=mean(x$latitude)
  n.leg=dim(x)[1]
  long=x$longitude
  longend=c(x$longitude[2:n.leg],x$longitude[n.leg])
  lat=x$latitude
  latend=c(x$latitude[2:n.leg],x$latitude[n.leg])
  df.leg=data.frame(long,longend,lat,latend)
  box <- make_bbox(longitude,latitude,data = x)
  map=get_map(location=c(reflong,reflat),source="google",zoom=15,maptype="satellite",crop=box)
  # x$fawn=as.factor(x$fawn)
  
  x=rbind(x,rep(NA,dim(x)[2]))
  x$latitude[dim(x)[1]]=d.part.sub$centerlat
  x$longitude[dim(x)[1]]=d.part.sub$centerlong
  x$combo=x$fawn
  x$combo[dim(x)[1]]=2
  x$order=1:dim(x)[1]
  
  
  pdf(file=paste("Maps_Allspring/hits/",x[1,1],"_map_allspring.pdf",sep=""))
    print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=as.factor(combo)),data = x, size = 2,alpha=.5)+
          scale_color_manual(values=cbPalette[1:3],
                             name="Observation",
                             labels=c("Adult - Spring Fixes", "Adult - Day of Capture", "Fawn Captured"))+
            geom_point(aes(longitude,latitude),data=x[x$fawn==1,],size=3,alpha=1,colour=cbPalette[2])+
            geom_path(aes(x = longitude, y = latitude),data = x[x$fawn==1,], size = 1,alpha=1,color=cbPalette[2])+
            geom_text(aes(x=longitude,y=latitude,label=order),data=x[x$fawn==1,],size=1)+
            geom_point(aes(x=birthingbedlong,y=birthingbedlat),data=d.part.sub,size=1,color=cbPalette[3])+
            geom_point(aes(x=centerlat,y=centerlong),data=d.part.sub,size=1,color=cbPalette[3])+
            ggtitle(paste(x[1,1]),"Hit"))
  dev.off()
  
})

###
### for bad individuals, who 
###

lapply(d.missed,function(x){

  reflong=mean(x$longitude)
  reflat=mean(x$latitude)
  n.leg=dim(x)[1]
  long=x$longitude
  longend=c(x$longitude[2:n.leg],x$longitude[n.leg])
  lat=x$latitude
  latend=c(x$latitude[2:n.leg],x$latitude[n.leg])
  df.leg=data.frame(long,longend,lat,latend)
  box <- make_bbox(longitude,latitude,data = x)
  map=get_map(location=c(reflong,reflat),source="google",zoom=15,maptype="satellite",crop=box)
  # x$fawn=as.factor(x$fawn)

  x$order=1:dim(x)[1]
  x$searched=as.integer(x$searched)

  pdf(file=paste("Maps_Allspring/misses/",x[1,1],"_map_allspring.pdf",sep=""))
  print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=as.factor(searched)),data = x, size = 2,alpha=.5)+
          scale_color_manual(values=cbPalette[c(1,7)],
                             name="Observation",
                             labels=c("Adult - Spring Fixes", "Day of Search"))+
          geom_point(aes(longitude,latitude),data=x[x$searched==1,],size=3,alpha=1,colour=cbPalette[7])+
          geom_path(aes(x = longitude, y = latitude),data = x[x$searched==1,], size = 1,alpha=1,color=cbPalette[7])+
          geom_text(aes(x=longitude,y=latitude,label=order),data=x[x$searched==1,],size=1)+
          ggtitle(paste(x[1,1]),"Missed"))
   dev.off()
  
})

#################################################################################################################
###
###  2 day maps
###
#################################################################################################################

lapply(d.fawned,function(x){
  
  x=x[x$julian %in% (x[x$fawn==1,]$julian[1]-0:2),]
  d.part.sub=d.part.fawned[d.part.fawned$lowtag==as.integer(x[1,1]),]
  
  reflong=mean(x$longitude)
  reflat=mean(x$latitude)
  n.leg=dim(x)[1]
  long=x$longitude
  longend=c(x$longitude[2:n.leg],x$longitude[n.leg])
  lat=x$latitude
  latend=c(x$latitude[2:n.leg],x$latitude[n.leg])
  df.leg=data.frame(long,longend,lat,latend)
  box <- make_bbox(longitude,latitude,data = x)
  map=get_map(location=c(reflong,reflat),source="google",zoom=16,maptype="satellite",crop=box)
  
  x=rbind(x,rep(NA,dim(x)[2]))
  x$latitude[dim(x)[1]]=d.part.sub$centerlat
  x$longitude[dim(x)[1]]=d.part.sub$centerlong
  x$combo=x$fawn
  x$combo[dim(x)[1]]=2
  x$combo
  x$order = 1:(dim(x)[1])
  
  pdf(file=paste("Maps_2day/hits/",x[1,1],"_map_2day.pdf",sep=""))
  print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=as.factor(combo)),data = x, size = 2,alpha=.5)+
          scale_color_manual(values=cbPalette[1:3],
                             name="Observation",
                             labels=c("Adult - Spring Fixes", "Adult - Day of Capture", "Fawn Captured"))+
          geom_path(aes(x = longitude, y = latitude),data = x, size = 1,alpha=1,color=cbPalette[1])+
          geom_point(aes(longitude,latitude),data=x[x$fawn==1,],size=3,alpha=1,colour=cbPalette[2])+
          geom_path(aes(x = longitude, y = latitude),data = x[x$fawn==1,], size = 1,alpha=1,color=cbPalette[2])+
          geom_text(aes(x=longitude,y=latitude,label=order),data=x,size=1)+
          geom_point(aes(x=birthingbedlong,y=birthingbedlat),data=d.part.sub,size=1,color=cbPalette[3])+
          geom_point(aes(x=centerlat,y=centerlong),data=d.part.sub,size=1,color=cbPalette[3])+
          ggtitle(x[1,1]))
  dev.off()
  
})


###
### maps 2 day missed
###

lapply(d.missed,function(x){
  x=x[x$julian %in% (x[x$searched==1,]$julian[1]-0:2),]
  d.part.sub=d.part.missed[d.part.missed$lowtag==as.integer(x[1,1]),]
  
  reflong=mean(x$longitude)
  reflat=mean(x$latitude)
  n.leg=dim(x)[1]
  long=x$longitude
  longend=c(x$longitude[2:n.leg],x$longitude[n.leg])
  lat=x$latitude
  latend=c(x$latitude[2:n.leg],x$latitude[n.leg])
  df.leg=data.frame(long,longend,lat,latend)
  box <- make_bbox(longitude,latitude,data = x)
  map=get_map(location=c(reflong,reflat),source="google",zoom=16,maptype="satellite",crop=box)

  x=rbind(x,rep(NA,dim(x)[2]))
  x$latitude[dim(x)[1]]=d.part.sub$birthingbedlat
  x$longitude[dim(x)[1]]=d.part.sub$birthingbedlong
  x$combo=x$searched
  x$combo[dim(x)[1]]=2
  x$order = 1:(dim(x)[1])
  
  pdf(file=paste("Maps_2day/missed/",x[1,1],"_map_2day.pdf",sep=""))
  print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=as.factor(combo)),data = x, size = 2,alpha=.5)+
    scale_color_manual(values=cbPalette[c(1,7,3)],
                       name="Observation",
                       labels=c("Adult - Spring Fixes", "Adult - Day of Search", "Birthing bed if there"))+
    geom_path(aes(x = longitude, y = latitude),data = x, size = 1,alpha=1,color=cbPalette[1])+
    geom_point(aes(longitude,latitude),data=x[x$searched==1,],size=3,alpha=1,colour=cbPalette[7])+
    geom_path(aes(x = longitude, y = latitude),data = x[x$searched==1,], size = 1,alpha=1,color=cbPalette[7])+
    geom_text(aes(x=longitude,y=latitude,label=order),data=x,size=1)+
    geom_point(aes(x=birthingbedlong,y=birthingbedlat),data=d.part.sub,size=1,color=cbPalette[3])+
    ggtitle(paste(x[1,1]),"Missed"))
    dev.off()
  
})


############################################################################################################
###
### Digging into the part search dataframe
###
############################################################################################################

library(plotrix)
library(gridExtra)

###
### Are number of searchers related to number of fawn found? (answer:no)
###

out.lm=lm(num_fawnscaptured~totalnumberofsearchers,data=d.part)
out=summary(out.lm)$coefficients

pdf("fawns_vs_searchers.pdf")
plot(d.part$totalnumberofsearchers,d.part$num_fawnscaptured,main="# fawns vs # searchers",xlab="Total # searchers",ylab="# fawns captured")
abline(coef(out.lm))
addtable2plot(15,2,round(out,3))
dev.off()


###
### Is time searched related to number of fawn found? (answer:no)
###

as.POSIXct(str_sub(as.character(d.part$starttime),10,19),format="%H:%M:%S")

d.part$searchtime=as.POSIXct(as.character(d.part$finishtime),format="%H:%M:%S")-as.POSIXct(str_sub(as.character(d.part$starttime),10,19),format="%H:%M:%S")
searchtime.lm=lm(num_fawnscaptured~searchtime,data=d.part)
searchtime.out=summary(searchtime.lm)$coefficients

pdf("fawns_vs_searchtime.pdf")
plot(d.part$searchtime,d.part$num_fawnscaptured,main="# fawns vs search time (min)",xlab="Search time (min)",ylab="# fawns captured")
abline(coef(searchtime.lm))
addtable2plot(40,3,round(searchtime.out,3))
dev.off()

###
### dates of searches
###

pdf("searches_julian.pdf")
hist(d.part$juliansearch,breaks=8)
dev.off()

###
### properties of searches
###

properties = data.frame(table(d.part$property))
names(properties)[1]="Property"
pdf("propterty_freq.pdf",height=11,width=8.5)
grid.table(properties)
dev.off()


###
### table of birthing beds found
###

names(d.part)
pdf("birthbed_numfawn.pdf",height=11,width=8.5)
grid.table(d.part[order(d.part$juliansearch),c(2,3,8,22)])
dev.off()



###
### false negatives
###

#names(d.part)

fn.notes=d.part$notes[d.part$alertdate==""]
fn.notes=sapply(lapply(fn.notes,strwrap,width=50),paste,collapse="\n")
fn.dates=d.part$searchdate[d.part$alertdate==""]
fn.lowtag=d.part$lowtag[d.part$alertdate==""]
fn.crew=d.part$fieldcrew[d.part$alertdate==""]

pdf("falseneg_notes.pdf",height=8.5,width=11)
grid.table(data.frame(fn.lowtag,fn.notes,fn.dates,fn.crew))
dev.off()

pdf("falseneg_table.pdf",height=8.5,width=11)
grid.table(d.part[d.part$alertdate=="",c(2,3,7,8,21,22)])
dev.off()

###
### Breakdown by searches from alerts
###

dim(d.part[d.part$alertdate!="" & d.part$num_fawnscaptured>0,])

pos.notes=d.part$notes[d.part$alertdate!="" & d.part$num_fawnscaptured>0]
pos.notes=sapply(lapply(pos.notes,strwrap,width=50),paste,collapse="\n")
pos.dates=d.part$searchdate[d.part$alertdate!="" & d.part$num_fawnscaptured>0]
pos.lowtag=d.part$lowtag[d.part$alertdate!="" & d.part$num_fawnscaptured>0]


pdf("pos_notes_page1.pdf",height=11,width=8.5)
grid.table(data.frame(pos.lowtag,pos.notes,pos.dates)[1:8,])
dev.off()
pdf("pos_notes_page2.pdf",height=11,width=8.5)
grid.table(data.frame(pos.lowtag,pos.notes,pos.dates)[9:13,])
dev.off()



#negatives



neg.notes=d.part$notes[d.part$alertdate!="" & d.part$num_fawnscaptured==0]
neg.notes=sapply(lapply(neg.notes,strwrap,width=50),paste,collapse="\n")
neg.dates=d.part$searchdate[d.part$alertdate!="" & d.part$num_fawnscaptured==0]
neg.lowtag=d.part$lowtag[d.part$alertdate!="" & d.part$num_fawnscaptured==0]


pdf("neg_notes_page1.pdf",height=11,width=8.5)
grid.table(data.frame(neg.lowtag,neg.notes,neg.dates)[1:8,])
dev.off()
pdf("neg_notes_page2.pdf",height=11,width=8.5)
grid.table(data.frame(neg.lowtag,neg.notes,neg.dates)[9:13,])
dev.off()




