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
### clean parturition search table
###
d.part$alertdate=str_sub(as.character(d.part$alertdate),1,8)
d.part$searchdate=str_sub(as.character(d.part$searchdate),1,8)
d.part$finishtime=str_sub(as.character(d.part$finishtime),9,17)
d.part$secondalertdate=str_sub(as.character(d.part$secondalertdate),1,8)
d.part$secondstarttime=str_sub(as.character(d.part$secondstarttime),9,17)
d.part$secondfinishtime=str_sub(as.character(d.part$secondfinishtime),9,17)
d.part$secondsearchdate=str_sub(as.character(d.part$secondsearchdate),1,8)

d.part$juliansearch=yday(mdy(d.part$searchdate))

# names(d.part)
# unique(d.part$lowtag)
n.part=length(d.part$lowtag)

#reorder part search data by lowtag/lowtag
d.part=d.part[order(d.part$lowtag),]

#extract the individual ids
individs=d.part$lowtag

###
### Load data GPS location data
###
read.table(paste("~/Documents/Data/GPS_Lotek/PositionData/",d.part$lowtag[1],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(i in 1:n.part){
  d.temp = read.table(paste("~/Documents/Data/GPS_Lotek/PositionData/",d.part$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
  d.temp$lowtag = d.part$lowtag[i]
  names(d.temp)=tolower(gsub('[[:punct:]]',"",names(d.temp)))
  d=rbind(d,as.matrix(d.temp))
}
d=d[-1,]
d=data.frame(d,stringsAsFactors = FALSE)
names(d)

for(j in 1:dim(d)[2]){
  d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
  d[,j]=as.numeric(d[,j])
}

d$lowtag=as.factor(d$lowtag)

head(d)

###
### double checking for duplicate locations
###

d=d[d$latitude!=0,]

summary(duplicated(d))
d[duplicated(d),]

d$datetimegmt=as.numeric(d$datetimegmt)

d$datetimelocal=as.numeric(d$datetimelocal)

#calculating julian day and omitting outside of parturition window
head(as.POSIXct( d$datetimelocal*24*60*60,origin="1899-12-30", tz="CST6CDT"))

d$julian=yday(mdy_hms(d$datetimelocal))


###
### increased fixes parturition window
###

start=yday(mdy("03/20/2017")) # May 6 beginning of parturition period
end=yday(mdy("07/07/2017")) #end of parturition period

###
### subset entire dataset to parturition window
###

d=d[d$julian>start & d$julian <= end,]

###
### increased fixes parturition window
###

start=yday(mdy("03/15/2018")) # May 6 beginning of parturition period
end=yday(mdy("06/20/2017")) #end of parturition period

