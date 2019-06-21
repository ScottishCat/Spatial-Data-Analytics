#Part-B
library(readr)
library(ggplot2)
library(gstat)
library(rspatial)
library(sp)
library(nlme)
library(scales)
library(dismo)
library(raster)
library(rgdal)
library(phylin)
library(sp)
library(rgdal)
library(spatstat)
library(maptools)
library(GISTools)
library(plyr)
library(gstat)
library(automap)


#Loading the dataset
community=readOGR("C:/Users/rvais/Wdirectory/Finals", "2_Community")
plot(community)
x<-community$S_Long
y<-community$S_Lat
centroidX<-coordinates(community)[,1]
centroidY<-coordinates(community)[,2]
nx<-length(x)
ncenter<-length(centroidX)
Rslt<-list()
ctr<-data.frame()
for(j in 1:ncenter){
  data<-data.frame()
  for(i in 1:nx){
    distance<-sqrt((x[i]-centroidX[j])^2+(y[i]-centroidY[j])^2)
    distance<-as.numeric(distance)
    data<-rbind(data,data.frame(x=x[i],y=y[i],distance=(distance)^(1),Intensity=community$Intensity[i],centroidX=centroidX[j],centroidY=centroidY[j]))
  }
  names(data)[4]<-"z"
  idx<-order(data$distance)
  data<-data[idx,]
  inverse_distance<-1/data$distance
  inverse_distance
  weight<-(inverse_distance)/sum(inverse_distance)
  weight
  weight_z<-data$z*weight
  weight_z
  IDW.table<-cbind(data,data.frame("Inversedistance"=inverse_distance,weight=weight,"weighted value"=weight_z))
  IDW.table
  cat('total z for each center',community@data$Intensity[j],':')
  z<-sum(IDW.table$weighted.value)
  print(z)
  ctr<-rbind(ctr,data.frame(x=centroidX[j],y=centroidY[j],z=z))
  Rslt[[j]]<-IDW.table
}
Rslt
ctr
View(ctr)
IDW<-as.data.frame(ctr)
coordinates(IDW) = ~x + y
x.range<-as.numeric(bbox(IDW)[1,])
y.range<-as.numeric(bbox(IDW)[2,])
grd <- expand.grid(x = seq(from = x.range[1], to =x.range[2], by = 0.01),
                   y = seq(from = y.range[1],to = y.range[2], by = 0.01)) #
#expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE
dat.idw <- gstat::idw(ctr$z ~ 1, locations=IDW,newdata=grd)
OP<- par( mar=c(0,0,0,0))
image(dat.idw,"var1.pred",col=terrain.colors(20))
contour(dat.idw,"var1.pred", add=TRUE, nlevels=10, col="#656565")
plot(IDW, add=TRUE, pch=16, cex=0.5)
text(coordinates(community), as.character(round(ctr$z,1)), pos=4, cex=0.8,
     col="blue")
par(OP)
title(main="IDW Interpolation Technique")

#OK Interpolation
x<-community$S_Long
y<-community$S_Lat
centroidX<-coordinates(community)[,1]
centroidY<-coordinates(community)[,2]
nx<-length(x)
ncenter<-length(centroidX)
Rslt<-list()
ctr<-data.frame()
for(j in 1:ncenter){
  dataset<-data.frame()
  for(i in 1:nx){
    dataset<-rbind(dataset,data.frame(x=x[i],y=y[i],Intensity=community$Intensity[i]))
     }
  } 
names(dataset)<-c("x","y","z")

community.length<-nrow(dataset)
centers.data<-data.frame(x=centroidX,y=centroidY,z=NA)
dataset<-rbind(dataset,centers.data)
rownames(dataset)<-1:nrow(dataset)
dataset
vgm1=dataset[1:community.length,]
coordinates(vgm1) = ~x + y
vgm1<-autofitVariogram(z~x+y,vgm1,model="Exp")
vgm1
plot(vgm1)
# function
gamma<-function(distance)
{
  result=14+(15-14)*(1-exp(-(abs(distance)/0.11)))
  return(result)
}
#calculate D
points<-nrow(dataset)
dist<-matrix(nrow = points+1,ncol =points+1)
dim(dist)
for(i in 1:points){
  for(j in 1:points){
    distance<-sqrt((dataset[i,1]-dataset[j,1])^2+(dataset[i,2]-
                                                    dataset[j,2])^2)
    distance<-as.numeric(distance)
    dist[i,j]=distance
    dist[j,i]=distance
  }
}
A<-dist
A<-A[1:(community.length+1),1:(community.length+1)]
dim(A)
A
# functions, need discussing
gama.a<-gamma(A)
gama.a[community.length+1,]=c(rep(1,community.length),0)
gama.a[,community.length+1]=c(rep(1,community.length),0)
diag(gama.a)<-0
gama.a
idx<-which(is.na(dataset$z))#get the index whose value is NA
idx
for(k in idx){
  #calcualte b and d
  d<-dist[k,1:community.length]
  d
  # functions, need discussing
  b<-gamma(d)
  b[length(b)+1]=1
  b
  w<-solve(gama.a,b)
  z<-w[1:community.length]*dataset$z[1:community.length]
  z<-sum(z)
  dataset$z[k]<-z
}
dataset
cat('the centers\' info are')
result<-dataset[idx,]
result
#OK map
dat<-as.data.frame(result)
coordinates(dat) = ~x + y
vgm2<-autofitVariogram(z~x+y,dat,model="Exp")
vgm2
plot(vgm2)
nugget2=0.01
sill2=0.02
range2=0.17
m <- vgm(range=range2, "Exp", nugget=nugget2, psill=sill2)
dat.krg <- krige(dat$z~1, dat, grd, model = m)
OP<- par( mar=c(0,0,0,0))
image(dat.krg ,"var1.pred",col=terrain.colors(20))
contour(dat.krg ,"var1.pred", add=TRUE, nlevels=10)
plot(dat, add=TRUE, pch=16, cex=0.5)
text(coordinates(dat), as.character(round(dat$z,1)), pos=4, cex=0.8,
     col="blue")
par(OP)
title(main="Ordinary Kriging Interpolation Technique")

