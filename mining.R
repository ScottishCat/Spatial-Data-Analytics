#Part B (Industrial Mining)
#loading all the packages required
library(raster)
library(rgdal)#read the shape file
library(spatstat)
library(sp)
library(maptools)

#setting the working directory

setwd("C:/Users/rvais/Wdirectory/SPA2/mining")
dir(getwd())

#Loading the shape file and processing the data
Indmin<- readOGR(dsn= getwd(), layer = "IndustrialMineralMiningOperations2014_10")

#finding the min and max co-ordinates of the shape file(x1min,x2min,x1max,x2max)
x1min=min(Indmin@coords[,1])
x2min=min(Indmin@coords[,2])
x1max=max(Indmin@coords[,1])
x2max=max(Indmin@coords[,2])

#Position plot and this is only for Map plotting.For each of G,F,K and L,I will have to comment this else it doesn't plot properly
par(mai=c(0,0,0.2,0)) 
#plotting map
plot(Indmin, main="Industrial Mining", pch=20)
#adding legend
legend(x1min, x2min, legend="Industrial Mining", col="black", pch=20, cex=0.8)
#Detaching the package GISTools and calling maps right before calling map.scale() because this function map.scale()is used for both maps and GISTools and hence getting masked. 
detach(package:GISTools)
library(maps)
#setting scale using map.scale
map.scale(x=x1max-200000, y=x2min+200, ratio=FALSE, metric=FALSE) 
library(GISTools)
#North arrow
north.arrow(xb=x1max, yb=x2min+100, len=8000, lab="N", col='Blue') 
#adding axes
map.axes(cex.axis=0.8)

#Store coordinates as list and make matrix for the coordinates
coord=list(x=Indmin@coords[,1], y=Indmin@coords[,2])
matx=matrix(unique(c(Indmin@coords[,1],Indmin@coords[,2])), ncol=2)

#Create an object of class "ppp" representing a point pattern dataset in the two-dimensional plane.
Indmin_owin=as.owin(c(x1min,x1max,x2min,x2max))
Indmin_unique=unique(Indmin@coords)
Indmin_ppp=ppp(Indmin_unique[,1], Indmin_unique[,2], Indmin_owin)

#G-function
G_FUNC=Gest(Indmin_ppp, correction="none")
plot(G_FUNC, ylim=c(0,1), main="IndustrialMining G-Function", xlab ="distance d")

#F-function
F_FUNC=Fest(Indmin_ppp, correction="none")
plot(F_FUNC, ylim=c(0,1), main="IndustrialMining F-Function",xlab="distance d")

#K-function
K_FUNC=Kest(Indmin_ppp, correction="none")
plot(K_FUNC, main="OilGasLoc K-Function", xlab ="distance d")

#L-function
L_FUNC=Lest(Indmin_ppp, correction="none")
plot(L_FUNC, main="OilGasLoc L-Function", xlab ="distance d")

