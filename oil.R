#Part-B (Oil-Gas Location PA)
#loading all the packages required
library(raster)
library(rgdal)#read the shape file
library(spatstat)
library(sp)
library(maptools)

#setting the working directory

setwd("C:/Users/rvais/Wdirectory/SPA2/gasoil")
dir(getwd())

#Loading the shape file and processing the data
locPA<- readOGR(dsn= getwd(), layer = "OilGasLocationPA")

#finding the min and max co-ordinates of the shape file(x1min,x2min,x1max,x2max)
x1min=min(locPA@coords[,1])
x2min=min(locPA@coords[,2])
x1max=max(locPA@coords[,1])
x2max=max(locPA@coords[,2])

#Position plot and this is only for Map plotting.For each of G,F,K and L,I will have to comment this else it doesn't plot properly
par(mai=c(0,0,0.2,0)) 
#plotting map
plot(locPA, main="Oil-Gas Location PA", pch=20) 
#adding legend
legend(x1min,x2min, legend="OilGas Loc", col="black", pch=20, cex=0.75)
#Detaching the package GISTools and calling maps right before calling map.scale() because this function map.scale()is used for both maps and GISTools and hence getting masked.
detach(package:GISTools)
library(maps)
#setting scale using map.scale after detach 
map.scale(x=x1max-100000, y=x2min+100, ratio=FALSE, metric=FALSE) 
library(GISTools)
#North arrow
north.arrow(xb=x1max/2, yb=x2min+100, len=8000, lab="N", col='Blue') 
#adding axes
map.axes(cex.axis=0.8)

#Store coordinates as list and make matrix for the coordinates
coord=list(x=locPA@coords[,1], y=locPA@coords[,2])
matx=matrix(unique(c(locPA@coords[,1],locPA@coords[,2])), ncol=2)

#Create an object of class "ppp" representing a point pattern dataset in the two-dimensional plane.
locPA_owin=as.owin(c(x1min,x1max,x2min,x2max))
locPA_unique=unique(locPA@coords)
locPA_ppp=ppp(locPA_unique[,1], locPA_unique[,2], locPA_owin)

#G-function
G_FUNC=Gest(locPA_ppp, correction="none")
plot(G_FUNC, ylim=c(0,1), main="OilGasLoc G-Function", xlab ="distance d")
#F-function
F_FUNC=Fest(locPA_ppp, correction="none")
plot(F_FUNC, ylim=c(0,1), main="OilGasLoc F-Function",xlab="distance d")
#K-function
K_FUNC=Kest(locPA_ppp, correction="none")
plot(K_FUNC,main="OilGasLoc K-Function", xlab ="distance d")
#L-function
L_FUNC=Lest(locPA_ppp, correction="none")
plot(L_FUNC,main="OilGasLoc L-Function", xlab ="distance d")

