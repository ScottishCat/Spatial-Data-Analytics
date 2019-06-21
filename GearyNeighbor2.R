library(raster)
library(rgdal)#read the shape file
library(spatstat)
library(sp)
library(maptools)
library(ggplot2)
library(ggmap)
library(maps)
library(GISTools)
library(spdep)

# Geary'C Index: You need to have a function for computing Geary'C Index,
gearyXc()
# This function takes a vector of values and a matrix of weights(W) measures
#spatial autocorrelation using Geary's C.
gearyXc <- function(values, Weights) {
  y <- values
  w <- Weights
  n <- length(y)
  m <- mean(y)
  DL <- 0;
  for(i in 1:n){
    DL <- DL + (y[i]-m)^2;
  }
  L <- (n-1)/DL
  NR <- 0;
  DR <- 0;
  for(i in 1:n){
    for(j in 1:n){
      NR <- NR + w[i,j]*(y[i]-y[j])^2
      DR <- DR + 2*w[i,j]
    }
  }
  C<-L * NR/DR
  return(C)
}
######################################################
#getArea() function : you should define a function for returning area of
#input polygon data
getArea<-function(PolygonData){
  n<-length(PolygonData)
  result<-vector();
  for(i in 1:n){
    v <-nrow(PolygonData@polygons[[i]]@Polygons[[1]]@coords)
    area=0
    for(j in 1:v){
      if(j==v){
        p=1
      }else{
        p=j+1
      }
      x1<-PolygonData@polygons[[i]]@Polygons[[1]]@coords[j,1]
      y1<-PolygonData@polygons[[i]]@Polygons[[1]]@coords[j,2]
      x2<-PolygonData@polygons[[i]]@Polygons[[1]]@coords[p,1]
      y2<-PolygonData@polygons[[i]]@Polygons[[1]]@coords[p,2]
      area<-area+(x1*y2-y1*x2)/2
    }
    result[i]<-abs(area)
  }
  return(result)
}
# Computing Area Vector for a shapefile type of polygon named X:
X <- readOGR("C:/Users/rvais/Wdirectory/Finals", "1_Neighbor2")


X.area <- getArea(X)
#calculate Geary C for polygon X
X.queen <- poly2nb(X,queen=TRUE)
X.rook <- poly2nb(X,queen=FALSE)
X.matrix.queen <- nb2mat(X.queen)
X.matrix.rook <- nb2mat(X.rook)
queen2<-gearyXc(X.area,X.matrix.queen)
queen2
rook2<-gearyXc(X.area,X.matrix.rook)
rook2
#Moran'I and Geary's C using R package
X.nb2listw.queen <- nb2listw(X.queen)
X.nb2listw.rook <- nb2listw(X.rook)
geary.test(X.area,X.nb2listw.queen)
geary.test(X.area,X.nb2listw.rook)

