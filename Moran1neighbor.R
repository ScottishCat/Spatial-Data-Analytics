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

moranXi<-function (values, Weights) {
  y <- values
  w <- Weights
  n <- length(y)
  m <- mean(y)
  sum <- rowSums(w)
  sum[sum == 0] <- 1
  w <- w/sum
  DL<-0;
  for(i in 1:n){
    DL<-DL+(y[i]-m)^2
  }
  I<-n/DL
  NR<-0
  DR<-0
  for(i in 1:n){
    for(j in 1:n){
      NR=NR+w[i,j]*(y[i]-m)*(y[j]-m)
      DR<-DR+w[i,j]
    }
  }
  I<-I*NR/DR
  I<-as.numeric(I)
  E_I_<- -1/(n-1)
  dy <- y-m
  S0 <- sum(w);
  S1 <- 0.5 * sum((w + t(w))^2)
  S2 <- sum((apply(w, 1, sum) + apply(w, 2, sum))^2)
  D <- (sum(dy^4))/(sum(dy^2))^2
  D <- D * n
  A <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
  B <- D * ((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)
  C <- (n - 1) * (n - 2) * (n - 3) * S0^2
  E_I2_<- (A-B)/C
  E_I_2<- E_I_^2
  V_I_<- E_I2_-E_I_2
  z <- (I-E_I_)/sqrt(V_I_)
  p_value <- 1-pnorm(z)
  data.frame(moran=I,p_value=p_value)
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
X <- readOGR("C:/Users/rvais/Wdirectory/Finals", "1_Neighbor1")
#calculate AlleghenyCounty_Municipal's Moran I
X.area <- getArea(X)
#calculateMoran'I and Geary C for polygon X
X.queen <- poly2nb(X,queen=TRUE)
X.rook <- poly2nb(X,queen=FALSE)
X.matrix.queen <- nb2mat(X.queen)
X.matrix.rook <- nb2mat(X.rook)
mi.q <- moranXi(X.area,X.matrix.queen)
mi.q
mi.r <-moranXi(X.area,X.matrix.rook)
mi.r
#Moran'I and Geary's C using R package
X.nb2listw.queen <- nb2listw(X.queen)
X.nb2listw.rook <- nb2listw(X.rook)
moran.test(X.area,X.nb2listw.queen)
moran.test(X.area,X.nb2listw.rook)
