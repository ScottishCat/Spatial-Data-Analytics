#PART A

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

setwd("C:/Users/rvais/Wdirectory/Finals")
getwd()


Neighbor1=readOGR("C:/Users/rvais/Wdirectory/Finals","1_Neighbor1")
Neighbor2=readOGR("C:/Users/rvais/Wdirectory/Finals", "1_Neighbor2")

#Representing shape file as Rook adjacency matrix
rooks = poly2nb(Neighbor1, row.names=Neighbor1$COUNTY, queen=FALSE)
class(rooks)
summary(rooks)

#Plotting the link between the polygons
library(GISTools)
par(mai=c(0,0,0.2,0))
plot(Neighbor1, col='gray', border='blue')
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
legend('bottomleft', legend = 'Neighbor1 Lines', lty =1, lwd = 1, bty = 'o')
xy = coordinates(Neighbor1)
plot(rooks, xy, col='black', lwd=2, add=TRUE)
map.axes(cex.axis=0.8)
title('Neighbor1_rooks')


#Calculating Moran's I using rooks adjacency
n=length(Neighbor1@data$STATE)
wm = nb2mat(rooks, style='B')#A spatial weights matrix reflects the intensity of the geographic relationship between observations
y = Neighbor1$POP_ARR02
ybar = mean(y)
dy = y - ybar
yi = rep(dy, each=n)
yj = rep(dy)
yiyj = yi * yj
pm = matrix(yiyj, ncol=n)#matrix of the multiplied pairs
pmw = pm * wm
spmw = sum(pmw)
smw = sum(wm)
sw  = spmw / smw
vr = n / sum(dy^2)
MI = vr * sw
MI
EI = -1/(n-1) ##Expected value of Moran's I
EI
#Create a 'listw' type spatial weights object to perform significance test
ww =  nb2listw(rooks, style='B')
moran(Neighbor1$POP_ARR02, ww, n=length(Neighbor1@data$STATE), S0=Szero(ww))
moran.test(Neighbor1$POP_ARR02, ww, randomisation=FALSE)
moran.mc(Neighbor1$POP_ARR02, ww, nsim=99)
rwm = mat2listw(wm, style='W')
mat = listw2mat(rwm)
apply(mat, 1, sum)[1:15]
moran.plot(y, rwm)

#queens for Neighbor1
queens = poly2nb(Neighbor1, row.names=Neighbor1$COUNTY, queen=TRUE)
class(queens)
summary(queens)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Neighbor1, col='gray', border='blue')
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Neighbor1_queens')
legend('bottomleft', legend = 'Neighbor1 Lines', lty = 1, lwd = 2, bty = 'o')
xy = coordinates(Neighbor1)
plot(queens, xy, col='black', lwd=2, add=TRUE)
#Calculating Moran's I using queen's adjacency
n=length(Neighbor1@data$STATE)
wm = nb2mat(queens, style='B')
y = Neighbor1$POP_ARR02
ybar = mean(y)
dy = y - ybar
yi = rep(dy, each=n)
yj = rep(dy)
yiyj = yi * yj
pm = matrix(yiyj, ncol=n)
pmw = pm * wm
spmw = sum(pmw)
smw = sum(wm)
sw  = spmw / smw
vr = n / sum(dy^2)
MI = vr * sw
MI

EI = -1/(n-1) ##Expected value of Moran's I
EI
#Create a 'listw' type spatial weights object to perform significance test
ww =  nb2listw(queens, style='B')
ww
moran(Neighbor1$POP_ARR02, ww, n=length(Neighbor1@data$STATE), S0=Szero(ww))
moran.test(Neighbor1$POP_ARR02, ww, randomisation=FALSE)
moran.mc(Neighbor1$POP_ARR02, ww, nsim=99)
rwm = mat2listw(wm, style='W')
mat = listw2mat(rwm)
apply(mat, 1, sum)[1:15]
moran.plot(y, rwm)

#For Neighbor 2
#Representing shape file as Rook adjacency matrix
rooks_neighbor2 = poly2nb(Neighbor2, row.names=Neighbor2$COUNTY, queen=FALSE)
class(rooks_neighbor2)
summary(rooks_neighbor2)

#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
graphics::plot(Neighbor2, col='gray', border='blue')
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Neighbor2_rooks ')
legend('bottomleft', legend = 'Neighbor2 Lines', lty = 1, lwd = 2, bty = 'o')
xy = coordinates(Neighbor2)
plot(rooks_neighbor2, xy, col='black', lwd=2, add=TRUE)

#Calculating Moran's I using rooks adjacency
n=length(Neighbor2@data$STATE)
wm = nb2mat(rooks_neighbor2, style='B')
n=length(Neighbor2@data$STATE)
y = Neighbor2$POP_ARR02
ybar = mean(y)
dy = y - ybar
yi = rep(dy, each=n)
yj = rep(dy)
yiyj = yi * yj
pm = matrix(yiyj, ncol=n)
pmw = pm * wm
spmw = sum(pmw)
smw = sum(wm)
sw  = spmw / smw
vr = n / sum(dy^2)
MI = vr * sw
MI
EI = -1/(n-1) ##Expected value of Moran's I
EI
#Create a 'listw' type spatial weights object to perform significance test
ww =  nb2listw(rooks_neighbor2, style='B')
moran(Neighbor2$POP_ARR02, ww, n=length(Neighbor2@data$STATE), S0=Szero(ww))
moran.test(Neighbor2$POP_ARR02, ww, randomisation=FALSE)
moran.mc(Neighbor2$POP_ARR02, ww, nsim=99)
rwm = mat2listw(wm, style='W')
mat = listw2mat(rwm)
apply(mat, 1, sum)[1:15]
moran.plot(y, rwm)

#queens for Neighbor2
queens_neighbor2 = poly2nb(Neighbor2, row.names=Neighbor2$COUNTY, queen=TRUE)
class(queens_neighbor2)
summary(queens_neighbor2)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Neighbor2, col='gray', border='blue')
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Neighbor2_queens')
legend('bottomleft', legend = 'Neighbor2 Lines', lty = 1, lwd = 2, bty = 'o')
xy = coordinates(Neighbor2)
plot(queens_neighbor2, xy, col='black', lwd=2, add=TRUE)
#Calculating Moran's I using queen's adjacency
n=length(Neighbor2@data$STATE)
wm = nb2mat(queens_neighbor2, style='B')
n=length(Neighbor2@data$STATE)
y = Neighbor2$POP_ARR02
ybar = mean(y)
dy = y - ybar
yi = rep(dy, each=n)
yj = rep(dy)
yiyj = yi * yj
pm = matrix(yiyj, ncol=n)
pmw = pm * wm
spmw = sum(pmw)
smw = sum(wm)
sw  = spmw / smw
vr = n / sum(dy^2)
MI = vr * sw
MI
EI = -1/(n-1) ##Expected value of Moran's I
EI
#Create a 'listw' type spatial weights object to perform significance test
ww =  nb2listw(queens_neighbor2, style='B')
ww
moran(Neighbor2$POP_ARR02, ww, n=length(Neighbor2@data$STATE), S0=Szero(ww))
moran.test(Neighbor2$POP_ARR02, ww, randomisation=FALSE)
moran.mc(Neighbor2$POP_ARR02, ww, nsim=99)
rwm = mat2listw(wm, style='W')
mat = listw2mat(rwm)
apply(mat, 1, sum)[1:15]
moran.plot(y, rwm)

