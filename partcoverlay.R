#PartC
library(raster)
library(rgdal)
library(phylin)
library(sp)
library(spatstat)
library(maptools)
library(GISTools)
library(plyr)
library(gstat)
library(automap)
library(rgeos)

p1=readOGR("C:/Users/rvais/Wdirectory/Finals", "3_Area1")
p2=readOGR("C:/Users/rvais/Wdirectory/Finals", "3_Area2")
p3=readOGR("C:/Users/rvais/Wdirectory/Finals", "3_TrailPoints")
p4=readOGR("C:/Users/rvais/Wdirectory/Finals", "3_State_Roads")


plot(p1,main="3Area1",col="grey",cex=0.8)
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
plot(p2,main="3Area2",col="orange",cex=0.8)
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
plot(p3,main="3_trailpoints",cex=0.8)
plot(p4,main="3_stateroads",cex=0.8)


#checking for the projetion of the 2 files to verify if theyare in the same co-ordinate system
p1=geometry(p1)
class(p1)
projection(p1)

p2=geometry(p2)
class(p2)
projection(p2)

####Yes,they belong to same CRS so we can overlay plot p2 on p1

##a
over(p2,p1)#Overlay p2 on p1

p12= intersect(p2, p1)
p12
sum(poly.areas(p12)) #1.132137
plot(p12,col="yellow",main="Intersection of p2 on p1 ")
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
p11=p1-p12 #difference of plot1-plot12
plot(p11,col="pink",main="p1-p12")
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242246, yb=480342, len=5000, lab="N",col='Red')
p22=p2-p12
plot(p22,col="light blue",main="p2-p12")
maps::map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
##bMap showing points
p3=geometry(p3)
class(p3)
projection(p3)# projection is different from that of p1,p2, and p12

p12=geometry(p12)
class(p12)
projection(p12)
summary(p12)


library(rgdal)
TA=CRS("+proj=longlat +datum=NAD27 +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat")
p3TA=spTransform(p3,TA)
projection(p3TA)
over(p3TA,p1)
over(p3TA,p2)
over(p3TA,p12)
bbox(p3TA)
head(p3TA@coords)
p31= intersect(p3TA, p1)
p31#428
sum(poly.counts(p31,p1))
p32= intersect(p3TA, p2)
p32#566
sum(poly.counts(p32,p2))
p33= intersect(p3TA, p12)
p33#281
sum(poly.counts(p33,p2))
par(mfrow=c(1,3))
plot(p1,col="yellow"); points(p31,col="red"); title(main = list("Plot of p3 points on p1", cex=0.8))
plot(p2,col="pink"); points(p32,col="purple"); title(main = list("Plot of p3 points on p2", cex=0.8))
plot(p12,col="light blue"); points(p33,col="blue"); title(main = list("Plot of p3 points on p12", cex=0.8))

#Map showing the lines
p4=geometry(p4)
class(p4)
projection(p4)# projection is different from that of p1,p2, and p12
p4TA=spTransform(p4,TA)#Transforming to the Co-ordinate system of p1,p2,p12
projection(p4TA)
over(p4TA,p1)
over(p4TA,p2)
over(p4TA,p12)
bbox(p4TA)
head(p4TA@lines)
p41= intersect(p4TA, p1)
p41#15358
p42= intersect(p4TA, p2)
p42#20453
p43= intersect(p4TA, p12)
p43#8581
par(mfrow=c(1,3))
plot(p1,col="yellow"); lines(p41,col="red"); title(main = list("Plot of p41 lines on p1", cex=0.8))
plot(p2,col="pink"); lines(p42,col="purple"); title(main = list("Plot of p42 lines on p2", cex=0.8))
plot(p12,col="light blue"); lines(p43,col="blue"); title(main = list("Plot of p43 points on p12", cex=0.8))
