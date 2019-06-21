#Part-B
library(sp)
library(ggplot2)
library(rgeos)
library(maptools)
library(ggplot2)
library(ggmap)
library(maps)
library(GISTools)
library(spdep)
#Loading the data
crime_data=readOGR("C:/Users/rvais/Wdirectory/PROJ3", "Crime_PA2002")
n=length(crime_data$COUNTY)
crime_centroid=gCentroid(crime_data, byid=T, id=crime_data@data$COUNTY)

#Plot of the crime dataset
plot(crime_data)
par(mai=c(0,0,0.2,0))
plot(crime_data, col='gray', border='blue')
map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
map.axes(cex.axis=0.8)
title('Crime dataset with Mifflin County')
points(crime_centroid[61,], col="blue", pch=20)
legend('bottomleft',legend="Mifflin County",pch=20, col="blue", cex=0.75)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')


#Heatmap of crime against each of population,agencies and county
spplot(aggregate(x=crime_data[c('INDEX01')], FUN=mean, by=crime_data), main="Crime by County")
spplot(aggregate(x=crime_data[c('POP_CRI01')], FUN=mean, by=crime_data), main="Crime by Population")
spplot(aggregate(x=crime_data[c('AG_CRI01')], FUN=mean, by=crime_data), main="Crime by Number of Crime Agencies")

#Global G-statistics
crime_distance=as.matrix(dist(cbind(crime_centroid@coords[,1], crime_centroid@coords[,2]))) #distance matrix
crime_inv=1/crime_distance #inverse distance matrix
diag(crime_inv)=0 #set diagonals to zero
crime_inv=crime_inv/rowSums(crime_inv) #divide crime_inv by row sums (maybe delete)
s=sum(crime_inv) #sum of inverse distance matrix
crime_i=matrix(rep(crime_data@data$INDEX01, length(crime_data@data$INDEX01)), nrow=length(crime_data@data$INDEX01))
diag(crime_i)=0
crime_j=t(crime_i)
product=crime_i*crime_j
global_gstatistics=sum(crime_inv*product)/sum(product)
global_gstatistics


#GWR technique
mifflin_xcord=c(rep(crime_centroid@coords[61,1], n))
mifflin_ycord=c(rep(crime_centroid@coords[61,2], n))
weights=sqrt((mifflin_xcord-crime_centroid@coords[,1])^2+(mifflin_ycord-crime_centroid@coords[,2])^2)#weighted
distance_matrix=matrix(0, nrow=n, ncol=n)
diag(distance_matrix)=weights
x=as.matrix(data.frame(crime_data@data$POP_CRI01, crime_data@data$AG_CRI01, crime_data@data$Area))
x_transpose=t(x)
y=crime_data@data$INDEX01
beta_value=solve(x_transpose%*%distance_matrix%*%x)%*%x_transpose%*%distance_matrix%*%y
beta_value#Gives beta values
Intercept_value=coefficients(lm(data=crime_data, INDEX01~POP_CRI01+AG_CRI01+Area))[1]
Intercept_value#Gives intercept
