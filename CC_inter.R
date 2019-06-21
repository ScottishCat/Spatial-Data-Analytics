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
#Loading the dataset
ozone_sensor=readOGR("C:/Users/rvais/Wdirectory/PROJ3", "Ozone_Sensor_Locs")
ozone_value=read_delim("C:/Users/rvais/Wdirectory/PROJ3/Ozone_Value.dat", "|",
                       escape_double=FALSE, trim_ws=TRUE)
ozone_val=ozone_value[ozone_value$parameter=="OZONE",]
ozone_val$id=ozone_val$id
pa_county=readOGR("C:/Users/rvais/Wdirectory/PROJ3", "PA_County_Select")

# Finding centroid_pacountys
centroid_pacounty=getSpPPolygonsLabptSlots(pa_county)
plot(pa_county)
points(ozone_sensor, col="Red")
points(centroid_pacounty, col="Black")

# Datasets are merged
merge_ozone=merge(ozone_sensor, ozone_val, by="id") #all.x=TRUE)
merge_ozone=merge_ozone[is.na(merge_ozone$value)==FALSE,]
new_coords.df=data.frame(x=centroid_pacounty[,1], y=centroid_pacounty[,2])
old_coords.df=data.frame(x=merge_ozone@coords[,1], y=merge_ozone@coords[,2])

# Inverse distance weight
distance_real=geo.dist(from=new_coords.df, to=old_coords.df) #distances
w=1/distance_real^2 # Inverse distance weights d=1,3 also tried
for(i in 1:10){
  row_median=median(w[i,])
  for(j in 1:11){
    w[i,j]=ifelse(w[i,j]>row_median, w[i,j], 0) #Five closest sensors
  }
}
weight_sum=rowSums(w) # Weight matrix 
weight_matrix=w %*% diag(merge_ozone$value)
idw_matrix=apply(weight_matrix/weight_sum, 1, sum, na.rm=TRUE) # Divide by rowsums
idw_values=idw_matrix
table_idw=data.frame(County=pa_county$COUNTY, idw_values) # IDW results table

#Interpolated surface plot
sp=SpatialPoints(data.frame(x=centroid_pacounty[,1], y=centroid_pacounty[,2]), proj4string=CRS("+proj=longlat +datum=NAD83"))
dsp=SpatialPointsDataFrame(sp, data.frame(x=centroid_pacounty[,1], y=centroid_pacounty[,2], value=idw_values))
dta=spTransform(dsp, CRSobj="+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
cata=spTransform(pa_county, CRSobj="+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
r=raster(cata, res=100)
ca=aggregate(cata)
library(gstat)
gs_idw=gstat(formula=value~1, locations=dta)
idw_interpolation=interpolate(r, gs_idw, xyOnly=T, xyNames=c("x","y"))
idw_result=mask(x=idw_interpolation, mask=ca)
plot(idw_result, main="IDW Interpolation of Ozone Levels")

#Applying semiveriogram technique and and plotting Interpolated surface  
K_variogram=variogram(log(value)~1, merge_ozone) # calculates sample variogram values 
variogram_fit=fit.variogram(K_variogram, vgm(psill=0.012, range=50, nugget=0, "Exp"))
plot(K_variogram, variogram_fit)
sp2=SpatialPoints(data.frame(x=centroid_pacounty[,1], y=centroid_pacounty[,2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
OrdinaryK=krige(log(value)~1, merge_ozone, newdata=sp2, model=variogram_fit)
ok_values=exp(OrdinaryK@data$var1.pred)
table_Ok=data.frame(County=pa_county$COUNTY, ok_values) # OK results table
dsp2=SpatialPointsDataFrame(sp2, data.frame(x=centroid_pacounty[,1], y=centroid_pacounty[,2], value=ok_values))
dta2=spTransform(dsp2, CRSobj="+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
cata=spTransform(pa_county, CRSobj="+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
r=raster(cata, res=100)
ca=aggregate(cata)
geo_kringing=gstat(formula=value~1, locations=dta2)
OK_interpolation=interpolate(r, geo_kringing, xyOnly=T, xyNames=c("x","y"))
OK_result=mask(x=OK_interpolation, mask=ca)
plot(OK_result, main="Ordinary Kriging Interpolation Technique")
