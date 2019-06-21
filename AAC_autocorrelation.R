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


#Set the working directory and get the directory

Allegh_council=readOGR("C:/Users/rvais/Wdirectory/PROJ3", "AlleghenyCounty_Council")
Allegh_municipal=readOGR("C:/Users/rvais/Wdirectory/PROJ3", "AlleghenyCounty_Municipal")


#Representing shape file as Rook adjacency matrix
rooks <- poly2nb(Allegh_council, row.names=Allegh_council$District, queen=FALSE)
class(rooks)
summary(rooks)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Allegh_council, col='gray', border='black')
map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Allegheny County Council')
legend('bottomleft', legend = 'Council Lines', lty = 1, lwd = 2, bty = 'o')

xy <- coordinates(Allegh_council)
plot(rooks, xy, col='black', lwd=2, add=TRUE)
n=length(council@data$District)
rooks_matrix=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    rooks_matrix[i,j]=ifelse(j %in% rooks[[i]], 1, 0)
  }
}
expected=-1/(n-1)
s=sum(rooks_matrix) #sum of weight matrix
m=mean(council@data$SHAPE_area)#mean of attribute
y=(council@data$SHAPE_area)-m #difference from mean
v=sum(y^2) #sum of squared differences
moran_rooks=(n/s)*(sum(rooks_matrix*y %o% y)/v) #Moran's I
moran_rooks
S1=0.5*sum((rooks_matrix+t(rooks_matrix))^2)
S2=sum((apply(rooks_matrix, 1, sum)+apply(rooks_matrix, 2, sum))^2)
s.sq=s^2
k=(sum(y^4)/n)/(v/n)^2
sdi=sqrt((n*((n^2-3*n+3)*S1-n*S2+3*s.sq)-
            k*(n*(n-1)*S1-2*n*S2+6*s.sq))/((n-1)*(n-2)*(n-3)*s.sq)-1/((n-1)^2))
p_value=pnorm(moran_rooks, mean = expected, sd = sdi)
if(moran_rooks<=expected)
{p_value=2*p_value
}else
{ p_value=2*(1-p_value)
}
p_value

#Council with queens
#Representing shape file as queens adjacency matrix
queens <- poly2nb(Allegh_council, row.names=Allegh_council$District, queen=TRUE)
class(queens)
summary(queens)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Allegh_council, col='gray', border='black')
map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Allegheny County Council')
legend('bottomleft', legend = 'Council Lines', lty = 1, lwd = 2, bty = 'o')

xy <- coordinates(Allegh_council)
plot(queens, xy, col='black', lwd=2, add=TRUE)
n=length(council@data$District)
queens_matrix=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    queens_matrix[i,j]=ifelse(j %in% rooks[[i]], 1, 0)
  }
}
expected=-1/(n-1)
s=sum(queens_matrix) #sum of weight matrix
m=mean(council@data$SHAPE_area) #mean of attribute
y=(council@data$SHAPE_area)-m #difference from mean
v=sum(y^2) #sum of squared differences
moran_queens=(n/s)*(sum(queens_matrix*y %o% y)/v) #Moran's I
moran_queens
S1=0.5*sum((queens+t(queens))^2)
S2=sum((apply(queens, 1, sum)+apply(queens, 2, sum))^2)
s.sq=s^2
k=(sum(y^4)/n)/(v/n)^2
sdi=sqrt((n*((n^2-3*n+3)*S1-n*S2+3*s.sq)-
            k*(n*(n-1)*S1-2*n*S2+6*s.sq))/((n-1)*(n-2)*(n-3)*s.sq)-1/((n-1)^2))
p_value=pnorm(moran_queens, mean = expected, sd = sdi)
if(moran_queens<=expected)
{
  p_value=2*p_value
}else
{
  p_value=2*(1-p_value)
}
p_value

#Municipal+Rooks

#Representing shape file as Rook adjacency matrix
rooks_municipal <- poly2nb(Allegh_municipal, row.names=Allegh_municipal@data$OBJECTID, queen=FALSE)
class(rooks_municipal)
summary(rooks_municipal)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Allegh_municipal, col='gray', border='black')
map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Allegheny County Municipal')
legend('bottomleft', legend = 'Municipal Lines', lty = 1, lwd = 2, bty = 'o')

xy <- coordinates(Allegh_municipal)
plot(rooks_municipal, xy, col='black', lwd=2, add=TRUE)
n=length(Allegh_municipal@data$OBJECTID)
rooks_matrix_municipal=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    rooks_matrix_municipal[i,j]=ifelse(j %in% rooks_municipal[[i]], 1, 0)
  }
}
expected=-1/(n-1)
s=sum(rooks_matrix_municipal) #sum of weight matrix
m=mean(Allegh_municipal@data$SHAPE_area) #mean of attribute
y=(Allegh_municipal@data$SHAPE_area)-m #difference from mean
v=sum(y^2) #sum of squared differences
moran_rooks_municipal=(n/s)*(sum(rooks_matrix_municipal*y %o% y)/v) #Moran's I
moran_rooks_municipal
S1=0.5*sum((rooks_matrix_municipal+t(rooks_matrix_municipal))^2)
S2=sum((apply(rooks_matrix_municipal, 1, sum)+apply(rooks_matrix_municipal, 2, sum))^2)
s.sq=s^2
k=(sum(y^4)/n)/(v/n)^2
sdi=sqrt((n*((n^2-3*n+3)*S1-n*S2+3*s.sq)-
            k*(n*(n-1)*S1-2*n*S2+6*s.sq))/((n-1)*(n-2)*(n-3)*s.sq)-1/((n-1)^2))
p_value=pnorm(moran_rooks_municipal, mean = expected, sd = sdi)
if(moran_rooks_municipal<=expected)
{p_value=2*p_value
}else
{ p_value=2*(1-p_value)
}
p_value

#Council with queens
#Representing shape file as queens adjacency matrix
queens_municipal <- poly2nb(Allegh_municipal, row.names=Allegh_municipal@data$OBJECTID, queen=TRUE)
class(queens_municipal)
summary(queens_municipal)
#Plotting the link between the polygons
par(mai=c(0,0,0.2,0))
plot(Allegh_municipal, col='gray', border='black')
map.scale(y = 360000, ratio=FALSE, relwidth=.2, metric=FALSE)
north.arrow(xb=1242280, yb=480342, len=5000, lab="N",col='Red')
map.axes(cex.axis=0.8)
title('Allegheny County Municipal')
legend('bottomleft', legend = 'Municipal Lines', lty = 1, lwd = 2, bty = 'o')
xy <- coordinates(Allegh_municipal)
plot(queens_municipal, xy, col='black', lwd=2, add=TRUE)
n=length(Allegh_municipal@data$OBJECTID)
queens_matrix_municipal=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
 for(j in 1:n){
   queens_matrix_municipal[i,j]=ifelse(j %in% queens_municipal[[i]], 1, 0)
 }
}
expected=-1/(n-1)
s=sum(queens_matrix_municipal) #sum of weight matrix
m=mean(Allegh_municipal@data$SHAPE_area) #mean of attribute
y=(Allegh_municipal@data$SHAPE_area)-m #difference from mean
v=sum(y^2) #sum of squared differences
moran_queens_municipal=(n/s)*(sum(queens_matrix_municipal*y %o% y)/v) #Moran's I
moran_queens_municipal
S1=0.5*sum((queens_matrix_municipal+t(queens_matrix_municipal))^2)
S2=sum((apply(queens_matrix_municipal, 1, sum)+apply(queens_matrix_municipal, 2, sum))^2)
s.sq=s^2
k=(sum(y^4)/n)/(v/n)^2
sdi=sqrt((n*((n^2-3*n+3)*S1-n*S2+3*s.sq)-
            k*(n*(n-1)*S1-2*n*S2+6*s.sq))/((n-1)*(n-2)*(n-3)*s.sq)-1/((n-1)^2))
p_value=pnorm(moran_queens_municipal, mean = expected, sd = sdi)
p_value=if(moran_queens_municipal<=expected)
{
  2*p_value
}else 
  {2*(1-p_value)
  }
p_value
#Geary's
geary_rooks=poly2nb(Allegh_council, row.names=Allegh_council@data$District, queen=F) #rooks adjacency
n=length(Allegh_council@data$District)
rooks_matrix_gearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
  for(i in 1:n){
    for(j in 1:n){
      rooks_matrix_gearyC[i,j]=ifelse(j %in% geary_rooks[[i]], 1, 0)
    }
  }
i_matrix=matrix(rep(Allegh_council@data$SHAPE_area, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Allegh_council@data$SHAPE_area)
geary=((n-1)/sum((i_matrix-mean)^2))*sum(rooks_matrix_gearyC*(diff_matrix)^2)/(2*sum(rooks_matrix_gearyC)) #Geary's C
geary

#Using Queens adjacency
geary_queens=poly2nb(Allegh_council,row.names=Allegh_council@data$District, queen=T) #queen's adjacency
n=length(Allegh_council@data$District)
queens_matrix_GearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
  for(i in 1:n){
    for(j in 1:n){
      queens_matrix_GearyC[i,j]=ifelse(j %in% geary_queens[[i]], 1, 0)
    }
  }
i_matrix=matrix(rep(Allegh_council@data$SHAPE_area, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Allegh_council@data$SHAPE_area)
geary_queens_county=((n-1)/sum((i_matrix-mean)^2))*sum(queens_matrix_GearyC*(diff_matrix)^2)/(2*sum(queens_matrix_GearyC)) #Geary's C
geary_queens_county

#Geary's Municipal
geary_rooks_municipal=poly2nb(Allegh_municipal, row.names=Allegh_municipal@data$OBJECTID, queen=F) #rooks adjacency
n=length(Allegh_municipal@data$OBJECTID)
rooks_matrix_gearyC_M=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    rooks_matrix_gearyC_M[i,j]=ifelse(j %in% geary_rooks_municipal[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Allegh_municipal@data$SHAPE_area, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Allegh_municipal@data$SHAPE_area)
geary_municipal=((n-1)/sum((i_matrix-mean)^2))*sum(rooks_matrix_gearyC_M*(diff_matrix)^2)/(2*sum(rooks_matrix_gearyC_M)) #Geary's C
geary_municipal


#queens
geary_queens_municipal=poly2nb(Allegh_municipal, row.names=Allegh_municipal@data$OBJECTID, queen=T) #rooks adjacency
n=length(Allegh_municipal@data$OBJECTID)
queens_matrix_gearyC_M=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    queens_matrix_gearyC_M[i,j]=ifelse(j %in% geary_queens_municipal[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Allegh_municipal@data$SHAPE_area, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
geary_municipal_queen=((n-1)/sum((i_matrix-mean)^2))*sum(queens_matrix_gearyC_M*(diff_matrix)^2)/(2*sum(queens_matrix_gearyC_M)) #Geary's C
geary_municipal_queen
