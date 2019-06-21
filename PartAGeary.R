#Part A-Gearys'C
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

Neighbor1=readOGR("C:/Users/rvais/Wdirectory/Finals","1_Neighbor1")
Neighbor2=readOGR("C:/Users/rvais/Wdirectory/Finals", "1_Neighbor2")

#Compute Gearys'C
geary_rooks=poly2nb(Neighbor1, row.names=Neighbor1@data$COUNTY, queen=F) #rooks adjacency
n=length(Neighbor1@data$STATE)
rooks_matrix_gearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    rooks_matrix_gearyC[i,j]=ifelse(j %in% geary_rooks[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Neighbor1@data$POP_ARR02, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Neighbor1@data$POP_ARR02)
geary=((n-1)/sum((i_matrix-mean)^2))*sum(rooks_matrix_gearyC*(diff_matrix)^2)/(2*sum(rooks_matrix_gearyC)) #Geary's C
geary

#GearyC queen neighbor1
geary_queens=poly2nb(Neighbor1,row.names=Neighbor1@data$COUNTY, queen=T) #queen's adjacency
n=length(Neighbor1@data$STATE)
queens_matrix_GearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    queens_matrix_GearyC[i,j]=ifelse(j %in% geary_queens[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Neighbor1@data$POP_ARR02, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Neighbor1@data$POP_ARR02)
geary_queens_neighbor1=((n-1)/sum((i_matrix-mean)^2))*sum(queens_matrix_GearyC*(diff_matrix)^2)/(2*sum(queens_matrix_GearyC)) #Geary's C
geary_queens_neighbor1

#Gearys'C rooks Neighbor 2
geary_rooks_neighbor2=poly2nb(Neighbor2, row.names=Neighbor2@data$COUNTY, queen=F) #rooks adjacency
n=length(Neighbor2@data$STATE)
rooks_matrix_gearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    rooks_matrix_gearyC[i,j]=ifelse(j %in% geary_rooks_neighbor2[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Neighbor2@data$POP_ARR02, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Neighbor2@data$POP_ARR02)
geary_neighbor2=((n-1)/sum((i_matrix-mean)^2))*sum(rooks_matrix_gearyC*(diff_matrix)^2)/(2*sum(rooks_matrix_gearyC)) #Geary's C
geary_neighbor2



#Gearys'C Neighbor 2 Queens
geary_queens_neighbor2=poly2nb(Neighbor2,row.names=Neighbor2@data$COUNTY, queen=T) #queen's adjacency
n=length(Neighbor2@data$STATE)
queens_matrix_GearyC=matrix(rep(0,(n^2)), nrow=n, ncol=n) #fill adjacency matrix
for(i in 1:n){
  for(j in 1:n){
    queens_matrix_GearyC[i,j]=ifelse(j %in% geary_queens_neighbor2[[i]], 1, 0)
  }
}
i_matrix=matrix(rep(Neighbor2@data$POP_ARR02, n), nrow=n) #attribute matrix
j_matrix=t(i_matrix) #transpose of attribute matrix
diff_matrix=i_matrix-j_matrix #yi-yj
mean=mean(Neighbor2@data$POP_ARR02)
geary_queens_neighbor2=((n-1)/sum((i_matrix-mean)^2))*sum(queens_matrix_GearyC*(diff_matrix)^2)/(2*sum(queens_matrix_GearyC)) #Geary's C
geary_queens_neighbor2


