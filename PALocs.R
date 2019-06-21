setwd("C:/Users/rvais/Wdirectory/SPA2/Midterm/PALocs")
getwd()
dir(getwd())
#Loading needed libraries
library(sp)
library(rgdal)
library(raster)
library(spatstat)
library(sp)
library(maptools)
#Reading the shapefile data and assign it to a variable
locPA <- readOGR(dsn= getwd(), layer = "PALocs")
#finding the min and max co-ordinates of the shape file(x1min,x2min,x1max,x2max)
x1min=min(locPA@coords[,1])
x2min=min(locPA@coords[,2])
x1max=max(locPA@coords[,1])
x2max=max(locPA@coords[,2])

#Store coordinates as list and make matrix for the coordinates
coord=list(x=locPA@coords[,1], y=locPA@coords[,2])
matx=matrix(unique(c(locPA@coords[,1],locPA@coords[,2])), ncol=2)

#Create an object of class "ppp" representing a point pattern dataset in the two-dimensional plane.
locPA_owin=as.owin(c(x1min,x1max,x2min,x2max))
locPA_unique=unique(locPA@coords)
locPA_ppp=ppp(locPA_unique[,1], locPA_unique[,2], locPA_owin)

#Performing quadrat sampling for regular and random sampling approaches

#Setting the dimensions of the quadrat
dim_x=30
dim_y=15

#Performing Random Quadrat sampling method
plot(locPA_ppp, pch=20, main="Random Quadrat Sampling") #Overlay points
legend(x1min, x2min, legend="PALocs", col="black", pch=20, cex=.75)
random_counts=sapply(1:(dim_x*dim_y), function(i){
  set.seed(i) #random quadrat samplings which are unique
  quadrat_w=(x1max-x1min)/dim_x #width of random quadrat 
  quadrat_h=(x2max-x2min)/dim_y #height of random quadrat 
  random_xaxis=runif(1, min=x1min, max=(x1max-quadrat_w)) #random X-coordinate 
  random_yaxis=runif(1, min=x2min, max=(x2max-quadrat_h)) #random Y-coordinate 
  locPArand_owin=as.owin(c(random_xaxis, (random_xaxis+quadrat_w), #create data window
                           random_yaxis, (random_yaxis+quadrat_h)))
  df=unique(data.frame(x=locPA@coords[,1],y=locPA@coords[,2])) #dataframe
  data_area=df[(df$x>random_xaxis)&(df$x<random_xaxis+quadrat_w) # points within random quadrat
               &(df$y>random_yaxis)&(df$y<random_yaxis+quadrat_h),]
  ppp_area_random=ppp(data_area[,1], data_area[,2],locPArand_owin) #reformatting the data to PPP format
  rect(xleft=random_xaxis, ybottom=random_yaxis, xright=(random_xaxis+quadrat_w), ytop=(random_yaxis+quadrat_h))
  quadrat_random=ifelse(is.na(data_area[1,1])==1, 0,
                        quadratcount(ppp_area_random, nx=1, ny=1)) #Points within the random quadrat are counted 
  quadrat_random
  return(quadrat_random)
})
#Calculating the mean
mu_random=mean(random_counts)
mu_random
total_quad_random=sum((random_counts-mu_random)^2)#total
#Calculating the variance
variance_random=total_quad_random/((dim_x*dim_y)-1) 
variance_random
#Calculating the VMR
VMR_random=variance_random/mu_random 
VMR_random



#Generating the table with attributes K,X,k-mu,(K-mu)^2,X(K-mu)^2

K_random=as.numeric(names(table(c(random_counts))))# Number of events K
X_random=c(table(c(random_counts)))#Number of quadrats X
table_random=data.frame(K_random,X_random)#Create data frame for the table with K and X
rownames(table_random)=NULL #Setting the fisrt row name to NULL for easy data collection
table_random$K_mu=table_random$K_random-mu_random #K-mu
table_random$K_mu2=table_random$K_mu^2 #(k-mu)^2
table_random$XK_mu2=table_random$X_random*table_random$K_mu2 #X(K-mu)^2
random_total=sum(table_random$XK_mu2)#total
write.csv(table_random,"C:/Users/rvais/Wdirectory/SPA2/Midterm/PAloc1")

G_FUNC=Gest(locPA_ppp, correction="none")
plot(G_FUNC, raw~r, ylim=c(0,1), main="PALoc G-Function",xlab=" distance d")



F_FUNC=Fest(locPA_ppp, correction="none")
plot(F_FUNC, raw~r, ylim=c(0,1), main="PALoc F-Function",xlab="distance d")
