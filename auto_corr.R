#Part-A(Morans'I technique implementation)

#Read the datafile
particlem=read.csv("http://gis40.exp.sis.pitt.edu/INFSCI2809_data/ParticulateMatter.csv", sep=",", header=T)
#reading longitude and latitude as distance matrix
particle_dists=as.matrix(dist(cbind(particlem$Lon, particlem$Lat))) 
#Taking inverse distance matrix
particle_dists_inv=1/particle_dists
#setting the diagonal to zero
diag(particle_dists_inv)=0 
#read obseravations
no_of_obs=length(particlem$PM25) 
#divide inverse matrix by rowsum
particle_dists_inv=particle_dists_inv/rowSums(particle_dists_inv)
#sum of inverse distance matrix
sum_of_inv=sum(particle_dists_inv) 
#mean of PM25 attribute
m=mean(particlem$PM25)
#difference of attribute from mean
y=particlem$PM25-m 
#sum of squared differences
v=sum(y^2) 
#Perform outer cross product
z=sum(particle_dists_inv*y %o% y)
#Apply autocorrelation using the formula 
autocorrelation=(no_of_obs/sum_of_inv)*((z)/v) 
print(autocorrelation)
