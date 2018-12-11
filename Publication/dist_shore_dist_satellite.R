##################################################################################################################################
## PROGRAM: dist_shore_dist_satellite.R
## SPECIFICATIONS: Se how distance between shore and satellite and satellite and in-situ affects r-squared
#####################################################################################################################################

rm(list=ls())

require(rgdal)
#####################################################################################################################################
#Load TSS/CHL with MODIS Data

filename <- "D:/Publication/TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt"#
#filename <- "D:/Publication/TSS and chl all stations Jul2002 to Dec2016 cleaned with MODISA.txt"

CB  <- read.table(filename, header=T, strip.white = F, sep="\t")
CB  <- na.omit(CB)

#Convert Distance to meters by converting lat lon to UTM 
convertToMeters <- function(lon, lat){
  xy              <- data.frame(X = c(lon), Y = c(lat))#,ID = c(EventId))
  
  coordinates(xy) <- c("X", "Y")
  
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  
  ##Transform to common data format
  coordStation             <- as.data.frame(spTransform(xy, CRS("+proj=utm +zone=18 ellps=WGS84")))
  #coordStation$ID <- NULL
  
  coordStation <-as.matrix(coordStation, ncol=2) #Coordinates of TSS stations
  return(coordStation)
}


#Station and Sat: Convert Distance to meters by converting lat lon to UTM
coordStation <- convertToMeters(CB$lon, CB$lat)
coordSat <- convertToMeters(CB$MODIS.lon, CB$MODIS.lat)

#Find distance
dist_sat_station <- sqrt(rowSums((coordStation-coordSat)^2))
######################################################################################################################
#Open Shapefile
setwd("D:/Publication/SatData/Baybox CB Shoreline")
shp = "baybox_dd.shp"
myshp = readOGR(shp, layer = basename(strsplit(shp, "\\.")[[1]])[1]) # This is a fancy way of being lazy, so I do not need to type the layer name in

#Dump shapefile coordinates into 'coords'
polys = attr(myshp,'polygons')
npolys = length(polys)
for (i in 1:npolys){
  poly = polys[[i]]
  polys2 = attr(poly,'Polygons')
  npolys2 = length(polys2)
  for (j in 1:npolys2){
    #do stuff with these values
    if (i == 1 & j==1) {
      coords =  coordinates(polys2[[j]])
    }else{
      coords =  rbind(coords, coordinates(polys2[[j]]))
    }
  }
}

#Convert coordinates to UTM
coordCB <- project(coords, "+proj=utm +zone=18 ellps=WGS84")  #Coordinates of figure
######################################################################################################################
#Find closest point to shore to satellite
nearestInd <- apply(coordSat, MARGIN =c(1), FUN = function(x) 
  which.min(rowSums(sweep(coordCB,2,x)^2)))
coordCoast <- coordCB[nearestInd,]

dist_sat_coast <- sqrt(rowSums((coordCoast-coordSat)^2))
######################################################################################################################

##Get rid of negativenot desired values
tssCol <- c(3, 10:19)
chlCol <- c(3, 10:20)

CB  <- CB[!apply(CB[,tssCol], MARGIN=1,FUN=function(x) { any(x < 0)}),]
#CB  <- CB[!apply(CB[,chlCol], MARGIN=1,FUN=function(x) { any(x < 0)}),]
#vec <- CB$chl
vec <- CB$TSS  
# Possible metrics to check
max_dist_sat_station <- seq(200, 500,  by=50)
min_dist_sat_coast   <- seq(400, 1200, by=50)

######################################################################################################################
# Sake of convenience
getRsquared <- function(vec, wavelength){
  return(summary(lm(vec~wavelength))$r.squared)
}

# Find accuracy (r2) from metrics
findAccuracy <- function(CB, vec, max_dist_sat_station, min_dist_sat_coast){
  CBnew <- CB[(dist_sat_station<max_dist_sat_station)&(dist_sat_coast>min_dist_sat_coast),]
  vec   <- vec[(dist_sat_station<max_dist_sat_station)&(dist_sat_coast>min_dist_sat_coast)]
  
  mat_rsquared   <- matrix(NA, nrow = 10)
  mat_wavelength <- as.matrix(cbind(CBnew$Rrs_412, CBnew$Rrs_443, CBnew$Rrs_469, 
                                    CBnew$Rrs_488, CBnew$Rrs_531, CBnew$Rrs_547, 
                                    CBnew$Rrs_555, CBnew$Rrs_645, CBnew$Rrs_667, CBnew$Rrs_678))

  for(i in 1:10){
    mat_rsquared[i] <- getRsquared(log10(vec), log10(mat_wavelength[,i]))
  }
  
  return(mat_rsquared)
}

# Create accuracy matrix
createAccuracyMat <- function(CB, vec, max_dist_sat_station, min_dist_sat_coast){
  CBmetrics <- matrix(NA, nrow=length(max_dist_sat_station), ncol=length(min_dist_sat_coast))
  
  for(sat_station in max_dist_sat_station){
    for(sat_coast in min_dist_sat_coast){
      CBmetrics[which(sat_station==max_dist_sat_station), 
                which(sat_coast==min_dist_sat_coast)] <-
        sum(findAccuracy(CB, vec, sat_station, sat_coast)) 
    }
  }
  return(CBmetrics)
}

# Flatten matrix and have corresponding x, y, and r^2
metrics <- as.vector(createAccuracyMat(CB, vec, max_dist_sat_station, min_dist_sat_coast))
x       <- rep(max_dist_sat_station, length(min_dist_sat_coast))
y       <- rep(min_dist_sat_coast, each=length(max_dist_sat_station))

write.csv(cbind(x,y,metrics), "D:/Publication/chl_dist_metrics.csv")
######################################################################################################################
## Plot contour
require(raster)

r_chl <- read.csv("D:/Publication/chl_dist_metrics.csv")
r_chl$X <- NULL
r_tss <- read.csv("D:/Publication/tss_dist_metrics.csv")
r_tss$X <- NULL

r_chl <- rasterFromXYZ(cbind(r_chl$x,r_chl$y,r_chl$metrics))
r_tss <- rasterFromXYZ(cbind(r_tss$x,r_tss$y,r_tss$metrics))

values(r_chl) <- values(r_chl)/max(values(r_chl))
values(r_tss) <- values(r_tss)/max(values(r_tss))

plot(r_chl)
plot(r_tss)
plot(r_chl+r_tss)
