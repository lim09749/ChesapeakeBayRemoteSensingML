##################################################################################################################################
## PROGRAM: extract and clean data.R
## SPECIFICATIONS: Extract and clean data
#####################################################################################################################################

rm(list=ls())

require(rgdal)
#####################################################################################################################################
## Load in data
setwd("D:/Publication/")
#filename <- "TSS and chl all stations Jul2002 to Dec2016 cleaned with MODISA.txt"
filename <- "TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt"

#####################################################################################################################################
#Load TSS/CHL with MODIS Data

CB  <- read.table(filename, header=T, strip.white = F, sep="\t")
CB  <- na.omit(CB)

#Convert Distance to meters by converting lat lon to UTM 
convertToMeters <- function(lon, lat){
  xy              <- data.frame(X = c(lon), Y = c(lat))#,ID = c(EventId))
  
  coordinates(xy) <- c("X", "Y")
  
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  
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
#Find closest point from shore to satellite
nearestInd <- apply(coordSat, MARGIN =c(1), FUN = function(x) 
  which.min(rowSums(sweep(coordCB,2,x)^2)))
coordCoast <- coordCB[nearestInd,]

dist_sat_coast <- sqrt(rowSums((coordCoast-coordSat)^2))


# .coord = project(cbind(CB$MODIS.lon[CB$EventId%in%outliers],
#                        CB$MODIS.lat[CB$EventId%in%outliers]),
#                  "+proj=utm +zone=18 ellps=WGS84")
# plot(coordCB[,1], coordCB[,2], col = "green", cex = 0.1)
# points(.coord[,1], .coord[,2], col = "red")

######################################################################################################################
##Combine info
CB <- cbind(CB, dist_sat_station, dist_sat_coast)

##Get rid of negative values
# cols_names <- c("chl", "Rrs_412", "Rrs_443", "Rrs_469", "Rrs_488", "Rrs_531",
#                 "Rrs_547", "Rrs_555","Rrs_645","Rrs_667","Rrs_678")
cols_names <- c("TSS", "Rrs_412", "Rrs_443", "Rrs_469", "Rrs_488", "Rrs_531",
                "Rrs_547", "Rrs_555","Rrs_645","Rrs_667","Rrs_678")

##Find relevant columns
relv_cols <- which(colnames(CB)%in%cols_names)

##Remove negative values
CB <- CB[!apply(CB[,relv_cols], MARGIN=1,FUN=function(x) { any(x < 0)}),]

##Adjust colnames
CB <- cbind(CB[,which(cols_names[1]==colnames(CB))],
               CB[,-which(cols_names[1]==colnames(CB))]) 
colnames(CB)[1] = cols_names[1]

##Clean data
CBnew <- CB[(CB$dist_sat_station<500)&(CB$dist_sat_coast>1000),]

##Save CBnew
#write.csv(CBnew, "CB_chl.csv")
write.csv(CBnew, "CB_tss.csv")
