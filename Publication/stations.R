##################################################################################################################################
## PROGRAM: stations.R
## SPECIFICATIONS: Plots the location of sampling stations
#####################################################################################################################################
#### Working Directory
rm(list=ls())
setwd("D:/Publication")

#####################################################################################################################################
#### Libraries
require(PBSmapping)
require(sp)
require(raster)
######################################################################################################################
#### EDIT THIS

# Main Plot Coordinates 

Nlat <- 39.7 #max(Sites$Y)+0.2
Slat <- 36.7 #min(Sites$Y)-0.2
Elon <- 284.4#max(Sites$X)+0.2
Wlon <- 282.5 #min(Sites$X)-0.2

######################################################################################################################

# Import base and inset maps
zoomMap <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), maxLevel=3)
#outMap  <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(36.7, 39.7), maxLevel=3)

# Import and prepare necessary points 
sites = read.csv(file="CB_tss.csv")[,c(5,6)]#Load data
sites$X <- sites$lon %% 360 #Make longitude positive, and load in value
sites$Y <- sites$lat #Load latitude
sites$EID <- 1:nrow(sites) #Random number as event ID
attr(sites, "projection")="LL" #Add number for projection type
sites = as.EventData(sites) #Create event ID

# HPL = rbind(c(1, (-76.1347496 %% 360), 38.5830642),c(2, (-76.1347496 %% 360), 38.5830642)) #Horn Point Lab location
# colnames(HPL) = c("EID", "X", "Y") #Change to right column names
# attr(HPL, "projection")="LL" #Have to create attribute "projection"
# HPL = as.EventData(HPL) #Create Event ID

# Parameterize and plot main map, add points listed above
tiff(file = "stations.tiff", width =3.75, height = 7, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(fig=c(0,1,0,1))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", ylab = "Latitude") 
addPoints(sites, col=2, pch=19, xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), lwd=5,cex=0.3)
#addPoints(HPL, col=4, pch=4, xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), lwd=5)
# Add necessary labels, legends, scale bar, compass rose

#text(283.957, 38.64, substitute(paste(italic("Choptank River"))), cex=.9) 
#text(283.56, 38.7, substitute(paste(italic("Chesapeake Bay"))), cex=1) 
legend("topleft", legend=c("Sampling Locations"), col=c(2), 
       lty=0, lwd=3, pch=c(19),cex=0.70)
scalebar(location="bottomright", y.min=Slat, y.max=Nlat, x.min=Wlon, x.max=Elon, dist=0.5, divs=2, type='bar', 
         lonlat=T, below="km")
dev.off()

#compassRose(282.9, 37.75, rot=0, cex=.9)
######################################################################################################################

# # Parameterize and plot inset map, replot border lines
# 
# par(fig=c(.745, 1, .1, .58), new=T, mar=c(0,0,0,0))
# plotMap(outMap, col="ivory2", bg="lightblue3", xlab="", ylab="", axe=F, frame.plot=T) 
# lines(c(283.4, 285), c(37.14, 37.14))
# lines(c(285, 285), c(37.1 39))
# lines(c(283.4, 285), c(39, 39))
# 
