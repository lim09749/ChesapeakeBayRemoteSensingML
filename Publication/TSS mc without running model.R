##################################################################################################################################
## PROGRAM: CHL mc without running models.R
## SPECIFICATIONS: Plots mc plots
#####################################################################################################################################

rm(list=ls())

require(rgdal)
require(raster)  # http://rpubs.com/etiennebr/visualraster
require(viridis)
require(rgeos)

################################################################################
setwd("D:/CB_MC_1km_from_L2") ##MC file location
load("D:/Publication/modis pretty.RData")

CBnew <- read.csv("D:/Publication/CB_tss.csv")
model_var <- c("TSS","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

##EDIT HERE
pattern <- c("rvm","svm", "glm", "ann")
mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
mymonths <- c("January","February","March", "April","May","June",
              "July","August","September", "October","November","December")

TSSmin = 0.01
TSSmax = 30

##Read in CB shapefile
basemap.r = raster("D:/GIS/r.tif")
basemap.g = raster("D:/GIS/g.tif")
basemap.b = raster("D:/GIS/b.tif")

basemap = stack(basemap.r, basemap.g, basemap.b)

# Function to plot data
plotMC <- function(MCraster, MClabel, TSSmin, TSSmax){
  par(mai=c(0,0,0,0)) # mai is the margin of individual plots
  # Blank plot
  plot(1,1, type="n", axes=F, xlab="", ylab="",
       xlim=c(-77, -75.6), ylim=c(36.8, 39.6))
  rect(-77, 36.8, -75.6, 39.6, col="grey50")
  plotRGB(basemap, add = T)
  image(MCraster, col=viridis(50), zlim=c(TSSmin,TSSmax), 
        axes=F, xlab="", ylab="", add=T)
  text(x=-75.8, y=38.75, labels=MClabel, cex=1.2, col="white", srt=90)
}
################################################################################
for(m in 1:4){
  myfile <- paste0("D:/Publication/TSS_",pattern[m],"_mc.tiff")
  tiff(file = myfile, width =7.5, height = 5.58, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  par(omi=c(0.25, 0.25, 0.05, 0.05)) 
  # omi is outer margin of the entire in inches 
  layout(matrix(c(1:12,rep(13, 6)),nrow=3,byrow=T), 
         heights = c(1,1, 0.2)) # Height of rows (relative) 
  ##Find files and plot them
  for(i in 1:12){
    pred <- raster(paste0("D:/Publication/TSS_",mon[i],".grd"), band = m)
    plotMC(pred, mymonths[i], TSSmin, TSSmax)
  }
  ##Color scale
  par(mai=c(0.25, 2, 0.02, 2))  #bottom, left, top, right
  barplot(rep(1,40), col=viridis(40), border=viridis(40), space=0, axes=F)
  axis(side=1, at = c(0, 10, 20, 30, 40), 
       labels=c("0", "7.5", "15", "22.5", "30"))
  mtext(side=1, "TSS (mg/L)", line=3, cex=1.1)
  dev.off()
}
