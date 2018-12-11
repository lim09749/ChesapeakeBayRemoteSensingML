##################################################################################################################################
## PROGRAM: CHL mc.R
## SPECIFICATIONS: Creates mc data (to speed up process)
#####################################################################################################################################

rm(list=ls())

require(rgdal)
require(raster)  # http://rpubs.com/etiennebr/visualraster
require(viridis)
require(rgeos)
require(e1071)      ## support vector regression
require(kernlab)    ## relevance vector machines
require(neuralnet)  ## neural networks

################################################################################
# Load 1 Year of Raster Data

setwd("D:/CB_MC_1km_from_L2") ##MC file location

CBnew <- read.csv("D:/Publication/CB_chl.csv")
model_var <- c("chl","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

##EDIT HERE
dir <- c("D:/Publication/CHL/RVM models/", "D:/Publication/CHL/SVM models/",
         "D:/Publication/CHL/GLM models/", "D:/Publication/CHL/ANN models/")
pattern <- c("rvm","svm", "glm", "ann")
wavelengths <- c("Rrs_412.tif", "Rrs_443.tif", "Rrs_469.tif", "Rrs_488.tif", "Rrs_531.tif",
                 "Rrs_547.tif", "Rrs_555.tif", "Rrs_645.tif", "Rrs_667.tif"," Rrs_678.tif")
mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
CHLmin = 0.01
CHLmax = 30

################################################################################
##Predict on mc dataset

##Acts as interface for different models
predictMC <- function(CBnew, CBmc, pattern, dir){
  if(pattern=="rvm" || pattern=="svm") return(vm_predictMC(CBnew, CBmc, pattern, dir))
  else if(pattern=="ann") ann_predictMC(CBnew, CBmc, pattern, dir)
  else if(pattern=="glm") glm_predictMC(CBnew, CBmc, pattern, dir)
}

##SVM and RVM
vm_predictMC <- function(CBnew, CBmc, pattern, dir){
  require(kernlab)
  require(e1071)
  ##Standardize data
  LogCBnew <- log10(CBnew)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBmc <- log10(CBmc)
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_center[2:11], FUN = "-")
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_scale[2:11], FUN = "/")
  
  ##Find the models
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBmc))
  
  ##Average their result
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern == "svm") {
      pred[j,] <- as.vector(10^(predict(chlsvm, StdLogCBmc)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(10^(predict(chlrvm, StdLogCBmc)*log_scale[1]+log_center[1]))
    }
    
  }
  pred
}


##GLM
glm_predictMC <- function(CBnew, CBmc, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBmc))
  
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <-  as.vector(predict(CHLglm,CBmc))
  }
  pred
}

##ANN
ann_predictMC <- function(CBnew, CBmc, pattern, dir){
  require(neuralnet)
  
  ##Standardize data
  LogCBnew <- log10(CBnew)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBmc <- log10(CBmc)
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_center[2:11], FUN = "-")
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_scale[2:11], FUN = "/")
  
  ##Iterate through different ann
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBmc))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(10^(compute(chlnet,
                                                 StdLogCBmc)
                                         $net.result*log_scale[1]+log_center[1]))
  }
  pred
}  


##Most common (helps us remove later)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.CHLmax(tabulate(match(v, uniqv)))]
}

##Get raster predictions
getRaster <- function(r, CBnew, df, pattern, dir, CHLmin, CHLmax){
  pmc <- rep(NA, nrow(df))
  pmc[!is.na(df$Rrs_412)] <- colMeans(predictMC(CBnew, na.omit(df), 
                                                pattern, dir))
  
  pred = raster(subset(r,1))
  values(pred) <- pmc
  values(pred)[values(pred)==getmode(values(pred))] <- NA
  values(pred)[values(pred) < CHLmin] <- CHLmin
  values(pred)[values(pred) > CHLmax] <- CHLmax
  pred
}
#####################################################################################
##Predictions
#Iterate through months
for(i in 1:12){
  myfiles = list.files(pattern = mon[i])
  
  r = stack(myfiles) # Load raster files into a stack
  # Need to set the coordinate system of the rasters
  crs(r)    = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  # Crop image to zoom in on the Bay
  r = crop(r, extent(-77, -75.6, 36.8, 39.6)) 
  
  ##Get rid of negative values
  df <- as.data.frame(values(r))
  colnames(df) = colnames(CBnew)[2:11]
  df[df <= 0] <- NA 
  
  ##Get us some values for different machine learning algorithms
  rrvm <- getRaster(r, CBnew, df, pattern[1], dir[1], CHLmin, CHLmax)
  rsvm <- getRaster(r, CBnew, df, pattern[2], dir[2], CHLmin, CHLmax)
  rglm <- getRaster(r, CBnew, df, pattern[3], dir[3], CHLmin, CHLmax)*10.41-12.82
  values(rglm)[values(rglm) < CHLmin] = CHLmin
  values(rglm)[values(rglm) > CHLmax] = CHLmax
  rann <- getRaster(r, CBnew, df, pattern[4], dir[4], CHLmin, CHLmax)*0.01678+9.85917
  values(rann)[values(rann) < CHLmin] = CHLmin
  values(rann)[values(rann) > CHLmax] = CHLmax
  
  ##stack into one raster
  pred <- stack(rrvm, rsvm, rglm, rann)
  names(pred) <- c("rvm", "svm", "glm", "ann")
  
  ##save it
  writeRaster(pred, filename = paste0("D:/Publication/CHL_",mon[i],".grd"), 
              overwrite = T)
}
