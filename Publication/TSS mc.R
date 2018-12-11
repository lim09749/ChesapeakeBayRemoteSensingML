##################################################################################################################################
## PROGRAM: TSS mc.R
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

CBnew <- read.csv("D:/Publication/CB_tss.csv")
model_var <- c("TSS","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

##directory, some parameters
dir <- c("D:/Publication/TSS/RVM models/", "D:/Publication/TSS/SVM models/",
         "D:/Publication/TSS/GLM models/", "D:/Publication/TSS/ANN models/")
pattern <- c("rvm","svm", "glm", "ann")
wavelengths <- c("Rrs_412.tif", "Rrs_443.tif", "Rrs_469.tif", "Rrs_488.tif", "Rrs_531.tif",
                 "Rrs_547.tif", "Rrs_555.tif", "Rrs_645.tif", "Rrs_667.tif"," Rrs_678.tif")
mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
TSSmin = 0.01 ##cut off for cleaner plots
TSSmax = 30

################################################################################
##Predict on mc dataset

##for TSS distribution
inverse_sqrt <- function(x){1/sqrt(x)}
inverse_square <- function(x){1/x^2}

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
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
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
      pred[j,] <- as.vector(inverse_square(predict(TSSsvm, StdLogCBmc)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(inverse_square(predict(TSSrvm, StdLogCBmc)*log_scale[1]+log_center[1]))
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
    pred[j,] <-  as.vector(predict(TSSglm,CBmc))
  }
  pred
}

##ANN
ann_predictMC <- function(CBnew, CBmc, pattern, dir){
  require(neuralnet)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBmc <- log10(CBmc)
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_center[2:11], FUN = "-")
  StdLogCBmc <- sweep(StdLogCBmc, 2L, log_scale[2:11], FUN = "/")
  
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBmc))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(inverse_square(compute(TSSnet,
                                                 StdLogCBmc)
                                         $net.result*log_scale[1]+log_center[1]))
  }
  pred
}  


##Most common
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.TSSmax(tabulate(match(v, uniqv)))]
}

getRaster <- function(r, CBnew, df, pattern, dir, TSSmin, TSSmax){
  pmc <- rep(NA, nrow(df))
  pmc[!is.na(df$Rrs_412)] <- colMeans(predictMC(CBnew, na.omit(df), 
                                                pattern, dir))
  
  pred = raster(subset(r,1))
  values(pred) <- pmc
  values(pred)[values(pred)==getmode(values(pred))] <- NA
  values(pred)[values(pred) < TSSmin] <- TSSmin
  values(pred)[values(pred) > TSSmax] <- TSSmax
  pred
}

################################################################################
##Predictions
for(i in 1:12){
  myfiles = list.files(pattern = mon[i])
  
  r = stack(myfiles) # Load raster files into a stack
  # Need to set the coordinate system of the rasters
  crs(r)    = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  # Crop image to zoom in on the Bay
  r = crop(r, extent(-77, -75.6, 36.8, 39.6)) 

  df <- as.data.frame(values(r))
  colnames(df) = colnames(CBnew)[2:11]
  df[df <= 0] <- NA 
  
  ##Get predictions for each ml algorithm
  rrvm <- getRaster(r, CBnew, df, pattern[1], dir[1], TSSmin, TSSmax)
  rsvm <- getRaster(r, CBnew, df, pattern[2], dir[2], TSSmin, TSSmax)
  rglm <- getRaster(r, CBnew, df, pattern[3], dir[3], TSSmin, TSSmax)*10.26-12.99
  rann <- getRaster(r, CBnew, df, pattern[4], dir[4], TSSmin, TSSmax)*-5.163e-06+8.591e+00
  values(rglm)[values(rglm) < TSSmin] = TSSmin
  values(rglm)[values(rglm) > TSSmax] = TSSmax
  values(rann)[values(rann) < TSSmin] = TSSmin
  values(rann)[values(rann) > TSSmax] = TSSmax
  
  pred <- stack(rrvm, rsvm, rglm, rann)
  names(pred) <- c("rvm", "svm", "glm", "ann")
  writeRaster(pred, filename = paste0("D:/Publication/TSS_",mon[i],".grd"), 
              overwrite = T)
}
