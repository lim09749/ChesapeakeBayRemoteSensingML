##################################################################################################################################
## PROGRAM: TSS line plots mc.R
## SPECIFICATIONS: Creates line plots monthly
#####################################################################################################################################

rm(list=ls())
#####################################################################################################################################
## Libraries
require(e1071)      ## support vector regression
require(kernlab)    ## relevance vector machines
require(neuralnet)  ## neural networks
require(raster)
require(doParallel)
#####################################################################################################################################
##Setting up programs
setwd("D:/CB_MC_1km_from_L2")

##Parameters
dir <- c("D:/Publication/TSS/RVM models/", "D:/Publication/TSS/SVM models/",
         "D:/Publication/TSS/GLM models/", "D:/Publication/TSS/ANN models/")
pattern <- c("rvm","svm", "glm", "ann")

# For standardization metrics
model_var <- c("TSS","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")

CBnew <- read.csv("D:/Publication/CB_tss.csv")
datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CBnew$DateTime),split = " "))
                                     [seq(1, nrow(CBnew)*2,2)],format = "%m/%d/%Y"), 6, 7))
##In-situ vector averages
real_upper <- rep(NA, 12)
real_middle <- rep(NA, 12)
real_lower <- rep(NA, 12)

##Read it in
for(i in 1:12){
  real_upper[i] = mean(CBnew$TSS[datenum==i & CBnew$lat >= 38.75])
  real_middle[i] = mean(CBnew$TSS[datenum==i & CBnew$lat < 38.75 & CBnew$lat >= 37.75])
  real_lower[i] = mean(CBnew$TSS[datenum==i & CBnew$lat < 37.75])
}  

##Frequencies (make sure enough values in each part in each month)
for(i in 1:12){
  print(c(length(which(datenum==i & CBnew$lat >= 38.75)), 
          length(which(datenum==i & CBnew$lat < 38.75 & CBnew$lat >= 37.75)),
          length(which(datenum==i & CBnew$lat < 37.75))
          ))
}

##Save important values
CBnew <- CBnew[,colnames(CBnew)%in%model_var]
####################################################################################################################################
## Functions for predicting

##for TSS distribution
inverse_sqrt <- function(x){1/sqrt(x)}
inverse_square <- function(x){1/x^2}

##Machine Learning Prediction Functions

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

#########################################################################################################################
##Predictions on mc
avg_mc <- function(CBnew, pattern, dir, mon, region){
  require(doParallel)
  avg <- matrix()
  avg <- foreach(i=1:length(mon), .combine=rbind) %dopar%{
    require(e1071)      ## support vector regression
    require(kernlab)    ## relevance vector machines
    require(neuralnet)  ## neural networks
    require(raster)
    myfiles = list.files(pattern = mon[i])
    
    r = stack(myfiles) # Load raster files into a stack
    
    ##mc rrs and corresponding latitutde
    CBmc <- as.data.frame(values(r))
    lat <- coordinates(r)[, 2]
    
    ##for TSS distribution
    inverse_sqrt <- function(x){1/sqrt(x)}
    inverse_square <- function(x){1/x^2}
    
    ##Machine Learning Prediction Functions
    
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
    
    ##Clean data
    lat <- lat[which(!is.na(rowSums(CBmc)) & apply(CBmc, 1, function(row) all(row   > 0 )))]
    CBmc <- CBmc[which(!is.na(rowSums(CBmc)) & apply(CBmc, 1, function(row) all(row   > 0 ))),]
    if(region == "Upper"){
      CBmc <- CBmc[lat >= 38.75,]
    }else if(region == "Middle"){
      CBmc <- CBmc[lat < 38.75 & lat >= 37.75,]
    }else if(region == "Lower"){
      CBmc <- CBmc[lat < 37.75,]
    }
    colnames(CBmc) = c("Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
                       "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")
    c(i,mean(predictMC(CBnew, CBmc, pattern = pattern[1], dir[1])),
      mean(predictMC(CBnew, CBmc, pattern = pattern[2], dir[2])),
      mean(predictMC(CBnew, CBmc, pattern = pattern[3], dir[3])),
      mean(predictMC(CBnew, CBmc, pattern = pattern[4], dir[4])))
  }
  avg
}

##Get dates
mc_dates <- function(filename){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CB$DateTime),split = " "))[seq(1, length(CB$DateTime)*2,2)], 
                                       format = "%m/%d/%Y"), 6, 7))
  datenum
}

##In-situ TSS for entire dataset (No need to have it cleaned)
mc_real <- function(filename, type){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- mc_dates(filename)
  
  mc_avg = rep(NA, 12)
  for(i in 1:12){
    if(type == "Upper"){
      mc_avg[i] = mean(CB$TSS[datenum==i & CB$lat >= 38.75])
    }
    else if(type == "Middle"){
      mc_avg[i] = mean(CB$TSS[datenum==i & CB$lat < 38.75 & CB$lat >= 37.75])
    }else if(type == "Lower"){
      mc_avg[i] = mean(CB$TSS[datenum==i & CB$lat < 37.75])
    }
    else{
      mc_avg[i] = mean(CB$TSS[datenum==i])
    }
  }
  mc_avg
}
#########################################################################################################################
##Get values
filename = "D:/Publication/TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt"
complete_real_upper <- mc_real(filename, type = "Upper")
complete_real_middle <- mc_real(filename, type = "Middle")
complete_real_lower <- mc_real(filename, type = "Lower")

mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

##Upper (Oligohaline)
upper_avg <- avg_mc(CBnew, pattern, dir, mon, region = "Upper") 
upper_avg <- upper_avg[order(upper_avg[,1]),]
upper_avg <- t(upper_avg)
upper_avg <- upper_avg[2:5,]
rownames(upper_avg) = c("RVM", "SVR", "GLM", "ANN")
colnames(upper_avg) = mon

##Middle (Mesohaline)
middle_avg <- avg_mc(CBnew, pattern, dir, mon, region = "Middle") 
middle_avg <- middle_avg[order(middle_avg[,1]),]
middle_avg <- t(middle_avg)
middle_avg <- middle_avg[2:5,]

rownames(middle_avg) = c("RVM", "SVR", "GLM", "ANN")
colnames(middle_avg) = mon

##Lower(Polyhaline)
lower_avg <- avg_mc(CBnew, pattern, dir, mon, region = "Lower") 
lower_avg <- lower_avg[order(lower_avg[,1]),]
lower_avg <- t(lower_avg)
lower_avg <- lower_avg[2:5,]

rownames(lower_avg) = c("RVM", "SVR", "GLM", "ANN")
colnames(lower_avg) = mon

stopCluster(cl)

##Write data
write.csv(upper_avg, file = "D:/Publication/TSS_upper_avg.csv")
write.csv(middle_avg, file = "D:/Publication/TSS_middle_avg.csv")
write.csv(lower_avg, file = "D:/Publication/TSS_lower_avg.csv")

##########################################################################################################################
