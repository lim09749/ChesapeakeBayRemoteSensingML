##################################################################################################################################
## PROGRAM: CHL scatterplots and metric.R
## SPECIFICATIONS: creates data for scatterplots
#####################################################################################################################################

rm(list=ls())
#####################################################################################################################################
## Libraries
require(e1071)      ## support vector regression
require(kernlab)    ## relevance vector machines
require(neuralnet)  ## neural networks

#####################################################################################################################################
##Setting up programs
setwd("D:/Publication")

##EDIT HERE
dir <- c("D:/Publication/CHL/RVM models/", "D:/Publication/CHL/SVM models/",
         "D:/Publication/CHL/GLM models/", "D:/Publication/CHL/ANN models/")
pattern <- c("rvm","svm", "glm", "ann")

##Model variables
model_var <- c("chl","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")

##Load test data
CBtest <- read.csv("chl_test.csv")
CBtest <- CBtest[,colnames(CBtest)%in%model_var]

##load entire dataset
CBnew <- read.csv("CB_chl.csv")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

#####################################################################################################################################
## Basic Functions

##Prediction metrics
mae <- function(error){return(mean(sqrt(error^2),na.rm=T))} 
mape <- function(error, values){
  return(mean(abs((error)/values) * 100, na.rm=T))
}
rmse <- function(error){return(sqrt(mean(error^2,na.rm=T)))} #Squares the error, finds the mean, 
r2  <- function(x, y){return(summary(lm(y~x))$r.squared)}

##log10 metrics, prints metrics together
predictMetrics <- function(pred, observed, modeltype){
  error <- log10(observed)-log10(pred)
  print(c(modeltype, "r2", "rmse", "mae", "mape", "mean"))
  c(r2(log10(pred), log10(observed)),
    rmse(error), mae(error), mape(error, log10(observed)),
    mean(error, na.rm=T))
}
#####################################################################################################################################
##Machine Learning Prediction Functions

##Predict on test dataset
##Acts as interface for different models
predictTest <- function(CBnew, CBtest, pattern, dir){
  if(pattern=="rvm" || pattern=="svm") return(vm_predictTest(CBnew, CBtest, pattern, dir))
  else if(pattern=="ann") ann_predictTest(CBnew, CBtest, pattern, dir)
  else if(pattern=="glm") glm_predictTest(CBnew, CBtest, pattern, dir)
}

##SVM and RVM
vm_predictTest <- function(CBnew, CBtest, pattern, dir){
  require(kernlab)
  require(e1071)
  
  ##Standardize data
  LogCBnew <- log10(CBnew)
  LogCBnew$chl <- log10(CBnew$chl)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBtest <- log10(CBtest)
  StdLogCBtest$chl <- log10(CBtest$chl)
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_center, FUN = "-")
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_scale, FUN = "/")
  
  ##Find the models
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBtest))
  
  ##Average their result
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern == "svm") {
      pred[j,] <- as.vector(10^(predict(chlsvm, StdLogCBtest)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(10^(predict(chlrvm, StdLogCBtest)*log_scale[1]+log_center[1]))
    }
    
  }
  pred
}

##ANN
ann_predictTest <- function(CBnew, CBtest, pattern, dir){
  require(neuralnet)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$chl <- log10(CBnew$chl)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBtest <- log10(CBtest)
  StdLogCBtest$chl <- log10(CBtest$chl)
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_center, FUN = "-")
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_scale, FUN = "/")
  
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBtest))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(10^(compute(chlnet,
                                                 StdLogCBtest[,-which(colnames(CBtest)=="chl")])
                                         $net.result*log_scale[1]+log_center[1]))
  }
  pred
}  

##GLM
glm_predictTest <- function(CBnew, CBtest, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBtest))
  
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <-  as.vector(predict(CHLglm,CBtest))
  }
  pred
}

##Predict on entire dataset
##Acts as interface for different models
predictTotal <- function(CBnew, pattern, dir){
  if(pattern=="rvm" || pattern=="svm") vm_predictTotal(CBnew, pattern, dir)
  else if(pattern=="ann") ann_predictTotal(CBnew, pattern, dir)
  else if(pattern=="glm") glm_predictTotal(CBnew, pattern, dir)
}

##SVM and RVM
vm_predictTotal <- function(CBnew, pattern, dir){
  require(kernlab)
  require(e1071)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$chl <- log10(CBnew$chl)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBnew))
  
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern == "svm") {
      pred[j,] <- as.vector(10^(predict(chlsvm, StdLogCBnew)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(10^(predict(chlrvm, StdLogCBnew)*log_scale[1]+log_center[1]))
    }
    
  }
  pred
}

##ANN
ann_predictTotal <- function(CBnew, pattern, dir){
  require(kernlab)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$chl <- log10(CBnew$chl)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBnew))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(10^(compute(chlnet,
                                                 StdLogCBnew[,-which(colnames(CBnew)=="chl")])
                                         $net.result*log_scale[1]+log_center[1]))
  }
  pred
}  

##GLM
glm_predictTotal <- function(CBnew, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBnew))
  
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(predict(CHLglm,CBnew))
  }
  pred
}


##Sparsity for RVM and SVM
vm_getSparsity <- function(CBnew, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  nvec <- rep(NA, length(files))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern=="svm") nvec[j] = chlsvm$tot.nSV
    else nvec[j] = chlrvm@nRV
  } 
  mean(nvec)
}

##Sparsity for GLM
glm_getSparsity <- function(CBnew, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  aic <- rep(NA, length(files))
  bic <- rep(NA, length(files))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    aic[j] <- AIC(CHLglm)
    bic[j] <- BIC(CHLglm)
  } 
  c(mean(aic), mean(bic))
}

##ANN
ann_getSparsity <- function(CBnew, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  aic <- rep(NA, length(files))
  bic <- rep(NA, length(files))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    aic[j] <- mean(chlnet$result.matrix["aic",])
    bic[j] <- mean(chlnet$result.matrix["bic",])
  } 
  c(mean(aic), mean(bic))
}
#####################################################################################################################################
##NASA Model

NASA_OCM3 <- function(CBnew){
  a0 <- rep(0.2424, nrow(CBnew))
  a1 <- rep(-2.7423, nrow(CBnew))
  a2 <- rep(1.8017, nrow(CBnew))
  a3 <- rep(0.0015, nrow(CBnew))
  a4 <- rep(-1.228, nrow(CBnew))
  blue = pmax(CBnew$Rrs_443, CBnew$Rrs_488)
  green = CBnew$Rrs_547
  
  predNASA = 10^(a0 + a1*(log10(blue/green))+a2*(log10(blue/green))^2+
                   a3*(log10(blue/green))^3+a4*(log10(blue/green))^4)
  return(predNASA)
}

#####################################################################################################################################
##Predict on test dataset

##This contains predictions (line below)
pred_test <- matrix(nrow = length(pattern)+1, ncol = nrow(CBtest))

##Fill in predictions
for(i in 1:length(pattern)){
  pred_test[i,] <- colMeans(predictTest(CBnew, CBtest, pattern = pattern[i], dir[i]))
}
pred_test[length(pattern)+1,] <- NASA_OCM3(CBtest)##Nasa algorithm

##Predict on entire dataset
pred_all <- matrix(nrow = length(pattern)+1, ncol = nrow(CBnew))

for(i in 1:length(pattern)){
  pred_all[i,] <- colMeans(predictTotal(CBnew, pattern = pattern[i], dir[i]))
}
pred_all[length(pattern)+1,] <- NASA_OCM3(CBnew)

##Send rownames 
rownames(pred_test) <- c(pattern, "NASA OCM3")
rownames(pred_all) <- c(pattern, "NASA OCM3")

##Get sparsity
#% Vec for rvm, svr
vm_sparsity <- c(vm_getSparsity(CBnew, pattern = pattern[1], dir[1]),
                  vm_getSparsity(CBnew, pattern = pattern[2], dir[2]))/(nrow(CBnew)-nrow(CBtest)-1)

##AIC and BIC for glm, ann
other_sparsity <- c(glm_getSparsity(CBnew, pattern = pattern[3], dir[3]),
                 ann_getSparsity(CBnew, pattern = pattern[4], dir[4]))

#####################################################################################################################################
##Metrics

##Test dataset metrics (Already log 10)
##Log 10 already applied in function
predictMetrics(pred_test[1,], CBtest$chl, pattern[1])
vm_sparsity[1]

predictMetrics(pred_test[2,], CBtest$chl, pattern[2])
vm_sparsity[2]

predictMetrics(pred_test[3,], CBtest$chl, pattern[3])
other_sparsity[1:2]

predictMetrics(pred_test[4,], CBtest$chl, pattern[4])
other_sparsity[3:4]

predictMetrics(pred_test[5,], CBtest$chl, "NASA OCM3")

mat <- rbind(c(predictMetrics(pred_test[1,], CBtest$chl, pattern[1]), vm_sparsity[1], NA, NA),
             c(predictMetrics(pred_test[2,], CBtest$chl, pattern[2]), vm_sparsity[2], NA, NA),
             c(predictMetrics(pred_test[3,], CBtest$chl, pattern[3]), NA, other_sparsity[1:2]),
             c(predictMetrics(pred_test[4,], CBtest$chl, pattern[4]), NA, other_sparsity[3:4]),
             c(predictMetrics(pred_test[5,], CBtest$chl, "NASA OCM3"), NA, NA, NA))

colnames(mat) <- c("r2", "rmse", "mae", "mape", "mean", "Sparsity ([%]) for RVM, SVM", "AIC for GLM, ANN", "BIC for GLM, ANN")
rownames(mat) <- c(pattern, "NASA OCM3")

##Save metrics and predictions
write.csv(mat, "chl_model_metrics.csv")
write.csv(rbind(pred_test, CBtest$chl), "chl_predict_test.csv")
write.csv(rbind(pred_all, CBnew$chl), "chl_predict_all.csv")
#####################################################################################################################################
##Scatterplot functions


plotscatter <- function(x, y, model, r2, lab){
  ##Set margin
  par(mai=c(0.1,0.1,0.1,0.1)) 
  ##plot it without axes, text, and then add them 
  plot(x, y, xlim = c(-1, 1.8), ylim = c(-1, 1.8), xaxt = "n", yaxt = "n", axes = FALSE)
  text(x  = -0.5, y = 1.5,labels=substitute(model))
  text(x  = -0.5, y = 1.3,labels=substitute(paste("r"[2]," = ", r2)))
  abline(a=0, b=1)
  axis(1, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[1])
  axis(2, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[2])
}


plotscatter_col <- function(x, y, model, r2, lab, col){
  ##Set margin
  par(mai=c(0.1,0.1,0.1,0.1)) 
  ##plot it without axes, text, and then add them 
  plot(x, y, xlim = c(-1, 1.8), ylim = c(-1, 1.8), xaxt = "n", yaxt = "n", col = col, axes = FALSE)
  text(x  = -0.5, y = 1.5,labels=substitute(model))
  text(x  = -0.5, y = 1.3,labels=substitute(paste("r"[2]," = ", r2)))
  abline(a=0, b=1)
  axis(1, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[1])
  axis(2, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[2])
}
#####################################################################################################################################
