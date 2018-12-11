##################################################################################################################################
## PROGRAM: CHL GLM cv with parallel.R
## SPECIFICATIONS: Develops ensemble generalized linear model for predicting chlorophyll
#####################################################################################################################################

rm(list=ls())
#####################################################################################################################################
## Libraries
##glm part of basic package in R
require(doParallel) ## parallel processing library

#####################################################################################################################################
##Load in data

##Set up directory
##Necessary files: CB_chl.csv
setwd("D:/Publication")

##read in file
CBnew <- read.csv("CB_chl.csv")

##Keep important information (chl and Rrs)
model_var <- c("chl","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]
######################################################################################################################
##Some parameters

##percent of tail to be trimmed when averaging predictions
##Not effective here, but may be effective in other settings
t <- 0

##percent as testing set
p <- 0.2

######################################################################################################################
## Functions

#Mean absolute error
mae <- function(error){return(mean(sqrt(error^2),na.rm=T))} 

#Mean absolute percentage error
mape <- function(error, values){
  return(mean(abs((error)/values) * 100, na.rm=T))
}

#Root mean square error
rmse <- function(error){return(sqrt(mean(error^2,na.rm=T)))} #Squares the error, finds the mean, 

#give r-squared
r2  <- function(x, y){return(summary(lm(y~x))$r.squared)}

#error matrices for ease of use
#predmat- predmat[i,j]- the jth prediction of the ith model
#val- observed in-situ values
create_error_mat <- function(predmat, val){
  errormat <- matrix(NA, nrow=nrow(predmat), ncol=length(val))
  for(i in 1:nrow(predmat)){
    errormat[i,] <- predmat[i,]-val
  }
  errormat
}

#Get vector of different r-squared values for different models
#predmat- predmat[i,j]- the jth prediction of the ith model
#val- observed in-situ values
r2_vec <- function(predmat, val){
  test_r2 <-  rep(NA, nrow(predmat))
  for(i in 1:nrow(predmat)){
    test_r2[i] = r2(predmat[i,],val)
  }
  test_r2
}

#Modified colmeans to account for trim
#In this program, acts as colMeans (because trim=0) but may be useful 
colMeans_trim <- function(predmat, trim){
  apply(predmat, MARGIN = 2, FUN= function(x) mean(x,trim=trim))
}
######################################################################################################################
##Data randomization and standardization

##Randomize data
CBnew<-CBnew[sample(nrow(CBnew)),]

##Standardization of data (center around 0 and set standard deviation to 1)
StdCBnew <- scale(CBnew)
no_log_center <- attr(StdCBnew,"scaled:center")
no_log_scale <- attr(StdCBnew,"scaled:scale")
attr(StdCBnew, "scaled:center") <- NULL
attr(StdCBnew, "scaled:scale") <- NULL
StdCBnew <- as.data.frame(StdCBnew) ##Only standardized

LogCBnew <- log10(CBnew)
StdLogCBnew <- scale(LogCBnew)
log_center <- attr(StdLogCBnew,"scaled:center")
log_scale <- attr(StdLogCBnew,"scaled:scale")
attr(StdLogCBnew, "scaled:center") <- NULL
attr(StdLogCBnew, "scaled:scale") <- NULL
StdLogCBnew <- as.data.frame(StdLogCBnew) ##Standardized and set standard deviation to 1
######################################################################################################################
##Seperate testing and training dataset

testIndexes <- sample.int(n = nrow(CBnew), size = floor(p*nrow(CBnew)), replace = F)

##Testing dataset
CBtest <- CBnew[testIndexes, ]
StdCBtest <- StdCBnew[testIndexes,]
LogCBtest <- LogCBnew[testIndexes, ]
StdLogCBtest <- StdLogCBnew[testIndexes, ]

##Training dataset
modCBnew <- CBnew[-testIndexes,]
modStdCBnew <- StdCBnew[-testIndexes,]
modLogCBnew <- LogCBnew[-testIndexes,]
modStdLogCBnew <- StdLogCBnew[-testIndexes,]
######################################################################################################################
##Train different models on training dataset

##Create cluster for parallel processing
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
##Train different models on training dataset

## Used parallel processing to make it faster 
perf <- foreach(i=1:nrow(modCBnew), .combine=rbind) %dopar%{
  valIndexes <- c(i)
  
  ##Create training and validation datasets
  CBtrain <- modCBnew[-c(i), ]
  StdCBtrain <- modStdCBnew[-c(i),]
  LogCBtrain <- modLogCBnew[-c(i), ]
  StdLogCBtrain <- modStdLogCBnew[-c(i), ]
  
  ##Train glm
  CHLglm<-glm(chl ~ ., data=CBtrain,family=Gamma(link="log"))
  
  ##Save the model 
  save(CHLglm, file = paste0("D:/Publication/CHL/GLM models/glm",i,".RData"))
  
  ##Adjoin the information given
  c(AIC(CHLglm), ##Akaike Information Criterion
    BIC(CHLglm), ##Bayesian Information Criterion
    as.vector(predict(CHLglm,CBtest))) ##predicted on test dataset
}
stopImplicitCluster()
##NOTE: ith row of perf IS NOT the performance ith model
##      sometimes the i+1th model finishes training faster than the ith model
predmat <- as.matrix(perf[,-c(1,2)])
errormat <- create_error_mat(predmat, CBtest$chl)
pred <- colMeans_trim(predmat, t)

##Some basic metrics
aic <- as.vector(perf[,1])
bic <- as.vector(perf[,2])
test_r2 <- r2_vec(log10(predmat), log10(CBtest$chl))

##Plot results
plot(log10(pred),log10(CBtest$chl))
r2(log10(pred),log10(CBtest$chl))

##Save training and testing dataset
write.csv(CBtest, "CHL_test.csv")
