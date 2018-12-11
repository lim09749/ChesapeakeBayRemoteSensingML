##################################################################################################################################
## PROGRAM: TSS SVM cv with parallel.R
## SPECIFICATIONS: Develops ensemble SVM
## for predicting TSS
#####################################################################################################################################

rm(list=ls())
#####################################################################################################################################
## Libraries
require(e1071)      ## support vector regression
require(doParallel) ## parallel processing library

#####################################################################################################################################
##Load in data

##Set up directory
##Necessary files: CB_tss.csv
setwd("D:/Publication")

##read in file
CBnew <- read.csv("CB_tss.csv")

##Keep important information (TSS and Rrs)
model_var <- c("TSS","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
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

##Standardization
StdCBnew <- scale(CBnew)
no_log_center <- attr(StdCBnew,"scaled:center")
no_log_scale <- attr(StdCBnew,"scaled:scale")
attr(StdCBnew, "scaled:center") <- NULL
attr(StdCBnew, "scaled:scale") <- NULL
StdCBnew <- as.data.frame(StdCBnew)

LogCBnew <- log10(CBnew)
LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
StdLogCBnew <- scale(LogCBnew)
log_center <- attr(StdLogCBnew,"scaled:center")
log_scale <- attr(StdLogCBnew,"scaled:scale")
attr(StdLogCBnew, "scaled:center") <- NULL
attr(StdLogCBnew, "scaled:scale") <- NULL
StdLogCBnew <- as.data.frame(StdLogCBnew)
######################################################################################################################
##Seperate testing and training dataset

testIndexes <- sample.int(n = nrow(CBnew), size = floor(p*nrow(CBnew)), replace = F)

#testIndexes <- which(do.call(paste0, CBnew) %in% do.call(paste0, z))
#^This was used to make sure the testing and training sets for the 
#different types of models are the same

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

## Used parallel processing to make it faster
perf <- foreach(i=1:nrow(modCBnew), .combine=rbind) %dopar%{
  require(e1071) ##Needs to load libraries in R
  valIndexes <- c(i)

  ##Create training and validation datasets
  CBtrain <- modCBnew[-valIndexes, ]
  StdCBtrain <- modStdCBnew[-valIndexes,]
  LogCBtrain <- modLogCBnew[-valIndexes, ]
  StdLogCBtrain <- modStdLogCBnew[-valIndexes, ]
  
  ##Train svm
  tuneResult <- tune(svm,TSS~ .,
                     data=StdLogCBtrain,
                     ranges = list(epsilon = seq(0.1,0.5,0.1),
                                   cost = 2^(0:3)),
                     type="eps-regression",kernel="linear")
  TSSsvm <- tuneResult$best.model
  
  
  ##Save the model 
  save(TSSsvm, file = paste0("D:/Publication/TSS/SVM models/svm",i,".RData"))
  
  ##Adjoin the information given
  c(TSSsvm$tot.nSV, ##Number of support vectors
    as.vector(inverse_square(
      predict(TSSsvm, StdLogCBtest)*log_scale[1]+log_center[1]))) ##predicted on test dataset
}
stopImplicitCluster()
##NOTE: ith row of perf IS NOT the performance ith model
##REASON: sometimes the i+1th model finishes training faster than the ith model
predmat <- as.matrix(perf[,-c(1)])
errormat <- create_error_mat(predmat, CBtest$TSS)
pred <- colMeans_trim(predmat, t)

##Some basic metrics
n_svr <- as.vector(perf[,1])
test_r2 <- r2_vec(log10(predmat), log10(CBtest$TSS))

##Plot results
plot(log10(pred),log10(CBtest$TSS))
r2(log10(pred),log10(CBtest$TSS))

##Save training and testing dataset
write.csv(CBtest, "TSS_test.csv")
