##################################################################################################################################
## PROGRAM: TSS scatterplots and metric.R
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
dir <- c("D:/Publication/TSS/RVM models/", "D:/Publication/TSS/SVM models/",
         "D:/Publication/TSS/GLM models/", "D:/Publication/TSS/ANN models/")
pattern <- c("rvm","svm", "glm", "ann")

##Model variables
model_var <- c("TSS","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")

##Load test data
CBtest <- read.csv("TSS_test.csv")
CBtest <- CBtest[,colnames(CBtest)%in%model_var]

##load entire dataset
CBnew <- read.csv("CB_tss.csv")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

#####################################################################################################################################
## Basic Functions

##for TSS distribution
inverse_sqrt <- function(x){1/sqrt(x)}
inverse_square <- function(x){1/x^2}

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
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBtest <- log10(CBtest)
  StdLogCBtest$TSS <- inverse_sqrt(CBtest$TSS)
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_center, FUN = "-")
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_scale, FUN = "/")
  
  ##Find the models
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBtest))
  
  ##Average their result
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern == "svm") {
      pred[j,] <- as.vector(inverse_square(predict(TSSsvm, StdLogCBtest)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(inverse_square(predict(TSSrvm, StdLogCBtest)*log_scale[1]+log_center[1]))
    }
    
  }
  pred
}

##ANN
ann_predictTest <- function(CBnew, CBtest, pattern, dir){
  require(neuralnet)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
  StdLogCBnew <- scale(LogCBnew)
  log_center <- attr(StdLogCBnew,"scaled:center")
  log_scale <- attr(StdLogCBnew,"scaled:scale")
  attr(StdLogCBnew, "scaled:center") <- NULL
  attr(StdLogCBnew, "scaled:scale") <- NULL
  StdLogCBnew <- as.data.frame(StdLogCBnew)
  
  StdLogCBtest <- log10(CBtest)
  StdLogCBtest$TSS <- inverse_sqrt(CBtest$TSS)
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_center, FUN = "-")
  StdLogCBtest <- sweep(StdLogCBtest, 2L, log_scale, FUN = "/")
  
  files <- list.files(path = dir, pattern=pattern)
  pred <- matrix(NA, nrow = length(files), ncol = nrow(CBtest))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    pred[j,] <- as.vector(inverse_square(compute(TSSnet,
                                                 StdLogCBtest[,-which(colnames(CBtest)=="TSS")])
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
    pred[j,] <-  as.vector(predict(TSSglm,CBtest))
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
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
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
      pred[j,] <- as.vector(inverse_square(predict(TSSsvm, StdLogCBnew)*log_scale[1]+log_center[1]))
    }
    else {
      pred[j,] <- as.vector(inverse_square(predict(TSSrvm, StdLogCBnew)*log_scale[1]+log_center[1]))
    }
      
  }
  pred
}

##ANN
ann_predictTotal <- function(CBnew, pattern, dir){
  require(kernlab)
  
  LogCBnew <- log10(CBnew)
  LogCBnew$TSS <- inverse_sqrt(CBnew$TSS)
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
    pred[j,] <- as.vector(inverse_square(compute(TSSnet,
                                                 StdLogCBnew[,-which(colnames(CBnew)=="TSS")])
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
    pred[j,] <- as.vector(predict(TSSglm,CBnew))
  }
  pred
}
#####################################################################################################################################
##Literature models (Ondrusek et al. 2012, Wang et al. 2009)

##Calculate normalized water leaving radiances
calc_nLw <- function(CBnew){
  fuller <-  	c(172.912, 187.622, 205.878, 194.933, 185.747, 
                186.539, 183.869, 157.811, 152.255, 148.052)
  
  rrs    <-   c("Rrs_412", "Rrs_443", "Rrs_469", "Rrs_488", "Rrs_531",
                "Rrs_547", "Rrs_555", "Rrs_645", "Rrs_667", "Rrs_678")
  for(i in 1:length(rrs)) CBnew[,colnames(CBnew)==rrs[i]] =CBnew[,colnames(CBnew)==rrs[i]] *fuller[i]
  return(CBnew) 
}

##Calculate R
calc_R <- function(CBnew){
  nlw_CBnew <- calc_nLw(CBnew)
  fuller <-  	c(172.912, 187.622, 205.878, 194.933, 185.747, 
                186.539, 183.869, 157.811, 152.255, 148.052)
  rrs    <-   c("Rrs_412", "Rrs_443", "Rrs_469", "Rrs_488", "Rrs_531",
                "Rrs_547", "Rrs_555", "Rrs_645", "Rrs_667", "Rrs_678")
  for(i in 1:length(rrs)){
    nlw_vec <- nlw_CBnew[,colnames(nlw_CBnew)==rrs[i]]
    nlw_CBnew[,colnames(nlw_CBnew)==rrs[i]] = 4*nlw_vec/(0.52*fuller[i]+1.7*nlw_vec)
  } 
  nlw_CBnew
}

##Ondrusek et al., 2012
Ondr_2012 <- function(CBnew){
  .CBnew <- calc_nLw(CBnew)
  3.8813*.CBnew$Rrs_645^3 - 13.822*.CBnew$Rrs_645^2 + 19.61*.CBnew$Rrs_645
}

##clear formula: Mueller, 2000
kd_clear <- function(CBnew){
  nlw_CB <- calc_nLw(CBnew)
  0.016+0.15645*(nlw_CB$Rrs_488/nlw_CB$Rrs_555)^(-1.5401) 
}

##Turbid formula: in wang paper
kd_turbid <- function(CBnew, eq){
  r_CBnew <- calc_R(CBnew)
  if(eq != 9 && eq != 12) print("no such equation exists")
  else if(eq == 9){
    2.697*10^(-4)/r_CBnew$Rrs_488+1.045*r_CBnew$Rrs_667/r_CBnew$Rrs_488+4.18*(7*10^(-4)+2.7135*r_CBnew$Rrs_667)*(1-0.52*exp(-2.533*10^(-3)/r_CBnew$Rrs_488-9.817*r_CBnew$Rrs_667/r_CBnew$Rrs_488))
  }
  else{
    -9.785*10^(-4)/r_CBnew$Rrs_488+0.8321*r_CBnew$Rrs_645/r_CBnew$Rrs_488+4.18*(-2.54*10^(-3)+2.1598*r_CBnew$Rrs_645)*(1-0.52*exp(9.19*10^(-3)/r_CBnew$Rrs_488-7.81*r_CBnew$Rrs_645/r_CBnew$Rrs_488))
  }
}

##Calc normalization for turbid and open ocean coefficients
calc_W <- function(CBnew){
  frac <- CBnew$Rrs_667/CBnew$Rrs_488
  w <- rep(NA, length(frac))
  w[frac < 0.2604] = 0
  w[frac>=0.2604 & frac <= 0.4821] = -1.175+4.512*frac[frac>=0.2604 & frac <= 0.4821]
  w[frac > 0.4821] = 1
  w
}

##Wang model
Wang_2009 <- function(CBnew, eq){
  w <- calc_W(CBnew)
  kd_comb = (1-w)*kd_clear(CBnew)+w*kd_turbid(CBnew, eq)
  kd_comb
}


##Sparsity for RVM and SVM
vm_getSparsity <- function(CBnew, pattern, dir){
  files <- list.files(path = dir, pattern=pattern)
  nvec <- rep(NA, length(files))
  for(j in 1:length(files)){
    load(paste0(dir,files[j]))
    if(pattern=="svm") nvec[j] = TSSsvm$tot.nSV
    else nvec[j] = TSSrvm@nRV
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
    aic[j] <- AIC(TSSglm)
    bic[j] <- BIC(TSSglm)
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
    aic[j] <- mean(TSSnet$result.matrix["aic",])
    bic[j] <- mean(TSSnet$result.matrix["bic",])
  } 
  c(mean(aic), mean(bic))
}
#####################################################################################################################################
##Predict on test dataset
b_9 <-  lm(CBnew$TSS~Wang_2009(CBnew, 9))$coefficients[1]
m_9 <-  lm(CBnew$TSS~Wang_2009(CBnew, 9))$coefficients[2]

b_12 <-  lm(CBnew$TSS~Wang_2009(CBnew, 12))$coefficients[1]
m_12 <-  lm(CBnew$TSS~Wang_2009(CBnew, 12))$coefficients[2]

pred_test <- matrix(nrow = length(pattern)+3, ncol = nrow(CBtest))

for(i in 1:length(pattern)){
  pred_test[i,] <- colMeans(predictTest(CBnew, CBtest, pattern = pattern[i], dir[i]))
}
pred_test[length(pattern)+1,] <- Ondr_2012(CBtest)
pred_test[length(pattern)+2,] <- Wang_2009(CBtest, 9)*m_9+b_9
pred_test[length(pattern)+3,] <- Wang_2009(CBtest, 12)*m_12+b_12


##Predict on entire dataset
pred_all <- matrix(nrow = length(pattern)+3, ncol = nrow(CBnew))

for(i in 1:length(pattern)){
  pred_all[i,] <- colMeans(predictTotal(CBnew, pattern = pattern[i], dir[i]))
}
pred_all[length(pattern)+1,] <- Ondr_2012(CBnew)
pred_all[length(pattern)+2,] <- Wang_2009(CBnew, 9)*m_9+b_9
pred_all[length(pattern)+3,] <- Wang_2009(CBnew, 12)*m_12+b_12

rownames(pred_test) <- c(pattern, "Ondrusek, 2012", "Wang et al., 2009 with Eq. 9", "Wang et al., 2009 with Eq. 12")
rownames(pred_all) <- c(pattern, "Ondrusek, 2012", "Wang et al., 2009 with Eq. 9", "Wang et al., 2009 with Eq. 12")

##Get sparsity
#% Vec for rvm, svr
vm_sparsity <- c(vm_getSparsity(CBnew, pattern = pattern[1], dir[1]),
                 vm_getSparsity(CBnew, pattern = pattern[2], dir[2]))/(nrow(CBnew)-nrow(CBtest)-1)

##AIC and BIC for glm, ann
other_sparsity <- c(glm_getSparsity(CBnew, pattern = pattern[3], dir[3]),
                    ann_getSparsity(CBnew, pattern = pattern[4], dir[4]))

#####################################################################################################################################
##Metrics

##Test dataset metrics
predictMetrics(pred_test[1,], CBtest$TSS, pattern[1])
vm_sparsity[1]
predictMetrics(pred_test[2,], CBtest$TSS, pattern[2])
vm_sparsity[2]
predictMetrics(pred_test[3,], CBtest$TSS, pattern[3])
other_sparsity[1:2]
predictMetrics(pred_test[4,], CBtest$TSS, pattern[4])
other_sparsity[3:4]
predictMetrics(pred_test[5,], CBtest$TSS, "Ondrusek, 2012")
predictMetrics(pred_test[6,], CBtest$TSS, "Wang et al., 2009 with Eq. 9") ##Best
predictMetrics(pred_test[7,], CBtest$TSS, "Wang et al., 2009 with Eq. 12")

mat <- rbind(c(predictMetrics(pred_test[1,], CBtest$TSS, pattern[1]), vm_sparsity[1], NA, NA),
             c(predictMetrics(pred_test[2,], CBtest$TSS, pattern[2]), vm_sparsity[2], NA, NA),
             c(predictMetrics(pred_test[3,], CBtest$TSS, pattern[3]), NA, other_sparsity[1:2]),
             c(predictMetrics(pred_test[4,], CBtest$TSS, pattern[4]), NA, other_sparsity[3:4]),
             c(predictMetrics(pred_test[5,], CBtest$TSS, "Ondrusek, 2012"), rep(NA, 3)),
             c(predictMetrics(pred_test[6,], CBtest$TSS, "Wang et al., 2009 with Eq. 9"), rep(NA, 3)), ##Best
             c(predictMetrics(pred_test[7,], CBtest$TSS, "Wang et al., 2009 with Eq. 12"), rep(NA, 3)))
colnames(mat) <- c("r2", "rmse", "mae", "mape", "mean", "Sparsity ([%]) for RVM, SVM", "AIC for GLM, ANN", "BIC for GLM, ANN")
rownames(mat) <- c(pattern, "Ondrusek, 2012", "Wang et al., 2009 with Eq. 9", "Wang et al., 2009 with Eq. 12")

write.csv(mat, "TSS_model_metrics.csv")
write.csv(rbind(pred_test, CBtest$TSS), "TSS_predict_test.csv")
write.csv(rbind(pred_all, CBnew$TSS), "TSS_predict_all.csv")
#####################################################################################################################################
##Scatterplot functions

plotscatter <- function(x, y, model, r2, lab){
  par(mai=c(0.1,0.1,0.1,0.1)) 
  plot(x, y, xlim = c(-1, 1.8), ylim = c(-1, 1.8), xaxt = "n", yaxt = "n", axes = FALSE)
  text(x  = -0.5, y = 1.5,labels=substitute(model))
  text(x  = -0.5, y = 1.3,labels=substitute(paste("r"[2]," = ", r2)))
  abline(a=0, b=1)
  axis(1, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[1])
  axis(2, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[2])
}


plotscatter_col <- function(x, y, model, r2, lab, col){
  par(mai=c(0.1,0.1,0.1,0.1)) 
  plot(x, y, xlim = c(-1, 1.8), ylim = c(-1, 1.8), xaxt = "n", yaxt = "n", col = col, axes = FALSE)
  text(x  = -0.5, y = 1.5,labels=substitute(model))
  text(x  = -0.5, y = 1.3,labels=substitute(paste("r"[2]," = ", r2)))
  abline(a=0, b=1)
  axis(1, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[1])
  axis(2, at =c(-1, -0.5, 0, 0.5, 1, 1.5), labels = lab[2])
}
#####################################################################################################################################
##Plotting 

##TSS
tiff(file = "TSS_test_only.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

#windows(width=7.48, height=8.60)
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
layout(matrix(1:6,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  
plotscatter(log10(pred_test[1,]), log10(CBtest$TSS), model = "RVM", r2 = "0.48", lab = c(FALSE, TRUE))
plotscatter(log10(pred_test[2,]), log10(CBtest$TSS), model = "SVR", r2 = "0.49", lab = c(FALSE, FALSE))
plotscatter(log10(pred_test[3,]), log10(CBtest$TSS), model = "GLM", r2 = "0.51", lab = c(FALSE, TRUE))
plotscatter(log10(pred_test[4,]), log10(CBtest$TSS), model = "ANN", r2 = "0.05", lab = c(FALSE, FALSE))
plotscatter(log10(pred_test[5,]), log10(CBtest$TSS), model = "Ondrusek et al., 2012", r2 = "0.37", 
            lab = c(TRUE, TRUE))
plotscatter(log10(pred_test[6,]), log10(CBtest$TSS), model = "Wang et al., 2009 with Eq. 9", r2 = "0.34", 
            lab = c(TRUE, FALSE))
mtext("Observed TSS", side = 1, outer= TRUE, line = 2)
mtext("Predicted TSS", side = 2, outer = TRUE, line = 2)

dev.off()

##Total
tiff(file = "TSS_total.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

#windows(width=7.48, height=8.60)
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
layout(matrix(1:6,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  
plotscatter(log10(pred_all[1,]), log10(CBnew$TSS), model = "RVM", r2 = "0.48", lab = c(FALSE, TRUE))
plotscatter(log10(pred_all[2,]), log10(CBnew$TSS), model = "SVR", r2 = "0.49", lab = c(FALSE, FALSE))
plotscatter(log10(pred_all[3,]), log10(CBnew$TSS), model = "GLM", r2 = "0.51", lab = c(FALSE, TRUE))
plotscatter(log10(pred_all[4,]), log10(CBnew$TSS), model = "ANN", r2 = "0.05", lab = c(FALSE, FALSE))
plotscatter(log10(pred_all[5,]), log10(CBnew$TSS), model = "Ondrusek et al., 2012", r2 = "0.37", 
            lab = c(TRUE, TRUE))
plotscatter(log10(pred_all[6,]), log10(CBnew$TSS), model = "Wang et al., 2009 with Eq. 9", r2 = "0.34", 
            lab = c(TRUE, FALSE))
mtext("Observed TSS", side = 1, outer= TRUE, line = 2)
mtext("Predicted TSS", side = 2, outer = TRUE, line = 2)

dev.off()

##Total w/ test dataset colored
testInd <- do.call(paste0, CBnew) %in% do.call(paste0, CBtest)
pointColors <- rep("black", nrow(CBnew))
pointColors[testInd] <- "red"

tiff(file = "TSS_total_colored.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

#windows(width=7.48, height=8.60)
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
layout(matrix(1:6,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  
plotscatter_col(log10(pred_all[1,]), log10(CBnew$TSS), model = "RVM", r2 = "0.48", lab = c(FALSE, TRUE), col = pointColors)
plotscatter_col(log10(pred_all[2,]), log10(CBnew$TSS), model = "SVR", r2 = "0.49", lab = c(FALSE, FALSE), col = pointColors)
plotscatter_col(log10(pred_all[3,]), log10(CBnew$TSS), model = "GLM", r2 = "0.51", lab = c(FALSE, TRUE), col = pointColors)
plotscatter_col(log10(pred_all[4,]), log10(CBnew$TSS), model = "ANN", r2 = "0.05", lab = c(FALSE, FALSE), col = pointColors)
plotscatter_col(log10(pred_all[5,]), log10(CBnew$TSS), model = "Ondrusek et al., 2012", r2 = "0.37", 
            lab = c(TRUE, TRUE), col = pointColors)
plotscatter_col(log10(pred_all[6,]), log10(CBnew$TSS), model = "Wang et al., 2009 with Eq. 9", r2 = "0.34", 
            lab = c(TRUE, FALSE), col = pointColors)
mtext("Observed TSS", side = 1, outer= TRUE, line = 2)
mtext("Predicted TSS", side = 2, outer = TRUE, line = 2)

dev.off()
