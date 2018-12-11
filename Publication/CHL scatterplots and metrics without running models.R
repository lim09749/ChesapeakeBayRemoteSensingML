##################################################################################################################################
## PROGRAM: CHL scatterplots and metrics without running models.R
## SPECIFICATIONS: plots scatterplots
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
#setwd("H:/My Drive/Manuscripts/Li 2017/Cross Validation Procedures/CHL")

pattern <- c("rvm","svm", "glm", "ann")

##Model variables
model_var <- c("chl","Rrs_412","Rrs_443","Rrs_469","Rrs_488","Rrs_531",
               "Rrs_547","Rrs_555","Rrs_645","Rrs_667","Rrs_678")

##Load test data
CBtest <- read.csv("CHL_test.csv")
CBtest <- CBtest[,colnames(CBtest)%in%model_var]

##load entire dataset
CBnew <- read.csv("CB_chl.csv")
CBnew <- CBnew[,colnames(CBnew)%in%model_var]

# Load Model Predictions
mat       = read.csv("CHL_model_metrics.csv", row.names = 1)
pred_test = read.csv("CHL_predict_test.csv", row.names = 1)
pred_all  = read.csv("CHL_predict_all.csv", row.names = 1)

#####################################################################################################################################
##Scatterplot functions

plotscatter <- function(x, y, model, r2, xlab, ylab){
  par(mai=c(0.1,0.1,0.1,0.1)) 
  plot(x, y, xlim = c(-0.5, 2), ylim = c(-0.5, 2), pch=16, axes=FALSE)
  text(x  = 0.1, y = 1.7,labels=substitute(model))
  text(x  = 0.1, y = 1.5,labels=substitute(paste("r"[2]," = ", r2)))
  clip(-0.5, 2, -0.5, 2)
  abline(a=0, b=1)
  axis(1, at =c(-0.5, 0, 0.5, 1, 1.5, 2), labels = xlab)
  axis(2, at =c(-0.5, 0, 0.5, 1, 1.5, 2), labels = ylab)
}


plotscatter_col <- function(x, y, model, r2, xlab, ylab, col){
  par(mai=c(0.1,0.1,0.1,0.1)) 
  plot(x, y, xlim = c(-0.5, 2), ylim = c(-0.5, 2), col = col, axes = FALSE, pch=16)
  text(x  = 0.1, y = 1.7,labels=substitute(model))
  text(x  = 0.1, y = 1.5,labels=substitute(paste("r"[2]," = ", r2)))
  clip(-0.5, 2, -0.5, 2)
  abline(a=0, b=1)
  axis(1, at =c(-0.5, 0, 0.5, 1, 1.5, 2), labels = xlab)
  axis(2, at =c(-0.5, 0, 0.5, 1, 1.5, 2), labels = ylab)
}
#####################################################################################################################################
##Plotting 

##Total w/ test dataset colored
testInd <- do.call(paste0, CBnew) %in% do.call(paste0, CBtest)
pointColors <- rep("black", nrow(CBnew))
pointColors[testInd] <- "red"

tiff(file = "CHL_total_colored.tiff", width =7.5, height = 10, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
#windows(width=7.48, height=8.60)
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
layout(matrix(c(rep(1:4,each=2),5,6,6,7),nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  
# layout(matrix(1:6,nrow=3,byrow=T), 
#        heights = c(1,1,1)) # Height of rows (relative)  
plotscatter_col(log10(as.numeric(pred_all[1,])), log10(CBnew$chl), model = "RVM", r2 = "0.51",
                xlab = c("", "", "", "", "", ""), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
points(x=1, y=0, pch=16)
points(x=1, y=-0.25, pch=16, col="red")
text(x=1.1, y=0, labels  = c("Training Data"), adj=0)
text(x=1.1, y=-0.25, labels  = c("Testing Data"), adj=0)
plotscatter_col(log10(CBnew$chl), log10(as.numeric(pred_all[2,])), model = "SVR", r2 = "0.44",
                xlab = c("", "", "", "", "", ""), ylab=c("", "", "", "", "", ""), col = pointColors)
plotscatter_col(log10(CBnew$chl), log10(as.numeric(pred_all[3,])), model = "GLM", r2 = "0.48",
                xlab = c("", "", "", "", "", ""), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
plotscatter_col(log10(CBnew$chl), log10(as.numeric(pred_all[4,])), model = "ANN", r2 = "0.40",
                xlab = c("", "", "", "", "", ""), ylab=c("", "", "", "", "", ""), col = pointColors)

plot.new()
plotscatter_col(log10(CBnew$chl), log10(as.numeric(pred_all[5,])), model = "NASA OCM3", r2 = "0.09", 
                xlab =  c(0.3, 1, 3.2, 10, 31.6, 100), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
mtext(expression(Observed~Chlorophyll~a~(mu*g/L)), side = 1, outer= TRUE, line = 2)
mtext(expression(Predicted~Chlorophyll~a~(mu*g/L)), side = 2, outer = TRUE, line = 2)
plot.new()
dev.off()

