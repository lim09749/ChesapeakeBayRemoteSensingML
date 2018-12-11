##################################################################################################################################
## PROGRAM: TSS scatterplots and metrics without running models.R
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
#setwd("D:/Publication")

setwd("H:/My Drive/Manuscripts/Li 2017/Cross Validation Procedures/TSS")

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

# Load Model Predictions
mat       = read.csv("TSS_model_metrics.csv", row.names = 1)
pred_test = read.csv("TSS_predict_test.csv", row.names = 1)
pred_all  = read.csv("TSS_predict_all.csv", row.names = 1)

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

##TSS

##Total w/ test dataset colored
testInd <- do.call(paste0, CBnew) %in% do.call(paste0, CBtest)
pointColors <- rep("black", nrow(CBnew))
pointColors[testInd] <- "red" ##Creates a vector with color names

##Save in tiff file
tiff(file = "TSS_total_colored.tiff", width =7.5, height = 10, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

#windows(width=7.48, height=8.60)
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
layout(matrix(1:6,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  
plotscatter_col(log10(as.numeric(pred_all[1,])), log10(CBnew$TSS), model = "RVM", r2 = "0.48",
                xlab = c("", "", "", "", "", ""), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
points(x=1, y=0, pch=16)
points(x=1, y=-0.25, pch=16, col="red")
text(x=1.1, y=0, labels  = c("Training Data"), adj=0)
text(x=1.1, y=-0.25, labels  = c("Testing Data"), adj=0)
plotscatter_col(log10(CBnew$TSS), log10(as.numeric(pred_all[2,])), model = "SVR", r2 = "0.49",
                xlab = c("", "", "", "", "", ""), ylab=c("", "", "", "", "", ""), col = pointColors)
plotscatter_col(log10(CBnew$TSS), log10(as.numeric(pred_all[3,])), model = "GLM", r2 = "0.51",
                xlab = c("", "", "", "", "", ""), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
plotscatter_col(log10(CBnew$TSS), log10(as.numeric(pred_all[4,])), model = "ANN", r2 = "0.05",
                xlab = c("", "", "", "", "", ""), ylab=c("", "", "", "", "", ""), col = pointColors)
plotscatter_col(log10(CBnew$TSS), log10(as.numeric(pred_all[5,])), model = "Ondrusek et al., 2012", r2 = "0.37", 
                xlab =  c(0.3, 1, 3.2, 10, 31.6, 100), ylab= c(0.3, 1, 3.2, 10, 31.6, 100), col = pointColors)
plotscatter_col(log10(CBnew$TSS), log10(as.numeric(pred_all[6,])), model = "Wang et al., 2009 with Eq. 9", r2 = "0.34", 
                xlab = c(0.3, 1, 3.2, 10, 31.6, 100), ylab=c("", "", "", "", "", ""), col = pointColors)
mtext("Observed TSS (mg/L)", side = 1, outer= TRUE, line = 2)
mtext("Predicted TSS (mg/L)", side = 2, outer = TRUE, line = 2)
dev.off()
#mtext(expression(Chlorophyll~a~(mu*g/L)), side=3, line=0)

