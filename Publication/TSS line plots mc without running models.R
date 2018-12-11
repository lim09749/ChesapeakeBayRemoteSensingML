##################################################################################################################################
## PROGRAM: TSS line plots mc without running models.R
## SPECIFICATIONS: Plots line plots for monthly averages
#####################################################################################################################################

rm(list=ls())
#####################################################################################################################################
## Libraries
require(e1071)      ## support vector regression
require(kernlab)    ## relevance vector machines
require(neuralnet)  ## neural networks

#####################################################################################################################################
##Setting up programs
#setwd("H:/My Drive/Manuscripts/Li 2017/Cross Validation Procedures/TSS")
setwd("D:/Publication/")
CBnew <- read.csv("D:/Publication/CB_tss.csv")

mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
upper_avg <- read.csv(file = "D:/Publication/TSS_upper_avg.csv", row.names = 1)
middle_avg <- read.csv(file = "D:/Publication/TSS_middle_avg.csv", row.names = 1)
lower_avg <- read.csv(file = "D:/Publication/TSS_lower_avg.csv", row.names = 1)
#####################################################################################################################################

##In-situ TSS for entire dataset (No need to have it cleaned)
##Get dates
mc_dates <- function(filename){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CB$DateTime),split = " "))[seq(1, length(CB$DateTime)*2,2)], 
                                       format = "%m/%d/%Y"), 6, 7))
  datenum
}

##Get in-situ data
mc_real <- function(filename, type){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- mc_dates(filename)
  
  ##Divide it into upper, middle, lower datasets
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

##Find in-situ upper, middle, lower averages
filename = "D:/Publication/TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt"
complete_real_upper <- mc_real(filename, type = "Upper")
complete_real_middle <- mc_real(filename, type = "Middle")
complete_real_lower <- mc_real(filename, type = "Lower")

##Create vec for data points
datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CBnew$DateTime),split = " "))
                                     [seq(1, nrow(CBnew)*2,2)],format = "%m/%d/%Y"), 6, 7))
real_upper <- rep(NA, 12)
real_middle <- rep(NA, 12)
real_lower <- rep(NA, 12)

##Read it in
for(i in 1:12){
  real_upper[i] = mean(CBnew$TSS[datenum==i & CBnew$lat >= 38.75])
  real_middle[i] = mean(CBnew$TSS[datenum==i & CBnew$lat < 38.75 & CBnew$lat >= 37.75])
  real_lower[i] = mean(CBnew$TSS[datenum==i & CBnew$lat < 37.75])
}  
#####################################################################################################################################
##Start plotting

#Cleaned Data
# 
# tiff(file = "TSS_line_plots_cleaned.tiff", width =7.5, height = 7.5, units = "in",
#      pointsize=10, res = 300, compression = c("lzw"))
# par(mfrow=c(3,1))
# 
# plot(real_upper, type = "l",
#      xaxt = "n", xlab='Months', ylim = c(0, 18), lwd = 2)
# points(real_upper, lty = "solid", pch = 15, cex = 1.5)
# axis(1, at = 1:12, labels=mon)
# 
# points(as.numeric(upper_avg[1,]), pch = 16, cex = 1.5)
# lines(as.numeric(upper_avg[1,]), lty = "dashed")
# 
# points(as.numeric(upper_avg[3,])*10.26-12.99, pch = 17, cex = 1.5)
# lines(as.numeric(upper_avg[3,])*10.26-12.99, lty = "dashed", pch = 17)
# 
# points(as.numeric(upper_avg[4,])*-5.163e-06+8.591e+00, pch = 18, cex = 1.5)
# lines(as.numeric(upper_avg[4,])*-5.163e-06+8.591e+00, 
#       lty = "dashed", pch = 18)
# 
# plot(real_middle, type = "l",
#      xaxt = "n", xlab='Months', ylim = c(0, 18), lwd = 2)
# points(real_middle, lty = "solid", pch = 15, cex = 1.5)
# axis(1, at = 1:12, labels=mon)
# 
# points(as.numeric(middle_avg[1,]), pch = 16, cex = 1.5)
# lines(as.numeric(middle_avg[1,]), lty = "dashed")
# 
# points(as.numeric(middle_avg[3,])*10.26-12.99, pch = 17, cex = 1.5)
# lines(as.numeric(middle_avg[3,])*10.26-12.99, lty = "dashed", pch = 17)
# 
# points(as.numeric(middle_avg[4,])*-5.163e-06+8.591e+00, pch = 18, cex = 1.5)
# lines(as.numeric(middle_avg[4,])*-5.163e-06+8.591e+00, 
#       lty = "dashed", pch = 18)
# 
# plot(real_lower, type = "l",
#      xaxt = "n", xlab='Months', ylim = c(0, 18), lwd = 2)
# points(real_lower, lty = "solid", pch = 15, cex = 1.5)
# axis(1, at = 1:12, labels=mon)
# 
# points(as.numeric(lower_avg[1,]), pch = 16, cex = 1.5)
# lines(as.numeric(lower_avg[1,]), lty = "dashed")
# 
# points(as.numeric(lower_avg[3,])*10.26-12.99, pch = 17, cex = 1.5)
# lines(as.numeric(lower_avg[3,])*10.26-12.99, lty = "dashed", pch = 17)
# 
# points(as.numeric(lower_avg[4,])*-5.163e-06+8.591e+00, pch = 18, cex = 1.5)
# lines(as.numeric(lower_avg[4,])*-5.163e-06+8.591e+00, 
#       lty = "dashed", pch = 18)
# dev.off()
#####################################################################################################################################

#Entire Dataset
#Cleaned Data
tiff(file = "TSS_line_plots_entire.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(omi=c(0.25, 0.25, 0.25, 0.25),
    mai = c(0.5, 0.5, 0.05, 0.2))
layout(matrix(1:3,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  

plot(complete_real_upper, type = "l",
     xaxt = "n", xlab='', ylab = "TSS (mg/L)", ylim = c(0, 18), lwd = 2,
     axes =F)
points(complete_real_upper, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(upper_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(upper_avg[1,]), lty = "dashed")
points(as.numeric(upper_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(upper_avg[2,]), lty = "dotted")
points(as.numeric(upper_avg[3,])*10.26-12.99, pch = 18, cex = 1.5)
lines(as.numeric(upper_avg[3,])*10.26-12.99, lty = "dotdash", pch = 17)
points(as.numeric(upper_avg[4,])*-5.163e-06+8.591e+00, pch = 0, cex = 1.5)
lines(as.numeric(upper_avg[4,])*-5.163e-06+8.591e+00, lty = "longdash", pch = 25)
text(6.5, 16, "Upper Chesapeake Bay", cex = 1.5)
legend(11, 7, legend= c("In-situ", "RVM", "SVR", "GLM", "ANN"),
       lty=c("solid","dashed", "dotted", "dotdash","longdash"),
       pch = c(15,16,17,18,0))

plot(complete_real_middle, type = "l",
     xaxt = "n", xlab='', ylab = "TSS (mg/L)", ylim = c(0, 18), lwd = 2,
     axes =F)
points(complete_real_middle, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(middle_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(middle_avg[1,]), lty = "dashed")
points(as.numeric(middle_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(middle_avg[2,]), lty = "dotted")
points(as.numeric(middle_avg[3,])*10.26-12.99, pch = 18, cex = 1.5)
lines(as.numeric(middle_avg[3,])*10.26-12.99, lty = "dotdash", pch = 17)
points(as.numeric(middle_avg[4,])*-5.163e-06+8.591e+00, pch = 0, cex = 1.5)
lines(as.numeric(middle_avg[4,])*-5.163e-06+8.591e+00, lty = "longdash", pch = 25)
text(6.5, 16, "Middle Chesapeake Bay", cex = 1.5)
plot(complete_real_lower, type = "l",
     xaxt = "n", xlab='Months', ylab = "TSS (mg/L)", ylim = c(0, 18), lwd = 2,
     axes =F)
points(complete_real_lower, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(lower_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(lower_avg[1,]), lty = "dashed")
points(as.numeric(lower_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(lower_avg[2,]), lty = "dotted")
points(as.numeric(lower_avg[3,])*10.26-12.99, pch = 18, cex = 1.5)
lines(as.numeric(lower_avg[3,])*10.26-12.99, lty = "dotdash", pch = 17)
points(as.numeric(lower_avg[4,])*-5.163e-06+8.591e+00, pch = 0, cex = 1.5)
lines(as.numeric(lower_avg[4,])*-5.163e-06+8.591e+00, lty = "longdash", pch = 25)
text(6.5, 16, "Lower Chesapeake Bay", cex = 1.5)
dev.off()
#####################################################################################################################################
