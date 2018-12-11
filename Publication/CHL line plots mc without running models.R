##################################################################################################################################
## PROGRAM: CHL line plots mc without running models.R
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
setwd("D:/Publication/")
CBnew <- read.csv("D:/Publication/CB_chl.csv")

##Read in data for upper, middle, and lower Bay
mon = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
upper_avg <- read.csv(file = "D:/Publication/CHL_upper_avg.csv", row.names = 1)
middle_avg <- read.csv(file = "D:/Publication/CHL_middle_avg.csv", row.names = 1)
lower_avg <- read.csv(file = "D:/Publication/CHL_lower_avg.csv", row.names = 1)
#####################################################################################################################################

##In-situ CHL for entire dataset (No need to have it cleaned)
##Get dates
mc_dates <- function(filename){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CB$DateTime),split = " "))[seq(1, length(CB$DateTime)*2,2)], 
                                       format = "%m/%d/%Y"), 6, 7))
  datenum
}

##Real data
mc_real <- function(filename, type){
  CB  <- na.omit(read.table(filename, header=T, strip.white = F, sep="\t"))
  datenum <- mc_dates(filename)
  
  mc_avg = rep(NA, 12)
  for(i in 1:12){
    ##Divide the data to different parts of the Bay
    if(type == "Upper"){
      mc_avg[i] = mean(CB$chl[datenum==i & CB$lat >= 38.75])
    }
    else if(type == "Middle"){
      mc_avg[i] = mean(CB$chl[datenum==i & CB$lat < 38.75 & CB$lat >= 37.75])
    }else if(type == "Lower"){
      mc_avg[i] = mean(CB$chl[datenum==i & CB$lat < 37.75])
    }
    else{
      mc_avg[i] = mean(CB$chl[datenum==i])
    }
  }
  mc_avg
}

##Save the data
filename = "D:/Publication/TSS and chl all stations Jul2002 to Dec2016 cleaned with MODISA.txt"
complete_real_upper <- mc_real(filename, type = "Upper")
complete_real_middle <- mc_real(filename, type = "Middle")
complete_real_lower <- mc_real(filename, type = "Lower")

##Save real data 
datenum <- as.numeric(substr(as.Date(unlist(strsplit(as.character(CBnew$DateTime),split = " "))
                                     [seq(1, nrow(CBnew)*2,2)],format = "%m/%d/%Y"), 6, 7))
real_upper <- rep(NA, 12)
real_middle <- rep(NA, 12)
real_lower <- rep(NA, 12)

##Load in real data
for(i in 1:12){
  real_upper[i] = mean(CBnew$chl[datenum==i & CBnew$lat >= 38.75])
  real_middle[i] = mean(CBnew$chl[datenum==i & CBnew$lat < 38.75 & CBnew$lat >= 37.75])
  real_lower[i] = mean(CBnew$chl[datenum==i & CBnew$lat < 37.75])
}  
#####################################################################################################################################

#Entire Dataset
#Cleaned Data
tiff(file = "chl_line_plots_entire.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(omi=c(0.25, 0.25, 0.25, 0.25),
    mai = c(0.5, 0.5, 0.05, 0.2))
layout(matrix(1:3,nrow=3,byrow=T), 
       heights = c(1,1,1)) # Height of rows (relative)  

##Upper Chesapeake Bay averages
plot(complete_real_upper, type = "l",
     xaxt = "n", xlab='', ylab = expression(Chlorophyll~a~(mu*g/L)), 
     ylim = c(0, 25), lwd = 2,
     axes =F)
points(complete_real_upper, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(upper_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(upper_avg[1,]), lty = "dashed")
points(as.numeric(upper_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(upper_avg[2,]), lty = "dotted")
points(as.numeric(upper_avg[3,])*10.41-12.82, pch = 18, cex = 1.5)
lines(as.numeric(upper_avg[3,])*10.41-12.82, lty = "dotdash", pch = 17)
points(as.numeric(upper_avg[4,])*0.01678+9.85917, pch = 0, cex = 1.5)
lines(as.numeric(upper_avg[4,])*0.01678+9.85917, lty = "longdash", pch = 25)
text(6.5, 22.5, "Upper Chesapeake Bay", cex = 1.5)
legend(11, 22.5, legend= c("In-situ", "RVM", "SVR", "GLM", "ANN"),
       lty=c("solid","dashed", "dotted", "dotdash","longdash"),
       pch = c(15,16,17,18,0))
##Middle Chesapeake Bay plots
plot(complete_real_middle, type = "l",
     xaxt = "n", xlab='', ylab = expression(Chlorophyll~a~(mu*g/L)), 
     ylim = c(0, 25), lwd = 2,
     axes =F)
points(complete_real_middle, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(middle_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(middle_avg[1,]), lty = "dashed")
points(as.numeric(middle_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(middle_avg[2,]), lty = "dotted")
points(as.numeric(middle_avg[3,])*10.41-12.82, pch = 18, cex = 1.5)
lines(as.numeric(middle_avg[3,])*10.41-12.82, lty = "dotdash", pch = 17)
points(as.numeric(middle_avg[4,])*0.01678+9.85917, pch = 0, cex = 1.5)
lines(as.numeric(middle_avg[4,])*0.01678+9.85917, lty = "longdash", pch = 25)
text(6.5, 22.5, "Middle Chesapeake Bay", cex = 1.5)

##Lower Chesapeake Bay averages
plot(complete_real_lower, type = "l",
     xaxt = "n", xlab='Months', ylab = expression(Chlorophyll~a~(mu*g/L)), 
     ylim = c(0, 25), lwd = 2,
     axes =F)
points(complete_real_lower, lty = "solid", pch = 15, cex = 1.5)
axis(1, at = 1:12, labels=mon)
axis(2)
points(as.numeric(lower_avg[1,]), pch = 16, cex = 1.5)
lines(as.numeric(lower_avg[1,]), lty = "dashed")
points(as.numeric(lower_avg[2,]), pch = 17, cex = 1.5)
lines(as.numeric(lower_avg[2,]), lty = "dotted")
points(as.numeric(lower_avg[3,])*10.41-12.82, pch = 18, cex = 1.5)
lines(as.numeric(lower_avg[3,])*10.41-12.82, lty = "dotdash", pch = 17)
points(as.numeric(lower_avg[4,])*0.01678+9.85917, pch = 0, cex = 1.5)
lines(as.numeric(lower_avg[4,])*0.01678+9.85917, lty = "longdash", pch = 25)
text(6.5, 22.5, "Lower Chesapeake Bay", cex = 1.5)
dev.off()
#####################################################################################################################################
