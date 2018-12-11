##################################################################################################################################
## PROGRAM: histograms.R
## SPECIFICATIONS: Plots histograms of Chlorophyll/TSS
## by chlorophyll and not stratified

##################################################################################################################################
#### Working Directory
rm(list=ls())
setwd("D:/Publication/")

require(e1071)
#####################################################################################################################################
## EDIT HERE
## Read in the entire dataset (cleaned) + Color Scale
load("modis pretty.RData")
CBchl <- read.csv("CB_chl.csv")
CBtss <- read.csv("CB_tss.csv")
CBchl$X <- NULL
CBtss$X <- NULL
######################################################################################################################
## Skewness for Chlorophyll and TSS

c(skewness(CBchl$chl),skewness(log10(CBchl$chl)))
c(skewness(CBtss$TSS),skewness(1/sqrt(CBtss$TSS)))

######################################################################################################################

##Sake of simplicity
portionGreater <- function(vec, limit, greater){
  if(greater){
    return(length(vec[vec > limit])/length(vec))
  }
  return(length(vec[vec <= limit])/length(vec))
}

##Approximately greater
maxCHL <- 60
maxTSS <- 60

portionGreater(CBchl$chl, maxCHL, greater=T)
portionGreater(CBtss$TSS, maxTSS, greater=T)

## Replace values greater than maxCHL/maxTSS
CBchl$chl[CBchl$chl > maxCHL] = maxCHL
CBtss$TSS[CBtss$TSS > maxTSS] = maxTSS

CBchl <- CBchl[order(CBchl$chl, decreasing = F),]
CBtss <- CBtss[order(CBtss$TSS, decreasing = F),]

#################################################################################################################################
## Plot Histograms of Chlorophyll

relFreq <- function(counts){
  return(counts/sum(counts))
}

# h1 <- hist(log10(CBchl$chl),  plot=F)
# h1$counts <- relFreq(h1$counts)
# 
# h2 <- hist(log10(CBtss$TSS), plot=F, breaks = h1$breaks)
# h2$counts <- relFreq(h2$counts)
h1 <- hist(CBchl$chl,  plot=F)
h1$counts <- relFreq(h1$counts)

h2 <- hist(CBtss$TSS, plot=F, breaks = h1$breaks)
h2$counts <- relFreq(h2$counts)

tiff(file = "histograms.tiff", width =7.5, height = 3.75, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1)) 
# plot(h1,xlab=expression(paste("Chl (log"["10"],"(mg/L))")),
#      ylab = "Relative Frequency", 
#      main="",cex.axis = 1.2,cex.lab=1.5,
#      ylim = c(0, 0.4), xlim=c(min(c(h1$breaks,h2$breaks)), 
#                                 max(c(h1$breaks,h2$breaks))))
# plot(h2,xlab=expression(paste("TSS (log"["10"],"(mg/L))")),
#      ylab = "Relative Frequency", 
#      main="",cex.axis = 1.2,cex.lab=1.5,
#      ylim = c(0, 0.4))
plot(h1,xlab=expression(Chlorophyll~a~(mu*g/L)),
     ylab = "Relative Frequency [%]", 
     main="",
     ylim = c(0, 0.6), xlim=c(min(c(h1$breaks,h2$breaks)), 
                              max(c(h1$breaks,h2$breaks))))
plot(h2,xlab="TSS (mg/L)",
     ylab = "", 
     main="",
     ylim = c(0, 0.6))
dev.off()
