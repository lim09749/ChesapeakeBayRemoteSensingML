##################################################################################################################################
## PROGRAM: spectrum plots.R
## SPECIFICATIONS: Plots histograms of Chlorophyll/TSS, log and without log; Plots Spectrum of remote sensing reflectance, stratified
## by chlorophyll and not stratified

##################################################################################################################################
#### Working Directory
rm(list=ls())
setwd("D:/Publication/")

#####################################################################################################################################
## EDIT HERE
spectrum <- c(412,443,469,488,531,547,555,645,667,678) # Wavelengths given
n<-255 # Number of subdivisions 

## Read in the entire dataset (cleaned) + Color Scale
load("modis pretty.RData")
CBchl <- read.csv("CB_chl.csv")
CBtss <- read.csv("CB_tss.csv")
CBchl$X <- NULL
CBtss$X <- NULL

##Use entire (uncleaned) dataset for avg rrs
rrs <- na.omit(read.table("TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt", 
                          header=T, strip.white = F, sep="\t"))
##
mean_rrs <- colMeans(rrs[,10:ncol(rrs)])
colsd <- function(CBnew){
  csd <- rep(NA, ncol(CBnew))
  for(i in 1:ncol(CBnew)){
    csd[i] <- sd(CBnew[,i])
  }
  csd
}
sd_rrs <- colsd(rrs[,10:ncol(rrs)])

# # Plot Mean +- sd
# plot(spectrum, as.matrix(mean_rrs),type="l",
#      xlab = expression(Wavelength~lambda~(nm)), 
#      ylab = expression(Remote~sensing~reflectance~R[rs](lambda)~(sr^-1)),
#      xlim = c(400,700),
#      ylim= c(0.00, max(mean_rrs+sd_rrs)),
#      axes = FALSE)
# lines(spectrum, mean_rrs+sd(mean_rrs),lty=2,lwd=1)
# lines(spectrum, mean_rrs-sd(mean_rrs),lty=2,lwd=1)
# legend(x=400, y=0.009, legend=c(expression(mu["rrs"]),expression(paste(mu["rrs"],"\u00B1",sigma["rrs"]))),
#        lty=c(1,2), cex=1, bty="n")
# axis(1, at = seq(400, 700, 50))
# axis(2)
# rect(400,0.008,440,0.009)

######################################################################################################################
##Sake of simplicity

##Approximately greater
maxCHL <- 60
maxTSS <- 60

## Replace values greater than maxCHL/maxTSS
CBchl$chl[CBchl$chl > maxCHL] = maxCHL
CBtss$TSS[CBtss$TSS > maxTSS] = maxTSS

CBchl <- CBchl[order(CBchl$chl, decreasing = F),]
CBtss <- CBtss[order(CBtss$TSS, decreasing = F),]

# ######################################################################################################################
# ## Subdivides the chlorophyll values based on their value: returns the cutted 
# # tss or chlorophyll 
# # Precondition: length(vec)>0, length(spectrum)>0, length(mycol)=n, vec=log10(chl or TSS)
# cutValues <- function(n, vec){
#   vec <- vec[order(vec, decreasing=F)]
#   c   <- cut(vec, n)
#   return(c)
# }
# #####################################################################################################################
# ## Return corresponding color vector
# returnColorsVec <- function(n, vec, mycol){# Prereq: length(myvol)=n, 
#   # vec in order (increasing)
#   if(!is.unsorted(vec)){
#     cutVal <- cutValues(n, vec)
#     colSeq <- rep("", length(vec)) ## Create a vector of the colors
#   
#     colSeq <- sapply(1:length(vec), 
#                      FUN=function(x) return(mycol[as.integer(cutVal[x])]),simplify= TRUE)
#     return(colSeq)
#   }
# }  
# ######################################################################################################################
# #### Plot the spectrum plot
# 
# plotSpectrum <- function(n, CBnew, vec, mycol){
#   colSeq <- returnColorsVec(n, vec, mycol)
#   plot(9001,1, 
#        xlim=c(412, 678), 
#        ylim=c(0, 0.025),
#        xlab = "Wavelength (nm)", 
#        ylab = "Remote sensing reflectance"
#        #,main=paste0("Spectrum Sub-divided by Chl concentrations (",n," groups)")
#   )
#   for(i in 1:nrow(CBnew)){
#     lines(spectrum, c(CBnew$rrs_412[i], CBnew$rrs_443[i], CBnew$rrs_469[i],
#                       CBnew$rrs_488[i], CBnew$rrs_531[i], CBnew$rrs_547[i],
#                       CBnew$rrs_555[i], CBnew$rrs_645[i], CBnew$rrs_667[i],
#                       CBnew$rrs_678[i]), col = colSeq[i])
#   }
# }
# ######################################################################################################################
# old.par <- par(mar = c(0, 0, 0, 0))
# par(old.par)
# plotSpectrum(n, CBchl, log10(CBchl$chl), mycol)
# plotKey(n, log10(CBchl$chl), mycol)
# 
# par(old.par)
# plotSpectrum(n, CBtss, log10(CBtss$TSS), mycol)
# plotKey(n, log10(CBtss$TSS), mycol)
# 
# #### Create a key
# 
#  #bottom, left, top, right
# barplot(rep(1,40), col=viridis(40), border=viridis(40), space=0, axes=F)
# axis(side=4, at = c(0, 10, 20, 30, 40), 
#      labels=c("0", "7.5", "15", "22.5", "30"))
# #return(labels)
# 
# 
# plotKey(n, log10(CBchl$chl), mycol)
# plotKey(n, log10(CBtss$TSS), mycol)
#######################################################################################################################
###Spectrum plots
tiff(file = "spectrum.tiff", width =7.5, height = 6.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(omi=c(0.25, 0.25, 0.25, 0.25),
    mai = c(0.5, 0.5, 0.05, 0.2))
layout(matrix(c(1,1:5),nrow=3,byrow=T), 
       heights = c(4,4,1)) # Height of rows (relative)  

##Spectrum rrs
# Plot Mean +- sd
plot(spectrum, as.matrix(mean_rrs),type="l",
     xlab = expression(Wavelength~lambda~(nm)), 
     ylab = expression(Remote~sensing~reflectance~R[rs](lambda)~(sr^-1)),
     xlim = c(400,700),
     ylim= c(0.00, max(mean_rrs+sd_rrs)),
     axes = FALSE)


lines(spectrum, mean_rrs+sd(mean_rrs),lty=2,lwd=1)
lines(spectrum, mean_rrs-sd(mean_rrs),lty=2,lwd=1)
legend(x=400, y=0.009, legend=c(expression(mu["rrs"]),expression(paste(mu["rrs"],"\u00B1",sigma["rrs"]))),
       lty=c(1,2), cex=1, bty="n")
axis(1, at = seq(400, 700, 50))
axis(2)
rect(400,0.007,440,0.009)

##TSS
tss_seq <-seq(from=0.01, to=maxTSS, length.out = n+1)
tss_ind <- unlist(lapply(CBtss$TSS, FUN = function(x) which(x-tss_seq <= 0)[1]-1))
tss_col <- mycol[tss_ind]

plot(9001,1, 
     xlim=c(400, 700), 
     ylim=c(0, 0.025),
     xlab = expression(Wavelength~lambda~(nm)), 
     ylab = expression(Remote~sensing~reflectance~R[rs](lambda)~(sr^-1)),
     axes = F)

axis(1, at = seq(400, 700, 50))
axis(2)

for(i in 1:nrow(CBtss)){
  lines(spectrum, c(CBtss$Rrs_412[i], CBtss$Rrs_443[i], CBtss$Rrs_469[i],
                    CBtss$Rrs_488[i], CBtss$Rrs_531[i], CBtss$Rrs_547[i],
                    CBtss$Rrs_555[i], CBtss$Rrs_645[i], CBtss$Rrs_667[i],
                    CBtss$Rrs_678[i]), col = tss_col[i])
}

##CHL
chl_seq <-seq(from=0.01, to=maxCHL, length.out = n+1)
chl_ind <- unlist(lapply(CBchl$chl, FUN = function(x) which(x-chl_seq <= 0)[1]-1))
chl_col <- mycol[chl_ind]

plot(9001,1, 
     xlim=c(400, 700), 
     ylim=c(0, 0.025),
     xlab = expression(Wavelength~lambda~(nm)), 
     ylab = "",
     #ylab = expression(Remote~sensing~reflectance~R[rs](lambda)~(sr^-1)),
     axes = F)


axis(1, at = seq(400, 700, 50))
axis(2)

for(i in 1:nrow(CBchl)){
  lines(spectrum, c(CBchl$Rrs_412[i], CBchl$Rrs_443[i], CBchl$Rrs_469[i],
                    CBchl$Rrs_488[i], CBchl$Rrs_531[i], CBchl$Rrs_547[i],
                    CBchl$Rrs_555[i], CBchl$Rrs_645[i], CBchl$Rrs_667[i],
                    CBchl$Rrs_678[i]), col = chl_col[i])
}

##Barplots
barplot(rep(1,255), col=mycol, border=mycol, space=0, axes=F)
axis(side=1, at = c(0, 64, 128, 192, 255), 
     labels=c("0", "15", "30", "45", "60"))
mtext(side=1, "TSS (mg/L)", line=3, cex=1.1)

barplot(rep(1,255), col=mycol, border=mycol, space=0, axes=F)
axis(side=1, at = c(0, 64, 128, 192, 255), 
     labels=c("0", "15", "30", "45", "60"))
mtext(side = 1, expression(Chlorophyll~a~(mu*g/L)), line = 3, cex = 1.1)
dev.off()
