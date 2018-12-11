rm(list=ls())

require(rgdal)
require(alphahull)
require(SearchTrees)
require(e1071)
#####################################################################################################################################
rse <- function(error){return(sqrt(error^2))} #Makes error values positive

rmse <- function(error){return(sqrt(mean(error^2)))} #Squares the error, finds the mean, 
#and then finds the square root of
#the mean (RMSE)

#Finds the interquartile mean
iqm <- function(values){return(mean(values[(values>=summary(values)[2])&
                                             (values <= summary(values)[5])]))}

######################################################################################################################
#Load TSS with MODIS Data

CB  <- read.table("E:/Marvin/Research 2016/SatData/TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt", header=T, strip.white = F, sep="\t")
CB<-na.omit(CB)


#Convert Distance to meters by converting lat lon to UTM (Station)
xy              <- data.frame(ID = c(CB$EventId), X = c(CB$lon), Y = c(CB$lat))

coordinates(xy) <- c("X", "Y")

proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example

coordStation             <- as.data.frame(spTransform(xy, CRS("+proj=utm +zone=18 ellps=WGS84")))
coordStation$ID <- NULL

coordStation <-as.matrix(coordStation, nrow=2, ncol=ncol(coordStation)) #Coordinates of TSS stations


#Convert Distance to meters by converting lat lon to UTM (Satellite)
xySat              <- data.frame(X = c(CB$MODIS.lon), Y = c(CB$MODIS.lat))

coordinates(xySat) <- c("X", "Y")

proj4string(xySat) <- CRS("+proj=longlat +datum=WGS84")  ## for example

coordStationSat             <- as.data.frame(spTransform(xySat, CRS("+proj=utm +zone=18 ellps=WGS84")))
coordStationSat <-as.matrix(coordStationSat, nrow=2, ncol=ncol(coordStationSat)) #Coordinates of TSS stations


######################################################################################################################
#Open Shapefile
setwd("E:/Marvin/Research 2016/SatData/Baybox CB Shoreline")
shp = "baybox_dd.shp"
myshp = readOGR(shp, layer = basename(strsplit(shp, "\\.")[[1]])[1]) # This is a fancy way of being lazy, so I do not need to type the layer name in

#Dump shapefile coordinates into 'coords'
polys = attr(myshp,'polygons')
npolys = length(polys)
for (i in 1:npolys){
  poly = polys[[i]]
  polys2 = attr(poly,'Polygons')
  npolys2 = length(polys2)
  for (j in 1:npolys2){
    #do stuff with these values
    if (i == 1 & j==1) {
      coords =  coordinates(polys2[[j]])
    }else{
      coords =  rbind(coords, coordinates(polys2[[j]]))
    }
  }
}

#Convert coordinates to UTM
coords <- project(coords, "+proj=utm +zone=18 ellps=WGS84")  #Coordinates of figure
######################################################################################################################
#Find closest point to shore to station
#I know the simplest codingway is just to use a lot of for loops, 
#but I decided to use a quadtree since it takes log(n) time instead
#of n time making it around X69,399.71 faster. It is basically a binary tree search 
#but with coordinates.
#install.packages("SearchTrees")

treeQuadCoords <- createTree(coords, treeType="quad",dataType="point",maxDepth=25)
#Create a quadtree
######################################################################################################################

inds         <- knnLookup(treeQuadCoords,newdat=coordStation,k=2)#Lookup 2 of the nearest neighbors
#In "inds", each ROW NUMBER corresponds to the TSS station ROW NUMBER while 
#the VALUE in the row (NOT the column number) corresponds to the ROW NUMBER

CoordCompare <- cbind(coordStation,coords[inds[,1],],coords[inds[,2],]) #Combine the TSS and nearest point
#The first two rows explain the coordinates of the TSS station, while the every two points
#represents the nearest points

dist1        <- sqrt((coordStation[,1]-coords[inds[,1],1])^2+(coordStation[,2]-coords[inds[,1],2])^2)
dist2        <- sqrt((coordStation[,1]-coords[inds[,2],1])^2+(coordStation[,2]-coords[inds[,2],2])^2)
#Find the distance

CoordDist   <- cbind(coordStation, dist1,dist2) #Combine them

######################################################################################################################

indsSat         <-  knnLookup(treeQuadCoords,newdat=coordStationSat,k=2)#Lookup 2 of the nearest neighbors
#In "inds", each ROW NUMBER corresponds to the Satellite Station ROW NUMBER while 
#the VALUE in the row (NOT the column number) corresponds to the ROW NUMBER

CoordCompareSat <- cbind(coordStationSat,coords[indsSat[,1],],coords[indsSat[,2],]) #Combine the TSS and nearest point
#The first two rows explain the coordinates of the Satellite Station, while the every two points
#represents the nearest points

dist1Sat        <- sqrt((coordStationSat[,1]-coords[indsSat[,1],1])^2+(coordStationSat[,2]
                                                                       -coords[indsSat[,1],2])^2)
dist2Sat        <- sqrt((coordStationSat[,1]-coords[indsSat[,2],1])^2+(coordStationSat[,2]
                                                                       -coords[indsSat[,2],2])^2)
#Find the distance

CoordDistSat    <- cbind(coordStationSat, dist1Sat,dist2Sat) #Combine them

######################################################################################################################
distSatStation <- sqrt((coordStationSat[,1]-coordStation[,1])^2
                       +(coordStationSat[,2]-coordStation[,2])^2) #Finds the distance 
#from satellite to station
######################################################################################################################
CBnew <- CB[,c(3,10:19)] #Only take the useful info
#Finds the points with a distance of less than 200 meters from satellite to station
#and 200 meters offshore
CBnew2 <- CBnew[((distSatStation<200)&(CoordDist[,1]>1000)),] 


#Divide into training and testing dataset
first_mark  <- round(0.8*nrow(CBnew2))
######################################################################################################################
plot(log10(CBnew2[,2]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,2])))
plot(log10(CBnew2[,3]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,3])))
plot(log10(CBnew2[,4]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,4])))
plot(log10(CBnew2[,5]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,5])))
plot(log10(CBnew2[,6]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,6])))
plot(log10(CBnew2[,7]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,7])))
plot(log10(CBnew2[,8]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,8])))
plot(log10(CBnew2[,9]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,9])))
plot(log10(CBnew2[,10]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,10])))
plot(log10(CBnew2[,11]), log10(CBnew2$TSS))
summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,11])))
######################################################################################################################
max<-0
k<-0
divide_list <- matrix(data=NA, nrow=16, ncol=3)

for(i in 2:ncol(CBnew2)){
  for (j in 2:ncol(CBnew2)){
    
    if ((as.double(as.character(summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,i]/CBnew2[,j]))))[9]))>0.20){
      divide_list[k+1,] <- c((as.double(as.character(summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,i]/CBnew2[,j]))))[9])),
                             colnames(CBnew2)[i],colnames(CBnew2)[j])
      print(colnames(CBnew2)[i])
      print(colnames(CBnew2)[j])
      plot(CBnew2$TSS,CBnew2[,i]/CBnew2[,j],log="xy")
      print(as.double(as.character(summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,i]/CBnew2[,j]))))[9]))
      k <- k+1
    }
    if((as.double(as.character(summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,i]/CBnew2[,j]))))[9]))>max){
      max <- (as.double(as.character(summary(lm(log10(CBnew2$TSS)~log10(CBnew2[,i]/CBnew2[,j]))))[9]))
    }
  }
}
print(max)

(divide_list[order(as.double(divide_list[,1]),decreasing=T),])
######################################################################################################################
CBnew2$Rrs_531Rrs_667 <- CBnew2$Rrs_531/CBnew2$Rrs_667

CBnew2$Rrs_547Rrs_667 <- CBnew2$Rrs_547/CBnew2$Rrs_667

CBnew2$Rrs_531Rrs_645 <- CBnew2$Rrs_531/CBnew2$Rrs_645

CBnew2$Rrs_547Rrs_645 <- CBnew2$Rrs_547/CBnew2$Rrs_645

CBnew2$Rrs_555Rrs_667 <- CBnew2$Rrs_555/CBnew2$Rrs_667

CBnew2$Rrs_555Rrs_645 <- CBnew2$Rrs_555/CBnew2$Rrs_645

CBnew2$Rrs_488Rrs_667 <- CBnew2$Rrs_488/CBnew2$Rrs_667

CBnew2$Rrs_531Rrs_678 <- CBnew2$Rrs_531/CBnew2$Rrs_678
######################################################################################################################

#Randomly arrange the row numbers of the pruned dataset
randomRows  <- sample.int(nrow(CBnew2), 
                          nrow(CBnew2))

#Organize into training and testing dataset
CBtrain     <- CBnew2[randomRows[1:first_mark],]
CBtest      <- CBnew2[randomRows[(first_mark+1):length(randomRows)],]


######################################################################################################################
#Set the max repetition
epochReal<-10^9
h2o.init(nthreads=-1)
h2o.removeAll()
h2o.no_progress()

CBtrain <- as.h2o(CBtrain)
CBtest <- as.h2o(CBtest)
######################################################################################################################

TSSnet2<-h2o.deeplearning(x=6:19,y=1, training_frame=CBtrain, #Train the ANN with all the data (total_data)
                          overwrite_with_best_model=FALSE, #Make sure it gets the best model
                          activation = "TanhWithDropout", #Specify the activation function 
                          #(Once weights are multiplied by the inputs
                          #and bias, the activation function is applied)
                          hidden = c(1000,500,250), #Specify the hidden layers (3 hidden 
                          #layers between input and output each with 500 neurons)
                          export_weights_and_biases=TRUE,#Set it so I can export weights 
                          variable_importances = TRUE,#Set it so I can extract variable improtances
                          #and neurons
                          epochs = epochReal) #Set the repetitions (How many times 
#it uses chain rule and gradient descent to adjust results)

#####################################################################################################################################
#Create a dataframe with all of the predicted and real data

predictANN = h2o.predict(object = TSSnet2, newdata = as.h2o(CBtest[c(-1)])) #Predict the results
comparisonANN <- data.frame(actual = (as.data.frame(as.h2o(CBtest[c(1)]))), prediction =(as.data.frame(predictANN))) #Put them in the same dataset
plot(comparisonANN$predict, comparisonANN$TSS, log="xy", xlab="Predicted TSS",ylab="Real TSS",main="Comparison of Real and Predicted TSS values")
abline(lm(TSS~predict,data=log10(comparisonANN)))
#####################################################################################################################################
plot(TSSnet2)
print(as.double(as.character(summary(lm(TSS~predict,data=log10(comparisonANN)))[9])))
#Print the r^2 value of the regression

#####################################################################################################################################
h2o.saveModel(TSSnet2, path = "E:\\Marvin\\Research 2016\\ANN", force = T)
#Save the model
#####################################################################################################################################
#Plots
plot(TSSnet2)#Plot the error function overtime

#Plot the predicted vs real TSS
plot(comparisonANN$predict, comparisonANN$TSS, xlab="Predicted TSS",ylab="Real TSS",main="Comparison of Real and Predicted TSS values")
abline(lm(TSS~predict,data=comparisonANN))

h2o.varimp_plot(TSSnet2) #Variable importance plots
#####################################################################################################################################
#Model analysis

#I went through the package and found every function that extracts the features of 
#ANNwork models and copied them to a txt file ("featureofNet.txt")
h2o.mae(TSSnet2)#Mean absolute error
h2o.mean_residual_deviance(TSSnet2)#Mean residual deviance
h2o.r2(TSSnet2)#Training r^2
h2o.rmse(TSSnet2)#Root mean squared error value
h2o.rmsle(TSSnet2)#Root mean squared log error value
h2o.scoreHistory(TSSnet2) #Scoring history
#More Data Analysis

#####################################################################################################################################
#Print basic info about the ANNwork
#Extract feature importance
layer_Features<-h2o.deepfeatures(TSSnet2, total_data, layer = 0)
write.csv(as.data.frame(layer_Features),file="E:\\Marvin\\Research 2016\\ANN\\featureImportance ANN.csv")
layer_Features<-h2o.deepfeatures(TSSnet2, total_data, layer = 1)
write.csv(as.data.frame(layer_Features),file="E:\\Marvin\\Research 2016\\ANN\\featureImportance1 ANN.csv")
layer_Features<-h2o.deepfeatures(TSSnet2, total_data, layer = 2)
write.csv(as.data.frame(layer_Features),file="E:\\Marvin\\Research 2016\\ANN\\featureImportance2 ANN.csv")
layer_Features<-h2o.deepfeatures(TSSnet2, total_data, layer = 3)
write.csv(as.data.frame(layer_Features),file="E:\\Marvin\\Research 2016\\ANN\\featureImportance3 ANN.csv")

#Save Variable Importances
write.csv(as.data.frame(h2o.varimp(TSSnet2)),file="E:\\Marvin\\Research 2016\\ANN\\varimp ANN.csv")

#Save the weights and biases
weights<-as.data.frame(h2o.weights(TSSnet2))
biases<-as.data.frame(h2o.biases(TSSnet2))
write.csv(biases,file="E:\\Marvin\\Research 2016\\ANN\\Biases2 ANN.csv")
write.csv(weights,file="E:\\Marvin\\Research 2016\\ANN\\Weights2 ANN.csv")
write.csv(comparisonANN,file="E:\\Marvin\\Research 2016\\ANN\\comparisonANN.csv")
#####################################################################################################################################
summary(TSSnet2) #Just some summarizing info
h2o.shutdown(prompt=FALSE)
#Shut down are return cpus back to their home
