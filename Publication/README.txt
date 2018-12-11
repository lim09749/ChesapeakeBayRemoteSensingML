##This file explains the purpose of the .R files (the code) in the folder
##Note that the machine learning algorithms occupied too much disk space, so they could not be stored here. Here is the google drive 
##link to the machine learning algorithms: https://drive.google.com/drive/folders/1fHipwjivyzlToDPpKg7R6Ttog5E7rB3s?usp=sharing
##Data assemblage
"dist_shore_dist_satellite.R"
##Measures effect of distance to shore from satellite and distance from satellite to in-situ on r-squared between reflectance and in-situ
"extract and clean data.R"

##Basic plots
"histograms.R"
##Histogram of in-situ chlorophyll and TSS concentrations
"spectrum plots.R"
##Spectrum plots of satellite reflectance
"stations.R"
##Map of station locations


##Models to Create Machine learning algorithms (ANN, GLM, RVM, SVM)
##Total Suspended Sediment
"TSS SVM cv with parallel.R"
"TSS RVM cv with parallel.R"
"TSS ANN cv with parallel.R"
"TSS GLM cv with parallel.R"

##Chlorophyll
"CHL SVM cv with parallel.R"
"CHL RVM cv with parallel.R"
"CHL ANN cv with parallel.R"
"CHL GLM cv with parallel.R"

##Scatterplots and metrics
"CHL scatterplots and metrics.R"
"TSS scatterplots and metrics.R"
"CHL scatterplots and metrics without running models.R"
"TSS scatterplots and metrics without running models.R"

##Line plots of monthly averages
"CHL line plots mc without running models.R"
"CHL line plots mc.R"
"TSS line plots mc without running models.R"
"TSS line plots mc.R"

##Monthly climatology images
"CHL mc.R"
"CHL mc without running model.R"
"TSS mc.R"
"TSS mc without running model.R"


