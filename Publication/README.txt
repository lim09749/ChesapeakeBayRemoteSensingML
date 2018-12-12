##This file explains the purpose of each of the files in the ./Publication/ folder.

################################################Section 1####################################################################
##This section explains the purpose of the .R files (the code) in the folder
###############################################################################################################################
##Note that the machine learning algorithms occupied too much disk space, so they could not be stored here. Here is the google drive 
##link to the machine learning algorithms: https://drive.google.com/drive/folders/1fHipwjivyzlToDPpKg7R6Ttog5E7rB3s?usp=sharing

##Data assemblage
"dist_shore_dist_satellite.R"
##Measures effect of distance to shore from satellite and distance from satellite to in-situ on r-squared between reflectance and in-situ
"extract and clean data.R"
##Cleans data (removes points too close to the shore or too far away from in-situ measurement)

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
###############################################################################################################################
################################################Section 2####################################################################
##This section explains the purpose of the .grd and .gri files in the folder
###############################################################################################################################
##These files act as the monthly climatology data. Note that all of the machine learning algorithms are stacked for each month.
##The first part indicates which water quality constituent and the second part indicates the month.

##Ex: "TSS_Jan.grd" contains the monthy climatology of the all of the TSS machine learning algorithms predictions for January. 

###############################################################################################################################
################################################Section 3####################################################################
##This section explains the purpose of some of the remaining miscellaneous files
###############################################################################################################################

##Raw in-situ and satellite reflectance data
"TSS all stations Jul2002 to Dec2016 cleaned with MODISA.txt" 
"TSS and chl all stations Jul2002 to Dec2016 cleaned with MODISA.txt"

##Cleaned in-situ and satellite reflectance data
"CB_chl.csv"
"CB_tss.csv"

##Machine learning predicitons on test and entire dataset
"chl_predict_all.csv"
"chl_predict_test.csv"
"TSS_predict_all.csv"
"TSS_predict_test.csv"

##Metrics on test data
"chl_model_metrics.csv"
"TSS_model_metrics.csv"

##Used for monthly climatology (color)
"modis pretty.RData"

##Averages in the lower, middle, and upper Bay. 
"CHL_lower_avg.csv"
"CHL_middle_avg.csv"
"CHL_upper_avg.csv"
"TSS_lower_avg.csv"
"TSS_middle_avg.csv"
"TSS_upper_avg.csv"

##Line plots of monthly climatology (Not used in this study)
"TSS_line_plots_entire.tiff"
"chl_line_plots_entire.tiff"

##Plots of histogram, spectrum, stations
"histograms.tiff"
"spectrum.tiff"
"stations.tiff"
