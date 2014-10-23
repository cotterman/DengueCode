
#Run with R 3.1.1
#To run remotely using Rstudio, go to URL: amold.lbl.gov:8787/

################################################################################
############### Prepare and describe data for Dengue Dx project ################
################################################################################

rm(list = ls()) #start with blank slate (clear everything from workspace)

library(lattice)
library(SuperLearner)
library(gtools) #enables smartbind
#library(leaps) #for best subset selection
library(bestglm) #for best subset selection
library(MESS) #for area under the curve
library(randomForest)
library(glmnet) #lasso
library(cvAUC) 
#library(xlsReadWrite) #need newer version of R for this package (NA for 3.1.1 also)

###############################################################################
########################## Establish directories ##############################
###############################################################################

#directory containing code
codeDir = "/srv/scratch/carolyn/Dengue_code/" #on Amold from Amold's perspective
#codeDir = "/home/carolyn/temp_Dengue_code/" #on myPC

#select main directory in which to find subfolders
homeDir = "/srv/scratch/carolyn/" #on Amold from Amold's perspective
#homeDir = "/home/carolyn/dengue_dx/" #on my PC

clinical_inputsDir = paste(homeDir, "lab_and_clinical_data/Cleaned/", sep="")
lcmsCO_inputsDir = paste(homeDir, "lcms_data_processed_in_CO/Cleaned/", sep="")
outputsDir = paste(homeDir, "intermediate_data/", sep="")
resultsDir = paste(homeDir, "Results/", sep="")

#directory containing code
codeDir = "/srv/scratch/carolyn/Dengue_code/" #on Amold from Amold's perspectiev
#codeDir = "/home/carolyn/temp_Dengue_code/" #on myPC


###############################################################################
########################### Load data to be used ##############################
###############################################################################

#info on clinical variables (data types and applicability to various analyses)
clinic_varsD = read.delim(paste(clinical_inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)

#this data was created in "create_data_for_analysis.R"
load(paste(outputsDir,"comboD1_filter50n.RData", sep="")) #loads comboD1_filter50n

#this data was created in "prepare_python_for_analysis.R"
load(paste(outputsDir,"df1_from_python_withRdata.RData", sep="")) #loads comboD1_bins50x50
load(paste(outputsDir,"df2_from_python_withRdata.RData", sep="")) #loads comboD2_bins50x50


###############################################################################
############################ Predict diagnosis ################################
###############################################################################

#### Specify location to place output ####
sink(paste(resultsDir,"prediction_output_DEN_v7.txt",sep=""), append=FALSE, split=TRUE)

#### Esablish functions to be run #### 
source(paste(codeDir,"prediction_functions.R",sep=""))

### Serum Samples from Nicaragua ###

resultsD1_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD1_bins50x50, "serum, D1", reduce_in_CV=T)
resultsD1_DEN[[1]] #view results
write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"resultsD1_DEN_bins_v7.txt",sep=""), row.names = FALSE)

resultsD1_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, "serum, D1", reduce_in_CV=T)
resultsD1_DEN[[1]] #view results
write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"resultsD1_DEN_CO_v7.txt",sep=""), row.names = FALSE)

resultsD2_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, "serum, D2", reduce_in_CV=T)
resultsD2_DEN[[1]] #view results
write.csv(x=resultsD2_DEN[[1]], file=paste(resultsDir,"resultsD2_DEN_v7.txt",sep=""), row.names = FALSE)


resultsD1_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD1_filter50n, "serum, D1", reduce_in_CV=T)
#selected_D1_MFs_SDEN = resultsD1_SDEN[[2]]
write.csv(x=resultsD1_SDEN[[1]], file=paste(resultsDir,"resultsD1_SDEN_v6a.txt"), row.names = FALSE)

resultsD2_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD2_filter50n, "serum, D2", reduce_in_CV=T)
selected_D2_MFs_SDEN = resultsD2_SDEN[[2]]
write.csv(x=resultsD2_SDEN[[1]], file=paste(resultsDir,"resultsD2_SDEN_v6a.txt"), row.names = FALSE)

## combine lists of MFs selected by dimension reduction methods
flist_DEN = merge(selected_D1_MFs_DEN, selected_D2_MFs_DEN,  by="Var1", all=TRUE)
flist_DEN
write.csv(x=flist_DEN, file=paste(resultsDir,"MFs_selected_by_RF_DEN_v6a.txt"), row.names = FALSE)
flist_SDEN = merge(selected_D1_MFs_SDEN, selected_D2_MFs_SDEN,  by="Var1", all=TRUE)
flist_SDEN
write.csv(x=flist_SDEN, file=paste(resultsDir,"MFs_selected_by_RF_SDEN_v6a.txt"), row.names = FALSE)


### Non-invasive samples from Nicaragua -- ND vs DEN only ###

resultsD3_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, "saliva, D3", reduce_in_CV=T)
#selected_D1_MFs_DEN = resultsD1_DEN[[2]]
resultsD3_DEN[[1]]
write.csv(x=resultsD3_DEN[[1]], file=paste(resultsDir,"resultsD3_DEN_v7.txt"), row.names = FALSE)

resultsD5_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, "urine, D5", reduce_in_CV=T)
#selected_D1_MFs_DEN = resultsD1_DEN[[2]]
resultsD5_DEN[[1]]
write.csv(x=resultsD5_DEN[[1]], file=paste(resultsDir,"resultsD5_DEN_v7.txt"), row.names = FALSE)



###############################################################################
#################### Variable Importance Analysis without CV ##################
###############################################################################

#### Esablish functions to be run #### 
source(paste(codeDir, "prediction_functions.R",sep=""))


D1_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, dim_reduce_covar=F, "D1_covarF")
D1_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, dim_reduce_covar=T, "D1_covarT")
D1_VIM = merge(D1_VIM_covarT, D1_VIM_covarF, by="row.names", all=T)

D2_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, dim_reduce_covar=F, "D2_covarF")
D2_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, dim_reduce_covar=T, "D2_covarT")
D2_VIM = merge(D2_VIM_covarT, D2_VIM_covarF, by="row.names", all=T)

D3_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, dim_reduce_covar=F, "D3_covarF")
D3_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, dim_reduce_covar=T, "D3_covarT")
D3_VIM = merge(D3_VIM_covarT, D3_VIM_covarF, by="row.names", all=T)

D5_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, dim_reduce_covar=F, "D5_covarF")
D5_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, dim_reduce_covar=T, "D5_covarT")
D5_VIM = merge(D5_VIM_covarT, D5_VIM_covarF, by="row.names", all=T)

## combine all variable importance lists ##
serum_VIM = merge(D1_VIM, D2_VIM, by="Row.names", all=T)
noninvasive_VIM = merge(D3_VIM, D5_VIM, by="Row.names", all=T)
VIM_RF = merge(serum_VIM, noninvasive_VIM, by="Row.names", all=T)
write.csv(x=VIM_RF, file=paste(resultsDir,"VIM_RF_DEN.txt"), row.names = FALSE)

