
#Run with R 3.1.1
#To run remotely using Rstudio, go to URL: amold.lbl.gov:8787/

################################################################################
############### Prepare and describe data for Dengue Dx project ################
################################################################################

rm(list = ls()) #start with blank slate (clear everything from workspace)

library(lattice)
library(ggplot2)
library(reshape2)
library(sjPlot)
#library(SuperLearner) #the CRAN version will sometimes throw errors
#library(devtools)
#install_github("ecpolley/SuperLearner") #this version corrects bug in CRAN version (last downloaded on 1-6-2015)
library(SuperLearner) #will load whichever SuperLearner was last installed
library(gtools) #enables smartbind
#library(leaps) #for best subset selection
library(bestglm) #for best subset selection
library(MESS) #for area under the curve
library(randomForest)
library(glmnet) #lasso
library(e1071) #svm
library(rpart) #Recursive partitioning for classification, regression and survival trees
library(party) #classification tree algorithm
library(caret) #classification package that runs with SL package
library(survival)
library(LogicReg) #Logic regression (to find predictors that are boolean combos of original predictors)
library(cvAUC) 
library(Hmisc) #allows labeling of data.frame columns (for documentation's sake), and has useful function contents()
library(xtable)
library(parallel) #for super learner to make use of multiple cores
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
codeDir = "/srv/scratch/carolyn/Dengue_code/" #on Amold from Amold's perspective
#codeDir = "/home/carolyn/temp_Dengue_code/" #on myPC

###############################################################################
##################### Esablish functions to be run ############################

source(paste(codeDir,"clean_data_functions.R",sep=""))
source(paste(codeDir,"clean_clinical_data_functions_v2.R",sep=""))
source(paste(codeDir,"summarize_clinical_functions.R",sep=""))


###############################################################################
########################### Create data to be used ############################
###############################################################################

#info on clinical variables (data types and applicability to various analyses)
clinic_varsD = read.delim(paste(clinical_inputsDir,"List of clinical variables for analysis_v5.txt", sep=""), header=TRUE, nrows=500)

#this takes a bit of time to run (~ 10 min) b/c of the imputing missing values step
#source(paste(codeDir,"create_data_for_analysis.R",sep="")) 


###############################################################################
########################### Load data to be used ##############################
###############################################################################

# this data was created in "create_data_for_analysis.R" - contains all clinical data (n=1624)
load(paste(outputsDir,"clinical_full_clean.RData", sep="")) #loads clinical_full_clean
load(paste(outputsDir,"clinical_D1_clean.RData", sep="")) #loads clinical_D1_clean

# this data was created in "create_data_for_analysis.R" - contains all clinical data with imputed values
load(paste(outputsDir,"clin_full_wImputedRF1.RData", sep="")) #loads clin_full_wImputedRF1

# this data was created in "create_data_for_analysis.R" - contains mass hunter LCMS combined with clinical
load(paste(outputsDir,"comboD1_filter50n.RData", sep="")) #loads comboD1_filter50n
#load(paste(outputsDir,"comboD3_filter50n.RData", sep="")) #loads comboD3_filter50n
#load(paste(outputsDir,"comboD5_filter50n.RData", sep="")) #loads comboD5_filter50n
# this data was created in "create_data_for_analysis.R" - contains mass hunter LCMS combined with imputed clinical
load(paste(outputsDir,"comboD1_filter50n_wImpRF1.RData", sep="")) #loads comboD1_filter50n_wImpRF1
#load(paste(outputsDir,"comboD3_filter50n_wImpRF1.RData", sep="")) #loads comboD3_filter50n_wImpRF1
#load(paste(outputsDir,"comboD5_filter50n_wImpRF1.RData", sep="")) #loads comboD5_filter50n_wImpRF1

# this data was created in "prepare_python_output_for_analysis.R" - contains binned LCMS combined with clinical
#load(paste(outputsDir,"df1_from_python_withRdata.RData", sep="")) #loads comboD1_bins50x50
load(paste(outputsDir,"df1_from_python_wImpRF1.RData", sep="")) #loads comboD1_bins50x50_wImpRF1
#load(paste(outputsDir,"df2_from_python_withRdata.RData", sep="")) #loads comboD2_bins50x50


###############################################################################
############################ Summarize data ###################################
###############################################################################

#source(paste(codeDir,"summarize_clinical_data.R",sep=""))


###############################################################################
############################ Predict diagnosis ################################
###############################################################################

#### Esablish functions to be run #### 
source(paste(codeDir,"prediction_functions_v2.R",sep=""))

#### Specify location to place output ####
#sink(paste(resultsDir,"predictions_clinical_DEN_v9.txt",sep=""), append=FALSE, split=TRUE)

#get list of clinical variables to include in prediction method
covarlist_all = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
covarlist_CohortRestrict = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                    eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                    XD=clinical_full_clean, restrict_to_cohort_vars=T, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
covarlist_noUltraX = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                         eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                         XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=F, BloodLab=T)
covarlist_genOnly = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                        eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                        XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=F, BloodLab=F)
covarlist_noMiss = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=T, 
                                        eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                        XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=F, UltraX=T, BloodLab=T)
#These are the 41 clinical vars included in original (R33) analysis -- the ones that are never missing among our 88 samples 
covarDEN_88noMiss = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=T, 
                                       eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                       XD=clinical_D1_clean, restrict_to_cohort_vars=T, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)


############### Predictions using clinical info only ##########################

#test
#testme = run_predictions(clinic_varsD, covarlist_noUltraX, "ND.vs.DEN", comboD1_filter50n_wImpRF1, paste("D1, noUltraX clinical"), 
#                         include_imp_dums="all")

#options(error=browser()) #view problem

run_runs = function(pdata, pdata_desc, imp){
  #pdata is data to use, pdata_desc is string that descibes data, imp is "all","study" or "none" for inclusion of imputation indicators
  #all other parameters are from namespace
  
  ##### ND vs DEN #####
  DEN_allClin = run_predictions(clinic_varsD, covarlist_all, "ND.vs.DEN", pdata, paste(pdata_desc,", all clinical"), include_imp_dums=imp)
  DEN_cohortRclin = run_predictions(clinic_varsD, covarlist_CohortRestrict, "ND.vs.DEN", pdata, paste(pdata_desc,", CohortRestrict clinical"), include_imp_dums=imp)
  #DEN_noUXclin = run_predictions(clinic_varsD, covarlist_noUltraX, "ND.vs.DEN", pdata, paste(pdata_desc,", noUltraX clinical"), include_imp_dums=imp)
  #DEN_genOnlyclin= run_predictions(clinic_varsD, covarlist_genOnly, "ND.vs.DEN", pdata, paste(pdata_desc,", genOnly clinical"), include_imp_dums=imp)
  #DEN_noMiss88clin= run_predictions(clinic_varsD, covarDEN_88noMiss, "ND.vs.DEN", pdata, paste(pdata_desc,", noMiss88 clinical"), include_imp_dums=imp)
  ##### DF vs DHF/DSS #####
  #DHF_allClin = run_predictions(clinic_varsD, covarlist_all, "DF.vs.DHF.DSS", pdata, paste(pdata_desc,", all clinical"), include_imp_dums=imp)
  #DHF_noUXclin = run_predictions(clinic_varsD, covarlist_noUltraX, "DF.vs.DHF.DSS", pdata, paste(pdata_desc,", noUltraX clinical"), include_imp_dums=imp)
  #DHF_cohortRclin = run_predictions(clinic_varsD, covarlist_CohortRestrict, "DF.vs.DHF.DSS", pdata, paste(pdata_desc, ", CohortRestrict clinical"), include_imp_dums=imp)
  #DHF_genOnlyclin= run_predictions(clinic_varsD, covarlist_genOnly, "DF.vs.DHF.DSS", pdata, paste(pdata_desc,", genOnly clinical"), include_imp_dums=imp)
  #DHF_noMissclin= run_predictions(clinic_varsD, covarlist_noMiss, "DF.vs.DHF.DSS", pdata, paste(pdata_desc,", noMiss clinical"), include_imp_dums=imp)
  
  #combine results and output to file
  #myruns = rbind(DEN_allClin[[1]], DEN_cohortRclin[[1]], DEN_noUXclin[[1]], DEN_genOnlyclin[[1]], DEN_noMissclin[[1]],
  #      DHF_allClin[[1]], DHF_cohortRclin[[1]], DHF_noUXclin[[1]], DHF_genOnlyclin[[1]], DHF_noMissclin[[1]])
  #myruns = rbind(DEN_allClin[[1]], DEN_cohortRclin[[1]],  DEN_noUXclin[[1]], DEN_genOnlyclin[[1]], DEN_noMiss88clin[[1]])
  myruns = rbind(DEN_allClin[[1]], DEN_cohortRclin[[1]])
  return(myruns)
}

Run_n1624_clinOnly=F  #need to re-rerun with imput dums at some point
if(Run_n1624_clinOnly==T){
  #pred with full data -- include imputed values and also indicators of imputation (need to rerun, though results look very similar to results from above)
  all_wImputdums_clinOnly = run_runs(clin_full_wImputedRF1, pdata_desc="all data",imp="all")
  write.csv(x=all_wImputdums_clinOnly, file=paste(resultsDir,"all_wImputdums_clinOnly.txt",sep=""), row.names = FALSE)
  #pred with full data --- include imputed values with no imputation indicators
  all_wImput_clinOnly = run_runs(clin_full_wImputedRF1, pdata_desc="all data", imp="none")
  write.csv(x=all_wImput_clinOnly, file=paste(resultsDir,"all_wImput_clinOnly.txt",sep=""), row.names = FALSE)
  #pred with full data --- include imputed values and indicator of study type (cohort), but no imputation indicators
  all_wImputstudy_clinOnly = run_runs(clin_full_wImputedRF1, pdata_desc="all data", imp="study")
  write.csv(x=all_wImputstudy_clinOnly, file=paste(resultsDir,"all_wImputstudy_clinOnly.txt",sep=""), row.names = FALSE)
}

Run_D1_clinOnly=F #finished running this (without age or DaysSick)
if(Run_D1_clinOnly==T){
  #create D1 data that has imputed info
  clin_D1_wImputedRF1 = clin_full_wImputedRF1[which(clin_full_wImputedRF1$serum==1),]
  #pred with full data --- include imputed values with no imputation indicators
  D1_wImput_clinOnly = run_runs(clin_D1_wImputedRF1, pdata_desc="D1", imp="none")
  write.csv(x=D1_wImput_clinOnly, file=paste(resultsDir,"D1_wImput_clinOnly.txt",sep=""), row.names = FALSE)
  #pred with D1 data --- include imputed values and indicator of study type (cohort), but no imputation indicators
  D1_wImputstudy_clinOnly = run_runs(clin_D1_wImputedRF1, pdata_desc="D1", imp="study")
  write.csv(x=D1_wImputstudy_clinOnly, file=paste(resultsDir,"D1_wImputstudy_clinOnly.txt",sep=""), row.names = FALSE)
  #pred with D1 data -- include imputed values and also indicators of imputation
  D1_wImputdums_clinOnly = run_runs(clin_D1_wImputedRF1, pdata_desc="D1", imp="all") #trimLogit error for DEN w noUltraX and w genOnly
  write.csv(x=D1_wImputdums_clinOnly, file=paste(resultsDir,"D1_wImputdums_clinOnly.txt",sep=""), row.names = FALSE)
}


############### Predictions using clinical and LCMS ###########################

## Results using just sample of 88.
# mass hunter
resultsD1_DEN = run_predictions(clinic_varsD, covarlist_all, "ND.vs.DEN", comboD1_filter50n_wImpRF1, "CohortRestrict, D1", include_imp_dums="none")
#write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"D1_wImput_LCMS_selectionMethods_MH.txt",sep=""), row.names = FALSE)
# binned
resultsD1_DEN_bins = run_predictions(clinic_varsD, covarlist_CohortRestrict, "ND.vs.DEN", comboD1_bins50x50_wImpRF1, "CohortRestrict, D1", include_imp_dums="none")
write.csv(x=resultsD1_DEN_bins[[1]], file=paste(resultsDir,"D1_wImput_LCMS_bins88b_selection.txt",sep=""), row.names = FALSE)
#serum batch 2 (n=75)
clin_D2_wImputedRF = clin_full_wImputedRF1[which(clin_full_wImputedRF1$serum==2),]
resultsD1_DEN = run_predictions(clinic_varsD, covarlist_all, "ND.vs.DEN", clin_D2_wImputedRF, "CohortRestrict, D2", include_imp_dums="none")

if(T==F){
  
### Non-invasive samples from Nicaragua -- ND vs DEN only ###

resultsD3_DEN = run_predictions(clinic_varsD, covarlist_CohortRestrict, "ND.vs.DEN", comboD3_filter50n_wImpRF1, "CohortRestrict, saliva, D3", include_imp_dums="none")
resultsD3_DEN[[1]]
write.csv(x=resultsD3_DEN[[1]], file=paste(resultsDir,"resultsD3_DEN_CohortRestrict.txt"), row.names = FALSE)

resultsD5_DEN = run_predictions(clinic_varsD, covarlist_CohortRestrict, "ND.vs.DEN", comboD5_filter50n_wImpRF1, "CohortRestrict, urine, D5", include_imp_dums="none")
resultsD5_DEN[[1]]
write.csv(x=resultsD5_DEN[[1]], file=paste(resultsDir,"resultsD5_DEN_CohortRestrict.txt"), row.names = FALSE)


### Serum Samples from Nicaragua ###

resultsD1_DEN = run_predictions(clinic_varsD, covarlist, "ND.vs.DEN", comboD1_bins50x50, "serum, D1")
resultsD1_DEN[[1]] #view results
write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir, "resultsD1_DEN_bins_v7.txt",sep=""), row.names = FALSE)

resultsD1_DEN = run_predictions(clinic_varsD, covarDEN_88noMiss, "ND.vs.DEN", comboD1_filter50n_wImpRF1, "noMiss88, serum, D1", include_imp_dums="none")
resultsD1_DEN[[1]] #view results
write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"resultsD1_DEN_CO_v9_last.txt",sep=""), row.names = FALSE)

resultsD2_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, "serum, D2")
resultsD2_DEN[[1]] #view results
write.csv(x=resultsD2_DEN[[1]], file=paste(resultsDir,"resultsD2_DEN_v7.txt",sep=""), row.names = FALSE)


resultsD1_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD1_filter50n, "serum, D1")
#selected_D1_MFs_SDEN = resultsD1_SDEN[[2]]
write.csv(x=resultsD1_SDEN[[1]], file=paste(resultsDir,"resultsD1_SDEN_v6a.txt"), row.names = FALSE)

resultsD2_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD2_filter50n, "serum, D2")
selected_D2_MFs_SDEN = resultsD2_SDEN[[2]]
write.csv(x=resultsD2_SDEN[[1]], file=paste(resultsDir,"resultsD2_SDEN_v6a.txt"), row.names = FALSE)

## combine lists of MFs selected by dimension reduction methods
flist_DEN = merge(selected_D1_MFs_DEN, selected_D2_MFs_DEN,  by="Var1", all=TRUE)
flist_DEN
write.csv(x=flist_DEN, file=paste(resultsDir,"MFs_selected_by_RF_DEN_v6a.txt"), row.names = FALSE)
flist_SDEN = merge(selected_D1_MFs_SDEN, selected_D2_MFs_SDEN,  by="Var1", all=TRUE)
flist_SDEN
write.csv(x=flist_SDEN, file=paste(resultsDir,"MFs_selected_by_RF_SDEN_v6a.txt"), row.names = FALSE)


###############################################################################
#################### Variable Importance Analysis without CV ##################
###############################################################################

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
}

