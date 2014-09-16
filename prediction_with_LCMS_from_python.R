
#############################################################################
## wrapper for using LCMS data processed in Python in prediction function  ##
#############################################################################

run_predictions_wrap = function(respD, lcms_run, load_from_R){
  
  print("Running R code using data matrix created in Python")
  
  library(lattice)
  library(SuperLearner)
  library(gtools) #enables smartbind
  #library(leaps) #for best subset selection
  library(bestglm) #for best subset selection
  library(MESS) #for area under the curve
  library(randomForest)
  library(glmnet) #lasso
  
  #### Establish directories ####
  
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
  
  #save(respD, file=paste(outputsDir,"df_from_python.RData", sep=""))
  load(paste(outputsDir,"df_from_python.RData", sep="")) #will load respD into R

      
  #### Specify location to place output ####
  sink(paste(resultsDir,"prediction_output_test2.txt", sep=""), append=FALSE, split=TRUE)
  
  #### Esablish functions to be run #### 
  source(paste(codeDir, "clean_data_functions.R", sep=""))
  source(paste(codeDir, "prediction_functions.R", sep=""))
  
  
  #### get variable naming convention compatible to what prediction code expects ####
  #LCMS data columns should be labeled with prefix "MZ_", and patient ID should be called "code"
  replace_X_with_MZ = function(x) sub("X", "MZ_", x)
  MZ_Names = sapply(X=colnames(respD), replace_X_with_MZ )
  newnames = as.character(c("code", MZ_Names[-c(length(MZ_Names))]))
  colnames(respD) = newnames #rename columns
  
  #convert intensity values from factors to numeric (not sure why they are ever factors...)
  Iindices = grep("MZ_",colnames(respD)) #column indices in dataset
  intensity_vars = colnames(respD[Iindices]) #names
  #convert these variables to numerics
  respD = convert_factor_to_numeric(respD, intensity_vars)

  #### add an indicator for type of study (will need this later when we merge) ####
  if(lcms_run==1){
    print("LCMS data batch 1")
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run)
  } 
  
  
  #### load other objects needed for prediction function (created by main_code.R) ####
  clinic_varsD = read.delim(paste(clinical_inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
  load(paste(outputsDir,"clinical_comboD_clean.RData", sep="")) 
  
  #### merge LCMS data with clinical data ####
  comboD1_filter50n <- merge(clinical_comboD_clean, respD, by=c("code","Study"), all=F) #75 obs
  save(comboD1_filter50n, file=paste(outputsDir,"df_from_python_withRdata.RData", sep=""))
  #load(paste(outputsDir,"df_from_python_withRdata.RData", sep="")) #will load comboD1_filter50n into R
  
  #make sure variable types are as expected
  print("Variable types right before run_predictions function")
  print( sapply(comboD1_filter50n[, 1:65] , class) )
  
  #run prediction algorithms and save results
  resultsD1_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, "serum, D1", reduce_in_CV=T)
  #selected_D1_MFs_DEN = resultsD1_DEN[[2]]
  resultsD1_DEN[[1]]
  write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"resultsD1_DEN_bins.txt", sep=""), row.names = FALSE)
  
}
