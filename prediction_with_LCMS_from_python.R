
#############################################################################
## wrapper for using LCMS data processed in Python in prediction function  ##
#############################################################################

(clinic_varsD, outcome, mydata, sample_name, reduce_in_CV=T)
run_predictions_wrap = function(respD, lcms_run, sample_name){
  
  print("Running R code using data matrix created in Python")


  #We might want to drop a couple of the R1 observations, as Natalia suggests (9-7-2014):
     # I had to drop the files 86 and 1208 because in my chromatography image seemed to have high abundances of PEG contamination
     # For patient 251, .....use the file that has the higher abundance and discard the other.
  save(respD, file=paste(outputsDir,"df1_from_python.RData", sep=""))
  load(paste(outputsDir,"df1_from_python.RData", sep="")) #will load respD into R

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
  respD = convert_factor_to_numeric(respD, intensity_vars) #takes a min to run

  #### add an indicator for type of study (will need this later when we merge) ####
  if(lcms_run==1){
    print("LCMS data batch 1")
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run, clinical_inputsDir)
  } 
  
  
  #### load other objects needed for prediction function (created by main_code.R) ####
  clinic_varsD = read.delim(paste(clinical_inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
  load(paste(outputsDir,"clinical_comboD_clean.RData", sep="")) 
  
  #### merge LCMS data with clinical data ####
  comboD1_filter50n <- merge(clinical_comboD_clean, respD, by=c("code","Study"), all=F) #75 obs
  save(comboD1_filter50n, file=paste(outputsDir,"df1_from_python_withRdata.RData", sep=""))
  load(paste(outputsDir,"df2_from_python_withRdata.RData", sep="")) #will load comboD1_filter50n into R
  #sanity check: check out results when ND/DEN labels are randomly permuted
  #comboD1_filter50n$DEN_dum = permute(comboD1_filter50n$DEN_dum)
  
  #combine D1 with D2 data
  #list the columns that are in common between serum runs 1 and 2
  commonvars = intersect(colnames(comboD1_filter50n),colnames(comboD2_filter50n)) 
  comboD1D2_filter50n = rbind(comboD1_filter50n[,commonvars], comboD2_filter50n[,commonvars])

  #make sure variable types are as expected
  print("Variable types right before run_predictions function")
  print( sapply(comboD1_filter50n[, 1:65] , class) )
  
  #### Specify location to place output ####
  sink(paste(resultsDir,"prediction_output_bins_D1D2.txt", sep=""), append=FALSE, split=TRUE)
  
  #run prediction algorithms and save results
  resultsD1_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD1D2_filter50n, "serum, D1D2", reduce_in_CV=T)
  selected_D1_MFs_DEN = resultsD1_DEN[[2]]
  resultsD1_DEN[[1]]
  write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"results_DEN_bins_D1D2.txt", sep=""), row.names = FALSE)
  
}

