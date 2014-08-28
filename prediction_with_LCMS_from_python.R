
#############################################################################
## wrapper for using LCMS data processed in Python in prediction function  ##
#############################################################################

run_predictions_wrap = function(respD, lcms_run){
  """respD should be a dataframe -- i.e., class = data.frame"""
    
  #### Specify location to place output ####
  sink(paste(resultsDir,"prediction_output_DEN_D2_bins.txt",sep=""), append=FALSE, split=TRUE)
  
  #### Esablish functions to be run #### 
  source("/home/carolyn/dengue_dx/Dengue_code/clean_data_functions.R")
  source("/home/carolyn/dengue_dx/Dengue_code/prediction_functions.R")
  
  #### Establish directory locations ####
  inputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_inputs/"
  outputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_outputs/"
  resultsDir = "/home/carolyn/dengue_dx/R_results/"
  
  #### add an indicator for type of study (will need this later when we merge) ####
  if(lcms_run==1){
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run)
  } 
  
  #### load other objects needed for prediction function (created by main_code.R) ####
  clinic_varsD = read.delim(paste(inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
  load(paste(outputsDir,"clinical_comboD_clean.RData")) 
  
  #### merge LCMS data with clinical data ####
  comboD2_filter50n <- merge(clinical_comboD_clean, respD, by=c("code","Study"), all=F) #75 obs
  
  #run prediction algorithms and save results
  resultsD2_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, "serum, D2", reduce_in_CV=T)
  #selected_D1_MFs_DEN = resultsD1_DEN[[2]]
  resultsD2_DEN[[1]]
  write.csv(x=resultsD2_DEN[[1]], file=paste(resultsDir,"resultsD2_DEN_bins.txt"), row.names = FALSE)

}