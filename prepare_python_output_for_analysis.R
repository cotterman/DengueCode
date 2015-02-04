###############################################################################
#################### prepare python output for analysis #######################
###############################################################################
# input: abundance data created in python
# output: abundance data reformated to work with prediction code and combined with clinical data 


prepare_python_output = function(respD, lcms_run, sample_name){
  
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
  
  # merge LCMS with clinical (no imputations)
  load(paste(outputsDir,"clinical_restricted_clean.RData", sep="")) 
  # using the restricted clinical data will correspond to dropping patients 86 and 1208 (the ones Natalia says are bad)
  comboD <- merge(clinical_restricted_clean, respD, by=c("code","Study"), all=F) 
  # keep only 1 observation for patient 251 (todo)
    #eliminating the ID0251 observation with MZ_1 value of 755341 resulted in better prediction than eliminating the one with MZ_1=773792
    #Natalia suggested going with the one with better mass spec data but it was not clear by observation which it should be -- see 9-7-2014 email.
  comboD = comboD[which(comboD$code!="ID0251" | comboD$MZ_1!=755341),]
  
  # merge LCMS with clinical data including imputations
  load(paste(outputsDir,"clin_full_wImputedRF1.RData", sep="")) 
  comboD_wImpRF1 <- merge(clin_full_wImputedRF1, respD, by=c("code","Study"), all=F) 
  # drop patients 86, 1208, and one of the 251 observations
  comboD_wImpRF1 = comboD_wImpRF1[which(comboD_wImpRF1$code!="ID0086" & 
                                        comboD_wImpRF1$code!="ID1208" &
                                       (comboD_wImpRF1$code!="ID0251" | comboD_wImpRF1$MZ_1!=755341)),] 
  
  
  return( list(comboD, comboD_wImpRF1) )

}

load(paste(outputsDir,"df1_from_python.RData", sep="")) #will load respD1_bins50x50 into R 
comboD1_bins50x50 = prepare_python_output(respD1_bins50x50, lcms_run=1)[[1]] #88 obs 
save(comboD1_bins50x50, file=paste(outputsDir,"df1_from_python_withRdata.RData", sep=""))
comboD1_bins50x50_wImpRF1 = prepare_python_output(respD1_bins50x50, lcms_run=1)[[2]] #88 obs
save(comboD1_bins50x50_wImpRF1, file=paste(outputsDir,"df1_from_python_wImpRF1.RData", sep=""))


load(paste(outputsDir,"df2_from_python.RData", sep="")) #will load respD2_bins50x50 into R 
comboD2_bins50x50 = prepare_python_output(respD2_bins50x50, lcms_run=2)[[1]]
save(comboD2_bins50x50, file=paste(outputsDir,"df2_from_python_withRdata.RData", sep=""))
comboD2_bins50x50_wImpRF1 = prepare_python_output(respD2_bins50x50, lcms_run=2)[[2]] 
save(comboD2_bins50x50_wImpRF1, file=paste(outputsDir,"df2_from_python_wImpRF1.RData", sep=""))



#combine D1 with D2 data
#list the columns that are in common between serum runs 1 and 2
commonvars = intersect(colnames(comboD1_bins50x50),colnames(comboD2_bins50x50)) 
comboD1D2_bins50x50 = rbind(comboD1_bins50x50[,commonvars], comboD2_bins50x50[,commonvars])

#make sure variable types are as expected
print("Variable types right before run_predictions function")
print( sapply(comboD1_bins50x50[, 1:65] , class) )



