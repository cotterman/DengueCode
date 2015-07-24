
library(multiPIM)
library(tmle) #works only for binary treatments
library(tmle.npvi) #works for continuous treatments but must provide reference value for X (or assumed 0)
    #and I think it is estimating a linear fit for VIM, which is not cool
  #I could force continuous variables to be on scale from 0 to 1....should find out more about how reference value of X is used in tmle.
  #I could also turn to binary by giving value of 0/1 if less than/greater than 75th percentile (or something)
  #For now, I could just do the variable importance stuff for general signs & symptoms, which I have already dichotomized
    #or, forget about tMLE altogether



get_rfVIMs = function(outcome, include_imp_dums, NoOFI){
  
  ### Parameter options ###
  #outcome = "is.DEN", "is.DHF_DSS"
  #include_imp_dums = "all", "None", "Only"
  #NoOFI = TRUE, FALSE
  #test
  outcome="is.DEN"
  include_imp_dums="Only"
  NoOFI=FALSE

  ### Parameters based on choosen options ###
  if (outcome=="is.DEN") {
    FileNamePrefix = "OFI.v.DEN" #use all samples; exclude PCR predictors
    NoInitialDHF = FALSE #whether to exclude samples with initial DHF/DSS diagnosis
    NoOFI = FALSE
  }
  if (outcome=="is.DHF_DSS"){
    NoInitialDHF = TRUE #whether to exclude samples with initial DHF/DSS diagnosis
    if (NoOFI == TRUE) {
      FileNamePrefix = "DF.v.DHFDSS"
    }
    else{
      FileNamePrefix = "OFIDF.v.DHFDSS"
    }
  }
  ## Suffix to indicate use of imputation dummies
  if (include_imp_dums=="Only") {
    FileNameSuffix = "_dumsOnly"
  }
  if (include_imp_dums=="all") {
      FileNameSuffix = "_impDums"
  } 
  if (include_imp_dums=="None") {
    FileNameSuffix = ""  # all covariates, but no missing indicators
  }
  
  ### Obtain data and variable list ###
  
  #Use just hospital data
  myhospit = clin12_full_wImputedRF1[which(clin12_full_wImputedRF1$Study=="Hospital"),]
  #drop patients with initial DHF/DSS Dx if specified
  if (NoInitialDHF==TRUE){
    myhospit = myhospit[which(myhospit$WHO_initial_given!="DHF" & myhospit$WHO_initial_given!="DSS"),]
  }
  #drop the DEN negative patients if specified
  if (NoOFI == TRUE){
    myhospit = myhospit[which(myhospit$is.DEN==TRUE),]
  }
  
  #obtain list of clinical variables to include.
    #this function is created in clean_clinical_data_functions_v2.R
  clinic_vars = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                    eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=10,
                                    XD=myhospit, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
  
  #Here is an opportunity to add imputation dummies (include_imp_dum="all")
    #this function is created in "prediction_functions_v2.R"
    #It will also eliminate PCR for OFI.v.DEN analysis
      #But we will run in this "OFI.v.DEN" mode b/c to avoid dropping the DHF/DSS indicators that now we want to keep
  myvars = get_clinical_predictors("OFI.v.DEN", clinic_varsD, clinic_vars, myhospit, include_imp_dums=include_imp_dums) 
  
  #in python we were not allowed factors so we used the binary version of the categorical variables  
  #since we are comparing VIMs in python with VIMs here, make all variables equivalent
  if (include_imp_dums!="Only"){
    prelim = myvars[myvars %in% c("Torniquete","Pulso","Sexo")==FALSE] #remove vars
    newlist = c(prelim, "is.torniquete20plus", "is.pulse_danger","is.female") #add vars to replace them
    #create the is.female variable
    myhospit[which(myhospit$Sexo=="female"), "is.female"] = 1
    myhospit[which(myhospit$Sexo=="male"), "is.female"] = 0
    #serotype will only be used for the DF vs DHF/DSS analysis
    if (FileNamePrefix == "DF.v.DHFDSS"){
      newlist = c(newlist, "is.serotype1", "is.serotype2", "is.serotype3")
    }
  }
  else{
    newlist = myvars
  }
  print("Final var list on which to do variable importance")
  print(newlist)
  
  #make outcome a factor so RF will do classification
  myhospit[which(myhospit[,outcome]==TRUE), "DEN_dum"] = "1"
  myhospit[which(myhospit[,outcome]==FALSE), "DEN_dum"] = "0"
  myhospit$DEN_dum = factor(myhospit[,outcome])
  
  ### Run random forest and output VIMs ###
  rf = randomForest(x=myhospit[,c(newlist)], y=myhospit[,"DEN_dum"], ntree=500, nPerm=1, importance=TRUE)
  VIM_OOB = importance(rf, type=1) #"mean decrease accuracy" from permuting variable in oob data, divided by SD
  oob = as.data.frame(cbind(as.numeric(VIM_OOB), rownames(VIM_OOB) ))
  VIM_Gini = importance(rf, type=2) #mean decrease Gini 
  gini = as.data.frame(cbind(as.numeric(VIM_Gini), rownames(VIM_Gini) ))
  VIM_rf = merge(oob, gini, by="V2")
  colnames(VIM_rf) = c("varname","RF_OOB","RF_Gini")
  #sort -- grr.....values are interpreted as chars not numerics for sorting.  strange.
  VIM_rf = VIM_rf[order(VIM_rf$RF_OOB),]
  
  ## Add proper variable names for displaying in tables/graphs and variable categories
  #first, modify the usual clinic_varsD to reflect the binary vars that we are now using
  mod_clinic_varsD = clinic_varsD
  mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="Sexo"),"variable.name.in.final.data"]="is.female"
  mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="Pulso"),"variable.name.in.final.data"]="is.pulse_danger"
  mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="Torniquete"),"variable.name.in.final.data"]="is.torniquete20plus"
  #include 1 row for each value of serotype
  add1 = mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="PCR"),]
  add1$CC_name = "Serotype 1"
  add1$variable.name.in.final.data = "is.serotype1"
  add2 = mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="PCR"),]
  add2$CC_name = "Serotype 2"
  add2$variable.name.in.final.data = "is.serotype2"
  add3 = mod_clinic_varsD[which(mod_clinic_varsD$variable.name.in.final.data=="PCR"),]
  add3$CC_name = "Serotype 3"
  add3$variable.name.in.final.data = "is.serotype3"
  mod_clinic_varsD = rbind(mod_clinic_varsD, add1, add2, add3)
  VIM_rf_plus = merge.data.frame(mod_clinic_varsD[,c("variable.name.in.final.data","CC_name","CC_category","CC_type",
                                                     "CC_description","in.cohort.data","CC_broadcat_sort", "CC_cat_sort")], 
                                 VIM_rf, by.x="variable.name.in.final.data", by.y="varname", in.y=all)
  #now export so we can read it into python and make graphs with other VIMs
  VIMout_from_R = paste("VIM_rf_", FileNamePrefix, FileNameSuffix, sep="")
  write.table(VIM_rf_plus, paste(outputsDir, VIMout_from_R, ".txt", sep=""), sep="\t") 
}

#no imputation dummies
get_rfVIMs(outcome="is.DEN", include_imp_dums="None", NoOFI=FALSE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="None", NoOFI=TRUE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="None", NoOFI=FALSE)

## Imputation code needs to be fixed -- problem with adding info from clinic_varsD
  #since merge is done on variable names without _imputed suffix
#imputation dummies plus clinical feature values
get_rfVIMs(outcome="is.DEN", include_imp_dums="all", NoOFI=FALSE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="all", NoOFI=TRUE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="all", NoOFI=FALSE)
#imputation dummies only
get_rfVIMs(outcome="is.DEN", include_imp_dums="Only", NoOFI=FALSE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="Only", NoOFI=TRUE)
get_rfVIMs(outcome="is.DHF_DSS", include_imp_dums="Only", NoOFI=FALSE)

#check results - empty dataset
checkme = read.delim(paste(outputsDir,"VIM_rf_DF.v.DHFDSS_dumsOnly.txt", sep=""), header=TRUE, nrows=200)







###############################################################################
############################## old analyses ###################################
###############################################################################

get_VIM = function(XD, method, suffix, clinic_vars, rf_type){
  
  #impute missing values using RF proximity measure (see RandomForest documentation)
  set.seed(101)
  XD_missimput = rfImpute(x=XD[,c(clinic_vars)], y=XD[,"DEN_dum"], ntree=500)
  colnames(XD_missimput) = c("DEN_dum",clinic_vars)
  
  #RF variable importance
  if(method=="rfVIM"){
    set.seed(102)
    rf = randomForest(x=XD_missimput[,c(clinic_vars)], y=XD_missimput[,"DEN_dum"], ntree=500, localImp=FALSE)
    if(rf_type=="oob"){
      rf_VIPS = importance(rf, type=1) #type=1 uses the permutation method.  
    }
    if(rf_type=="gini"){
      rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  
    }
    rfVIP = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
    myVIM = as.data.frame(cbind(names(rfVIP), rfVIP))
    colnames(myVIM) = c("VarName",paste("rfVIM",suffix,sep="_"))
  }
  #convert variables into numerics
  YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                             "Variable.Name.in.Hospital.data"]
  YesNo_varsH_char = sapply(YesNo_varsH, as.character)
  XD_tmle = convert_factor_to_numeric(XD_missimput, c("is.DEN","Torniquete","Sexo",YesNo_varsH_char))
  
  #TMLE variable importance
  set.seed(102)
  #obtain marginal treatment effects of each dichotomized point treatment
  #tmleOut = tmle(Y=XD_tmle[,"DEN_dum"], A=XD_tmle[,"Sexo"], W=XD_tmle[,c("Cefalea","Dolor_ocular","Rash")], 
  #               Q.SL.library = c("SL.mean", "SL.knn", "SL.randomForest"), family="binomial")
  #summary(tmleOut)
  #tmleOut$estimates$RR$psi #estimate of relative risk (could alternatively request OR or ATE)
  #tmleOut$estimates$RR$pvalue #corresponding pvalue
  
  #Lasso variable importance
  if(method=="lasso"){
    set.seed(102)
    cvLasso = cv.glmnet(x=as.matrix(XD_tmle[,c(clinic_vars)]),
                        y=as.matrix(XD_tmle[,"DEN_dum"]),
                        family='binomial', alpha=1, nfolds=10, type.measure="auc")
    cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
    MFpreds = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs
    myVIM = as.data.frame(cbind(rownames(cvLcoefs), cvLcoefs))
    colnames(myVIM) = c("VarName",paste("Lasso",suffix,sep="_"))
  }
  
  return(myVIM)
}


#find most important variables
rfVIM_n88 = get_VIM(XD=comboD1_filter50n, method="rfVIM", suffix="n88", clinic_vars)
rfVIM_n1624 = get_VIM(XD=clinical_prelim, method="rfVIM", suffix="n1624", clinic_vars)
lassoVIM_n88 = get_VIM(XD=comboD1_filter50n, method="lasso", suffix="n88", clinic_vars)
lassoVIM_n1624 = get_VIM(XD=clinical_prelim, method="lasso", suffix="n1624", clinic_vars)

#combine variable importance measures into 1 chart
VIM_merged = merge(rfVIM_n88, rfVIM_n1624,  by="VarName", all=TRUE)
VIM_merged = merge(VIM_merged, lassoVIM_n88,  by="VarName", all=TRUE)
VIM_merged = merge(VIM_merged, lassoVIM_n1624,  by="VarName", all=TRUE)
VIM_merged
write.csv(x=VIM_merged, file=paste(resultsDir,"VIM_merged.txt"), row.names = FALSE)
