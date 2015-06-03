
### function to reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode = function(x, var) sprintf("ID%04d",as.numeric(as.character(x[var])))

### function to give us the duplicated values in input vector
get_dupIDs = function(vec, vec_name){
  #verify there are no duplicate code values in either dataset
  t1 = table(vec)
  v1dups = table(t1) 
  if(dim(v1dups)>1){
    print(paste("Counts of duplicates in", vec_name))
    print(v1dups)
    print(paste("Duplicated values in", vec_name))
    print(t1[which(t1>1)])
  }else{
    print(paste("No duplicated values in", vec_name))
  }
}

### function to give us the values which are not found in both vectors
    # also checks for duplicate values in each dataset
get_mismatches = function(v1, v2, v1_name, v2_name){
  #v1 and v2 should be vectors --- ex: respD[,"code"]
  #v1_name should be a string describing v1 --- ex: "LCMS"
  
  get_dupIDs(v1, v1_name)
  get_dupIDs(v2, v2_name)
  
  #check out values that are not in both datasets
  v1_not_v2 = v1[v1 %in% v2==FALSE] #v1 values not found in v2
  v2_not_v1 = v2[v2 %in% v1==FALSE] #v2 values not found in v1
  if(length(v1_not_v2)>0){
    print(paste(v1_name,"values not found in",v2_name))
    print(v1_not_v2)
  }else{
    print(paste("All values in",v1_name,"are also in",v2_name))
  }
  if(length(v2_not_v1)>0){
    print(paste(v2_name,"values not found in",v1_name))
    print(v2_not_v1)
  }else{
    print(paste("All values in",v2_name,"are also in",v1_name))
  }
}

### function to read in the file indicating which study each sample belongs to ###
get_study = function(respD_only, lcms_run, clinical_inputsDir){
  
  # process file that contains indicator of study (cohort versus hospital)
  if(lcms_run==2){
    study_info = read.delim(paste(clinical_inputsDir,"Study classifications_batches 3 and 4.txt", sep=""), header=TRUE, nrows=200)
  } else if(lcms_run==3){
    study_info = read.delim(paste(clinical_inputsDir,"Study classifications_saliva.txt", sep=""), header=TRUE, nrows=200)
  } else if(lcms_run==5){
    study_info = read.delim(paste(clinical_inputsDir,"Study classifications_urine.txt", sep=""), header=TRUE, nrows=200)
  }
  newnames = as.character(c("Study","code","Res_Final2","OMS"))
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  fcode =    function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
  code_list = apply(X=study_info, MARGIN = 1, FUN=fcode)
  study_info["code"] = code_list
  
  #merge with lc-ms data
  mismatches = get_mismatches(v1=study_info[,"code"], v2=respD_only[,"code"], v1_name="study info", v2_name="LCMS")
  respD_new = merge(study_info[,c("Study","code")], respD_only, by="code") #this returns the intersection of each dataset
  get_dupIDs(respD_new[,"code"], "LCMS") #may still have duplicates
  
  return(respD_new)
}

#function takes variables in "predictors" to imput "var" in dataframe "XD"
fit_and_imput = function(var, method, XD, predictors){
  cat(paste("imputation for", var,"\n"))
  myformula = as.formula(paste(var, " ~ ", paste(predictors, collapse= "+")))
  #fit model (use obs with non-missing Y value, of course)
  #SL doesn't deal well with factors - it will classify binary outcomes but otherwise will only do regression. 
      #todo: learn how to change parameter settings of functions called by SL (should be able to do categorical classifications) 
      #for now, just use RF for the 2 categorical variables in the data (actually, only Pulso has missing values to impute)
  if(method=="RF" | var=="Pulso" | var=="Torniquete"){
    #With RF, if Y is a factor, classification is assumed, otherwise regression is assumed
       #na.action=na.omit means to omit obs that have a missing value for one or more predictors
       #na.action=na.fail means RF will not run if any of the predictors have a missing value
    rf = randomForest(myformula, data=XD[which(!is.na(XD[var])),], ntree=501, na.action=na.fail) 
    #fill var with predicted value when var is missing
      #type=response says to choose class with majority vote (as opposed to fitted probabilities) if RF is a classification object
        #will give the usual regression result if RF is a classification object 
    XD[which(is.na(XD[var])),var] = predict(rf, XD[which(is.na(XD[var])),predictors], type="response") 
  }
  if(method=="SL"){
    if(is.factor(XD[,var])){
      
    }
    
  }
  return(XD[var])
}

#function takes vector "mycolumn" and returns a vector of 0 and 1s to indicate where mycolumn has missing values
create_missing_indicator = function(mycolumn){
  #dummy will equal 1 when var is missing (values will be imputed for these observations) and 0 otherwise
  dummy=as.vector(rep(x=0, length=length(mycolumn)))
  dummy[which(is.na(mycolumn))] = 1
  return(dummy)
}

#function creates dummy variables to indicate missingness and then replaces missings with imputed values
impute_missings = function(XD, vars_to_impute, predictor_vars_prelim, exclude_miss_vars, method){
  #XD is dataframe. This is the dataframe that will be altered.
  #vars_to_impute is list of variable names. If any of them contain missing values, then values will be imputed.
  #predictor_vars_prelim is preliminary list of names of variables to consider as predictors for missing values.
  #exclude_miss_vars:
  # = True (implemented) if you wish to use the covaiate vars which have no missing vals in XD to imput missing values.
  #So variable list for prediction will be fixed regardless of which variable we are trying to imput
  #also, all observations will be used (as opposed to dropping those with missing covariates)
  # = False (not yet implemented) if you wish to use all covaiate vars which have no missing vals for given obs to imput missing values.
  #So variable list for prediction will vary from observation to observation, depending on what is non-missing.
  #If prediction method cannot use obs with missing values, then prediction model will not use obs with missing values
  #method:
  # = "RF" (implemented) if want to use random forest to impute values
  # = "SL" (not yet implemented - todo -- see scratch code for start) if want to use super learner to impute values
  set.seed(200)
  
  #testing
  #XD = clinical_full_clean
  #vars_to_impute = varlist
  #predictor_vars_prelim = varlist
  #method="RF"
  
  #get list of vars to impute (of vars_to_impute, find those which contain missings)
  to_impute_missCount = sapply(XD[,c(as.character(vars_to_impute))], function(x) sum(is.na(x)))
  vars_to_impute = names(which(to_impute_missCount>0))
  
  #get list of vars to use as predictors (of predictor_vars_prelim, find nomiss_vars (to use as predictors), starting from vars_to_use list
  if(exclude_miss_vars==T){
    missCount = sapply(XD[,c(as.character(predictor_vars_prelim))], function(x) sum(is.na(x)))
    predictors = names(which(missCount==0))
  }
  
  #create dummy variables to indicate when vars_to_impute are missing
  dum_var_names = lapply(vars_to_impute, function(x) paste(x,"imputed",sep="_"))
  XD[as.character(dum_var_names)]  = lapply(X = XD[vars_to_impute], FUN=create_missing_indicator)  
  
  #fit model and predict values of var when missing
    #todo: make exception for IR and PCR such that they are only imputed using DENV patients
  XD[vars_to_impute] = sapply(X=vars_to_impute, FUN=fit_and_imput, method, XD, predictors)
  return(XD)
}
  
### knn (run thru SuperLearner) does not like factor variables with 2 levels, so convert to numeric.
#convert all yes/no questions into factors (regression algorithms may run suboptimal if coded as numeric)
convert_factor_to_numeric = function(data, vars_to_convert){
  #print(vars_to_convert)
  for (name in vars_to_convert) {
    data[,name] = as.numeric(as.character(data[,name]))
  } 
  return(data)
}