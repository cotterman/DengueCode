
#create a version of the data which replaces values for missing compounds with zeros
#other options: fill with min value observed for that compound,
#fill with 1/2 of the min value, input value
fill_LCMS_blanks = function(mydata){
  cindices = grep("MZ_",colnames(mydata))
  processedD = mydata
  #note: must have a return statement b/c assignment operator will not return anything (and default is 0)
  FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
  processedD[cindices] = lapply(mydata[cindices],FUN=FillMissingsWithZeros)
  return(processedD)
}

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


### function to read in the LC-MS text file and organize it into dataframe ###
clean_LCMS = function(infile, lcms_run, roundme=FALSE, decimals=2, printme=FALSE) {
  
  #infile="LCMS_saliva_Nica_50percent.txt"
  #infile = "LCMS_urine_Nica_50percent.txt"
  #lcms_run=5
  #printme=TRUE
  #decimals=2
  #roundme=FALSE
  
  abundancesD1 = read.delim(paste(lcmsCO_inputsDir, infile, sep=""), header=TRUE, nrows=200)
  nColumns = dim(abundancesD1)[2] #number of columns in abundance data
  #MZ = mass to charge of the compound, it usually is approximately +1 of the mass.  
  #This value distinguishes the type of metabolite.  
  #RT = retention times.
  #Resp = response (aka abundance).
  
  #obtained list of compounds (in same order in which they appear in spreadsheet columns)
  # take first row of values for variables MZ, MZ.1, MZ.2.,,,,MZ.2232
  MZ_Cols = seq(from=2, to=nColumns-2, by=3)
  #MZ_Cols #these are the column numbers that we will want to keep
  MZ_Nums = abundancesD1[1,MZ_Cols] #take the first row of data containing MZ values
  #MZ_Nums #these are the values in the specified columns
  
  #view MZ values in histogram
  MZ_Nums_reshaped1 = t(MZ_Nums) #min value is 107, max is 1427
  if(roundme==FALSE) {
    MZ_Names = paste("MZ_", sapply(MZ_Nums, as.character), sep="") #no rounding
  }
  takeDdecimals = function(x) round(x, digits=decimals)
  MZ_Nums_short = sapply(X=MZ_Nums, takeDdecimals)
  if(roundme==TRUE) {
    MZ_Names = paste("MZ_", sapply(MZ_Nums_short, as.character), sep="") #with rounding
  }
  
  #drop the RT and MZ columns and relabel the Resp columns with MZ names
  RespCols = seq(from=4, to=nColumns, by=3)
  respD = abundancesD1[,c(1,RespCols)]
  newnames = as.character(c("code", MZ_Names))
  #Note: when MZ values are repeated, R will give the repeated MZs a .1, .2, etc. suffix
     #ex: MZ_336.32513 and MZ_336.32513.1 and MZ_336.32513.2 will be in the data if MZ_336.32513 occurs 3 times
  colnames(respD) = newnames #rename columns
  
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  formatID = function(x, start, stop) sprintf("ID%04d",as.numeric(substr(as.character(x["code"]), start=start, stop=stop)))
  newIDs_1 = apply(X=respD, MARGIN = 1, FUN=formatID, start=13, stop=17)
  respD["code"] = newIDs_1
  
  #add an indicator of LCMS_run (to keep track once merged with other round of LC-MS)
  respD["LCMS_run"] = lcms_run
  
  #add an indicator for type of study (will need this later when we merge)
  if(lcms_run==1){
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run, clinical_inputsDir)
  } 
  
  #count  the number of duplicated MZ values
  if(printme==TRUE){
    print("Number of duplicated MZ values with no rounding of MZ values:")
    print(table(table(MZ_Nums_reshaped1))) #without rounding
    print(paste("Number of duplicated MZ values with rounding to", decimals, "decimals"))
    print(table(table(MZ_Nums_short))) #with rounding
  }
  
  #create a version of the data which replaces values for missing compounds with zeros
  respD = fill_LCMS_blanks(respD)
  
  return(respD)
  
}


### function takes abundance data as input and returns matrix of numeric MZ values (including duplicates if present)
get_MZ_vals = function(mydata){
  #obtain values of MZ numbers from processed abundance data
  temp = colnames(mydata) #all column names
  MZ_Names = temp[grep("MZ_",temp)] #take just the columns that start with "MZ_"
  #eliminate the extensions that indicate duplicate values --- these are always after the second decimal
  #(TODO: should verify there are no whole-number MZ values -- or reformat to eliminate whole numbers
  #will only be a problem if there are repeats of whole-number MZ values, as I wont detect them)
  #knock out the "MZ_" prefix and convert to number
  get_mz_nums = function(x){
    LocateDecimal2 = gregexpr(pattern ='\\.', x)[[1]][2]
    mystop = ifelse(is.na(LocateDecimal2), 50, LocateDecimal2-1)
    as.numeric(substr(x, start=4, stop=mystop))
  }
  MZ_Nums = lapply(X=MZ_Names, FUN=get_mz_nums) 
  MZ_Nums_reshaped = as.matrix(do.call(rbind.data.frame, MZ_Nums)) #convert list to matrix
  return(MZ_Nums_reshaped)
}

### function get list of patient IDs that are present in specified data
get_IDs_in_common = function(resp_D1D2, respD3_filter50n, respD5_filter50n){
  IDs_in_resp_D1_D2 = resp_D1D2[,c("code","Study","LCMS_run")]
  temp = merge(IDs_in_resp_D1_D2[,c("code","Study","LCMS_run")], 
               respD3_filter50n[,c("code","Study","LCMS_run")], by=c("code","Study"), all=T)
  temp$serum = temp$LCMS_run.x
  temp$LCMS_run.x = NULL
  temp$saliva = temp$LCMS_run.y
  temp$LCMS_run.y = NULL
  IDs_in_resp_all = merge(temp, respD5_filter50n[,c("code","Study","LCMS_run")], by=c("code","Study"), all=T)
  IDs_in_resp_all$urine = IDs_in_resp_all$LCMS_run
  IDs_in_resp_all$LCMS_run = NULL
  #are all saliva and urine samples are from same people? yes
  table(IDs_in_resp_all$urine, IDs_in_resp_all$saliva, useNA="ifany")
  #are all saliva and serum samples from the same people? just 1 (might be mistake -- ID1408 in cohort, which was already problematic)
  table(IDs_in_resp_all$serum, IDs_in_resp_all$saliva, useNA="ifany")
  return(IDs_in_resp_all)
}

### Function to create cateogrical variables using continuous counterparts ###
create_binary_variables = function(mydata){
  
  #define liver enlargement to be > 2cm
  mydata$is.LiverEnlarged = mydata$Higado>2
  mydata[which(is.na(mydata$Higado)), "is.LiverEnlarged"] = NA
  
  #Hypotension for age: 
    #Systolic blood pressure <80 mmHg for children <5 yrs old; <90 mmHg for children >=5 years old
  mydata$is.hypotension = (mydata$Presion_Arterial_Sist<80 & mydata$age<5) |
                          (mydata$Presion_Arterial_Sist<90 & mydata$age>=5)
  mydata[which( is.na(mydata$age) | is.na(mydata$Presion_Arterial_Sist) ), "is.hypotension"] = NA

  #Hypertention for age and gender 
    #approximate - source: http://emedicine.medscape.com/article/889877-overview, taken from NHLBI, 2004
    #todo: include height info and get more precise age cutoffs (I took midpoints)
  mydata$is.hypertension = (((mydata$age<3 &                 (mydata$Presion_Arterial_Sist>104 | mydata$Presion_Arterial_Dias>58)) | 
                          (mydata$age>=3 & mydata$age<9 &  (mydata$Presion_Arterial_Sist>111 | mydata$Presion_Arterial_Dias>74)) | 
                          (mydata$age>=9 & mydata$age<15 & (mydata$Presion_Arterial_Sist>123 | mydata$Presion_Arterial_Dias>80)) | 
                          (mydata$age>=15 &                (mydata$Presion_Arterial_Sist>129 | mydata$Presion_Arterial_Dias>84)) ) & mydata$Sexo=="female") |
                         (((mydata$age<3 &                 (mydata$Presion_Arterial_Sist>103 | mydata$Presion_Arterial_Dias>56)) | 
                         (mydata$age>=3 & mydata$age<9 &  (mydata$Presion_Arterial_Sist>114 | mydata$Presion_Arterial_Dias>74)) | 
                         (mydata$age>=9 & mydata$age<15 & (mydata$Presion_Arterial_Sist>123 | mydata$Presion_Arterial_Dias>81)) | 
                         (mydata$age>=15 &                (mydata$Presion_Arterial_Sist>136 | mydata$Presion_Arterial_Dias>87)) ) & mydata$Sexo=="male") 
  mydata[which(is.na(mydata$age) | (is.na(mydata$Presion_Arterial_Sist) & is.na(mydata$Presion_Arterial_Dias))), "is.hypertension"] = NA

  #define fever as temperature over 37.8
  mydata$is.fever = mydata$Temperatura>37.8
  mydata[which(is.na(mydata$Temperatura)), "is.fever"] = NA

  #freq_card can be dropped since we are including tachycardia and braquicardia, which is based on freq_card (listed in Lionels cutpoint list)
    #note: my age variable is continuous, so I have slightly different coding as compared to Lionel
  mydata$is.braquicardia = (mydata$age<1 & mydata$freq_card<109) | 
                           (mydata$age>=1 & mydata$age<4 & mydata$freq_card<90) | 
                           (mydata$age>=4 & mydata$age<6 & mydata$freq_card<65) | 
                           (mydata$age>=6 & mydata$age<9 & mydata$freq_card<60) | 
                           (mydata$age>=9 & mydata$age<12 & mydata$freq_card<60) | 
                           (mydata$age>=12 & mydata$age<17 & mydata$freq_card<60) #tod: add adult cutoff in case we get ppl above 16 yrs old
  mydata[which(is.na(mydata$age) | is.na(mydata$freq_card)), "is.braquicardia"] = NA
  
  #respiration cut point from Lionel's cutpoint list (different from "difficulty breathing"?)
  mydata$is.resp_normal = (mydata$age<1 & mydata$freq_resp<50) | 
                          (mydata$age>=1 & mydata$age<6 & mydata$freq_card<40) | 
                          (mydata$age>=6 & mydata$age<13 & mydata$freq_card<30) | 
                          (mydata$age>=13 & mydata$age<22 & mydata$freq_card<22) 
  mydata[which(is.na(mydata$age) | is.na(mydata$freq_card)), "is.resp_normal"] = NA

  #this is the one categorical "general hemorrhaging" variable
  mydata$is.torniquete10plus = (as.numeric(mydata$Torniquete)>=2)
  mydata$is.torniquete20plus = (as.numeric(mydata$Torniquete)>=3)
  mydata[which(is.na(mydata$Torniquete)), "is.torniquete20plus"] = NA #error
  mydata[which(is.na(mydata$Torniquete)), "is.torniquete10plus"] = NA #error

  #this is the one categorical "general symptom" variable
  mydata$is.pulse_rapid = (mydata$Pulso=="rapid")
  mydata$is.pulse_strong = (mydata$Pulso=="strong")
  mydata[which(is.na(mydata$Pulso)), "is.pulse_rapid"] = NA
  mydata[which(is.na(mydata$Pulso)), "is.pulse_strong"] = NA

  #age breaks inspired by CART (1.2, 2.5, 3.3, 4.6, 5.5, 6.5, 7.1, 9.5, 11, 12)
    #and by ctree (3.3, 4.6)
  #age is never missing
  mydata$is.age_under1 = (mydata$age<1) #infant
  mydata$is.age_1or2 = (mydata$age>=1 & mydata$age<3) #toddler
  mydata$is.age_3or4 = (mydata$age>=3 & mydata$age<5) #young child
  mydata$is.age_5to8 = (mydata$age>=5 & mydata$age<9) #old child
  mydata$is.age_9to12 = (mydata$age>=9 & mydata$age<13) #preteen
  mydata$is.age_13to15 = (mydata$age>=13 & mydata$age<16) #teenager
  
  #DaysSick breaks inspired by ctree (3)
  #DaysSick is never missing
  mydata$is.DaysSick_under4 = (mydata$DaysSick<4) #febrile
  mydata$is.DaysSick_4to6 = (mydata$DaysSick>=4 & mydata$DaysSick<=6) #critical
  mydata$is.DaysSick_7plus = (mydata$DaysSick>=7) #recovery

  #combine hemarrhage variables into 1 indicator(?)
  
  #is rash the same as hematoma?
  
  ## add variables that other methods say are important
  #meet definition of thrombocytopenia if <100,000 platelets per mm3.  (Remember Plaquetas variable is /1000)
  mydata$is.thrombocytopenia = mydata$Plaquetas<100
  mydata[which(is.na(mydata$Plaquetas)), "is.thrombocytopenia"] = NA
  #define leukopenia to be when white blood count <= 5000 cells/mm3
  mydata$is.leukopenia = mydata$Leucocitos<=5
  mydata[which(is.na(mydata$Leucocitos)), "is.leukopenia"] = NA
  
  #age categories inspired by Hope's paper
  mydata$age_cat = cut(mydata$age, c(0,1,4,9,16), right=FALSE)
  #day of illness categories
  mydata$DaysSick_cat = cut(mydata$DaysSick, c(6,15), right=FALSE)

  #define gallbladder enlargement to be > 2mm (reasonable?)
  #mydata$is.GallbladderEnlarged = mydata$Engrosamiento_mm>2

  return(mydata)
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
  XD = clinical_full_clean
  vars_to_impute = varlist
  predictor_vars_prelim = varlist
  method="RF"
  
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
  

######################################################

### knn (run thru SuperLearner) does not like factor variables with 2 levels, so convert to numeric.
#convert all yes/no questions into factors (regression algorithms may run suboptimal if coded as numeric)
convert_factor_to_numeric = function(data, vars_to_convert){
  #print(vars_to_convert)
  for (name in vars_to_convert) {
    data[,name] = as.numeric(as.character(data[,name]))
  } 
  return(data)
}