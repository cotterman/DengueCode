
###############################################################################
####################### Data from Mexico ######################################
###############################################################################  

### function explores the limited info we have on patients from Mexico
clean_clinical_mexico = function(){
  # This is the meta data that I extracted from the .d file names 
  clinical_mexicoD = read.csv(paste(clinical_inputsDir,"meta data for mexico serum samples.csv", sep=""), header=TRUE, nrows=2000)
  #drop the strange samples
  char_counter = function(x) {nchar(as.character(x))}
  char_counts = sapply(X=clinical_mexicoD$code, FUN=char_counter)
  temp = clinical_mexicoD[which(char_counts==5 & clinical_mexicoD$code!="Beaty"),]
  #obtain number of unique samples
  print(table(table(temp$code)))
  #remove repeated patient codes 
  temp = temp[!duplicated(temp$code),]
  #summarize de-dupped data
  summary(temp)
  table(temp$DaysSick)
  #estimate number of patients who end with DHF but did not have DHF
  #when they appeared to clinic (so days 0 - 2 at most) -- 10 patients
  table(temp[which(temp$DxFinal4cat=="DHF"),"DaysSick"])
  return(temp)
}


####################################################################
################ Clean clinical data from Nicaragua ################
####################################################################

### function that returns either min or max (=specified by myfunc) of the myvar.x and myvar.y values
return_min_or_max = function(mydata, myvar, myfunc){
  #myfunc can be min or max
  myvar.y = c(paste(myvar,".y",sep=""))
  myvar.x = c(paste(myvar,".x",sep=""))
  mydata[,myvar] = apply(mydata[, c(myvar.x, myvar.y)], 1, function(x) myfunc(x, na.rm=T))
  #replace -Inf and Inf values with NAs
  mydata[,myvar][mydata[,myvar]==Inf | mydata[,myvar]==-Inf] = NA
  #exclude redundant vars
  col_myvars.x = names(mydata) %in% c(myvar.x)
  mydata = mydata[!col_myvars.x] 
  col_myvars.y = names(mydata) %in% c(myvar.y)
  mydata = mydata[!col_myvars.y] 
  return(mydata)
}

### function that takes the min or max of myvar.x and myvar.y if they exist
test_then_return_min_or_max = function(mydata, myvar, myfunc){
  #take the min/max of the .y and .x columns if they exist
  myvar.y = c(paste(myvar,".y",sep=""))
  myvar.x = c(paste(myvar,".x",sep=""))
  if(length(grep(myvar.x, names(mydata))) > 0 & 
       length(grep(myvar.y, names(mydata))) > 0 ){
    return(return_min_or_max(mydata, myvar, myfunc))
  }
  #if variable does not have .x and .y versions then just return data unchanged
  else{
    return(mydata)
  }
}

### Function to obtain list of clinical variables based on user-supplied criteria
get_clinic_var_list = function(clinic_varsD, outcome, eliminate_vars_with_missings=F, 
                               eliminate_constant_vars=F, eliminate_vars_with_minXnomiss=0, restrict_to_binary=F,
                               XD, restrict_to_cohort_vars=F, restrict_to_hospit_vars=F, UltraX=T, BloodLab=T){
  #clinic_varsD: spreadsheet containing universe of clinical/lab vars and info about them
  #outcome: "ND.vs.DEN" or "DF.vs.DHF.DSS" or "either"
  #eliminate_vars_with_missings: T if you wish to restrict to variables that are never missing (else F)
  #eliminate_vars_with_minXnomiss: set equal to the minimum number of non-missing values that must be present 
  #restrict_to_binary: set to T if you want only binary variables included in output
  #eliminate_constant_vars: T if you wish to eliminate variables take the same value for observations
  #(prediction methods will give error if trying to predict variable that is always 1 particular value)
  #XD: dataset containing patient clinical info (relevant for determining missingness)
  #restrict_to_cohort_vars: T if you wish to only include variables that are (generally) present in cohort data
  #restrict_to_hospit_vars: T if you wish to only include variables that are (generally) present in hospital data
  #UltraX: T if you wish to include the ultrasound and x-ray variables.  F otherwise.
  #BloodLab: T if you wish to include the blood lab variables.  F otherwise.
  
  if(outcome=="ND.vs.DEN"){
    #obtain list of candidate clinical variables to include in ND vs. DEN prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                              "variable.name.in.final.data"]))
    clinicCount = length(clinic_vars)
    cat("\n Number of clinical variables eligible for ", outcome, " analysis before drops:")
    print(clinicCount) #view 
  }
  if(outcome=="DF.vs.DHF.DSS"){
    #obtain list of candidate clinical variables to include in DF vs. DHF/DSS prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),
                                              "variable.name.in.final.data"])) #30
    #drop ND observations (will affect which vars we consider missing)
    XD = XD[!is.na(XD$DHF_DSS_dum),]
    clinicCount = length(clinic_vars)
    cat("\n Number of clinical variables eligible for ", outcome, " analysis before drops:")
    print(clinicCount) #view 
  }
  if(outcome=="either"){
    #obtain list of candidate clinical variables to include in either DF vs. DHF/DSS or ND vs. DEN prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1
                                                    | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),
                                              "variable.name.in.final.data"])) 
    clinicCount=length(clinic_vars)
    cat("\n Number of clinical variables eligible for ", outcome, " analysis before drops:")
    print(clinicCount) #view 
  }
  if(eliminate_constant_vars==T){
    variances = sapply(XD[,clinic_vars], var_after_conversion)
    cvars_constants = names(variances[which(variances==0 | is.na(variances)==T)])
    cat("\n Clinical variables that take same value for all observations \n")
    print(cvars_constants) #view 
    cvars_noConstants = names(variances[which(variances!=0)])
    cat("\n Clinical variables that take have some variation \n")
    print(cvars_noConstants) #view 
    clinic_vars = cvars_noConstants
  }
  if(eliminate_vars_with_missings==T){
    #obtain list of clinical variables that have no missing values in data
    clinic_vars_missCount = sapply(XD[clinic_vars], function(x) sum(is.na(x)))
    cvars_noMiss = clinic_vars[clinic_vars_missCount==0] #40 
    clinicCount = length(cvars_noMiss)
    cat("\n Clinical variables with no missings \n")
    print(cvars_noMiss) #view 
    cat("\n Clinical variables with missings \n")
    print(clinic_vars[clinic_vars_missCount>0]) #view 
    clinic_vars = cvars_noMiss
  }
  if(eliminate_vars_with_minXnomiss>0){
    #obtain list of clinical variables that have fewer than 50 non-missing values
    clinic_vars_nomissCount = sapply(XD[clinic_vars], function(x) sum(!is.na(x)))
    cvars_Xp = clinic_vars[clinic_vars_nomissCount>=eliminate_vars_with_minXnomiss] 
    clinicCount = length(cvars_Xp)
    #cat("\n Clinical variables with at least X non-missing values \n")
    #print(cvars_Xp) #view 
    cat("\n Clinical variables with fewer than ",eliminate_vars_with_minXnomiss, " non-missing values \n")
    print(clinic_vars[clinic_vars_nomissCount<eliminate_vars_with_minXnomiss]) #view 
    clinic_vars = cvars_Xp
  }
  if(restrict_to_cohort_vars==T){
    cohort_vars = c(as.character(clinic_varsD[which(clinic_varsD$in.cohort.data==1 & clinic_varsD$variable.name.in.final.data!=0),
                                              "variable.name.in.final.data"]),"age","DaysSick") 
    #take the variables that have satisfied above criteria and are in cohort data
    clinic_vars_inCohort = clinic_vars[clinic_vars %in% cohort_vars==TRUE]
    v1_not_v2 = clinic_vars[clinic_vars %in% clinic_vars_inCohort==FALSE] #clinic_vars values not found in clinic_vars_inCohort
    cat("\n Clinical variables dropped due to not being collected for cohort samples \n")
    print(v1_not_v2)
    clinic_vars = clinic_vars_inCohort
  }
  if(restrict_to_hospit_vars==T){
    hospit_vars = c(as.character(clinic_varsD[which(clinic_varsD$in.hospit.data==1 & clinic_varsD$variable.name.in.final.data!=0),
                                              "variable.name.in.final.data"]),"age","DaysSick") 
    #take the variables that have satisfied above criteria and are in hosptial data
    clinic_vars_inHospit = clinic_vars[clinic_vars %in% hospit_vars==TRUE]
    v1_not_v2 = clinic_vars[clinic_vars %in% clinic_vars_inHospit==FALSE] #clinic_vars values not found in clinic_vars_inCohort
    cat("\n Clinical variables dropped due to not being collected for hospital samples \n")
    print(v1_not_v2)
    clinic_vars = clinic_vars_inHospit
  }
  if(UltraX==F){
    UltraXvars = c(as.character(clinic_varsD[
      which(clinic_varsD$Variable.Category=="Laboratory-Ultrasound" | 
              clinic_varsD$Variable.Category=="Laboratory-Ultrasound and X-ray" | 
              clinic_varsD$Variable.Category=="Laboratory-X ray"),c("variable.name.in.final.data")]))
    clinic_vars_notUltraX = clinic_vars[clinic_vars %in% UltraXvars==FALSE]
    v1_and_v2 = clinic_vars[clinic_vars %in% UltraXvars==TRUE] #clinic_vars values found in UltraXvars
    cat("\n Clinical variables dropped due to being an ultrasound of Xray variable \n")
    print(v1_and_v2)
    clinic_vars = clinic_vars_notUltraX
  }
  if(BloodLab==F){
    BloodLabvars = c(as.character(clinic_varsD[
      which(clinic_varsD$Variable.Category=="Laboratory-Clinical lab" | 
              clinic_varsD$Variable.Category=="Laboratory-Urine analysis" | 
              clinic_varsD$Variable.Category=="Laboratory-Virological"),c("variable.name.in.final.data")]))
    clinic_vars_notBloodLab = clinic_vars[clinic_vars %in% BloodLabvars==FALSE]
    v1_and_v2 = clinic_vars[clinic_vars %in% BloodLabvars==TRUE] #clinic_vars values found in BloodLabvars
    cat("\n Clinical variables dropped due to being an blood lab variable \n")
    print(v1_and_v2)
    clinic_vars = clinic_vars_notBloodLab
  }
  if(restrict_to_binary==T){
    Binaryvars = c(as.character(clinic_varsD[
      which(clinic_varsD$CC_type=="binary"),c("variable.name.in.final.data")]))
    clinic_vars_notBinary = clinic_vars[clinic_vars %in% Binaryvars==FALSE]
    v1_and_v2 = clinic_vars[clinic_vars %in% Binaryvars==TRUE] #clinic_vars values found in BloodLabvars
    cat("\n Clinical variables dropped due to not being a binary variable \n")
    print(clinic_vars_notBinary)
    clinic_vars = v1_and_v2
  }
  cat("\n Final list of ", length(clinic_vars), " clinical variables:\n")
  print(clinic_vars)
  return(clinic_vars)
}

### Function to create cateogrical variables using continuous counterparts -- applicable to cohort and hospital datas ###
create_binary_variables = function(mydata){
  
  #define liver enlargement to be > 2cm
  mydata$is.LiverEnlarged = mydata$Higado>2
  #mydata[which(is.na(mydata$Higado)), "is.LiverEnlarged"] = NA #not necessary -- NAs already created as desired
  
  
  #Hypotension for age -- used for DSS classification
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
  
  #define fever as temperature over 37.5
  mydata$is.fever = mydata$Temperatura>=37.5
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
  mydata$is.torniquete20plus = (as.numeric(mydata$Torniquete)>=3) #used in DHF classification
  mydata[which(is.na(mydata$Torniquete)), "is.torniquete20plus"] = NA #error
  mydata[which(is.na(mydata$Torniquete)), "is.torniquete10plus"] = NA #error
  
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
  
  #is rash the same as hematoma?
  
  ## add variables that other methods say are important
  #meet definition of thrombocytopenia if <100,000 platelets per mm3.  (Remember Plaquetas variable is /1000)
  mydata$is.thrombocytopenia = mydata$Plaquetas<=100 #used in DHF classification
  mydata[which(is.na(mydata$Plaquetas)), "is.thrombocytopenia"] = NA
  #define leukopenia to be when white blood count <= 5000 cells/mm3
  mydata$is.leukopenia = mydata$Leucocitos<=5
  mydata[which(is.na(mydata$Leucocitos)), "is.leukopenia"] = NA
  
  #vars used for DSS definition
  mydata$is.Estrechamiento = with(mydata,ifelse(Presion_Arterial_Sist-Presion_Arterial_Dias<=20,1,0))

  #age categories inspired by Hope's paper
  mydata$age_cat = cut(mydata$age, c(0,1,4,9,16), right=FALSE)
  #day of illness categories
  mydata$DaysSick_cat = cut(mydata$DaysSick, c(6,15), right=FALSE)
  
  mydata$is.FirstInfection = mydata$IR=="primary"
  
  mydata$is.serotype1 = mydata$PCR==1
  mydata$is.serotype2 = mydata$PCR==2
  mydata$is.serotype3 = mydata$PCR==3
  
  #define gallbladder enlargement to be > 2mm (reasonable?)
  #mydata$is.GallbladderEnlarged = mydata$Engrosamiento_mm>2
  
  return(mydata)
}


### function produces clean clinical data based on the "24 hr" clinical data
clean_clin24_data = function(clinic_varsD, IDs_in_resp_all){
  
  
  ####### List of clinical/lab variables to consider for analysis ########
  
  #names of variables from hospital data to include in either analysis (these are what will mostly become the standardized names)
  clinic_varsH = clinic_varsD[which(clinic_varsD$in.hospit.data==1 & 
                                   (clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1))
                                   ,"Variable.Name.in.Hospital.data"]
  clinic_varsH_char = sapply(clinic_varsH, as.character) #convert factor to character
  
  #names of variables in cohort data (these names will be changed to match the hospital data variable names)
  clinic_varsC = clinic_varsD[which(clinic_varsD$in.cohort.data==1 & 
                                   (clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1))
                                   ,"Variable.Name.in.Cohort.data"]
  clinic_varsC_char = sapply(clinic_varsC, as.character)


  
  ###############################################################################
  ##################### Data from the hospital study ############################
  ###############################################################################
  
  # This should contain "first 24 hr" data from the first batch and from some of the 3rd/4th batches
  clinical_hospitD = read.delim(paste(clinical_inputsDir,"Clinical data_hospital study.txt", sep=""), header=TRUE, nrows=2000)
  clinical_hospitD["code"] = apply(X=clinical_hospitD, MARGIN = 1, FUN=fcode, var="code")
  
  #code the yes/no variables to take numeric values 1/0 (to be consistent with cohort data)
  #in future, we convert all yes/no variables to logical vectors (booleans)
  YesNo_varsH = clinic_varsD[which((clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1) & 
                                     clinic_varsD$Variable.Output=="Yes, No" & clinic_varsD$in.hospit.data==1),
                             "Variable.Name.in.Hospital.data"]
  YesNo_varsH_char = sapply(YesNo_varsH, as.character)
  for (name in YesNo_varsH_char) {
    clinical_hospitD[which(clinical_hospitD[name]=="Si"),"temp"]=1
    clinical_hospitD[which(clinical_hospitD[name]=="No"),"temp"]=0
    class(clinical_hospitD$temp)
    #drop old version and rename new version
    clinical_hospitD[name] = NULL
    clinical_hospitD[name] = clinical_hospitD$temp
    clinical_hospitD$temp = NULL
  }
  
  other_varsH = clinic_varsD[which((clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1)
                                   & clinic_varsD$Variable.Output!="Yes, No"),
                             c("Variable.Name.in.Hospital.data","Variable.Output")]
  other_varsH #see what is left
  other_varsH_char = sapply(other_varsH["Variable.Name.in.Hospital.data"], as.character)
  #make sure these variables are compatible between datasets
  #for (name in other_varsH_char) {
  #  print(name)
  #  print(summary(clinical_hospitD[name]))
  #  print(summary(clinical_cohortD_clean[name]))
  #}
  
  #make the PCR variable a factor variable (else it will mismatch with cohort data and generate error when merged)
  #also, give ND patients a PCR value of 0.  (NAs will be given to indeterminants)
  clinical_hospitD[which(clinical_hospitD$DENV=="Negativo"),"PCR"]=0
  clinical_hospitD$PCR = factor(x = as.character(clinical_hospitD$PCR))
  table(clinical_hospitD$PCR)
  
  #make IR variable consistent with cohort data
  clinical_hospitD[which(clinical_hospitD$IR=="Secundario"),"new_IR"]="S"
  clinical_hospitD[which(clinical_hospitD$IR=="Primario"),"new_IR"]="P"
  clinical_hospitD[which(clinical_hospitD$IR=="Indeterminado"),"new_IR"]="I"
  clinical_hospitD["IR"] = NULL
  clinical_hospitD$IR = factor(x = as.character(clinical_hospitD$new_IR))
  table(clinical_hospitD$IR)
  class(clinical_hospitD$IR)
  
  #turn Torniquete into factor (rather than numeric) variable
  clinical_hospitD$temp = factor(x = as.character(clinical_hospitD$Torniquete))
  clinical_hospitD["Torniquete"] = NULL
  clinical_hospitD["Torniquete"] = clinical_hospitD$temp
  clinical_hospitD$temp = NULL
  table(clinical_hospitD$Torniquete)
  class(clinical_hospitD$Torniquete)
  
  #add an indicator for type of study (will need this later when we merge)
  clinical_hospitD["Study"] = "Hospital"
  
  #calculate age of patient as of sampling date (first reformat variables into dates)
  clinical_hospitD$FTM = as.Date(clinical_hospitD$FTM, format = "%m/%d/%Y")
  clinical_hospitD$Fecha_Nac = as.Date(clinical_hospitD$Fecha_Nac, format = "%m/%d/%Y")
  clinical_hospitD$FIS = as.Date(clinical_hospitD$FIS, format = "%m/%d/%Y")
  clinical_hospitD$FIF = as.Date(clinical_hospitD$FIF, format = "%m/%d/%Y")
  clinical_hospitD$age = as.numeric(clinical_hospitD$FTM - clinical_hospitD$Fecha_Nac)/365.25
  #calculate days from fever onset to sampling date (changed from days since symptom onset on 11-25-2014)
  clinical_hospitD$DaysSick = as.numeric(clinical_hospitD$FTM - clinical_hospitD$FIF) + 1
  #use FIS if FIF is missing
  clinical_hospitD[which(is.na(clinical_hospitD$FIF)),]$DaysSick = 
    as.numeric(clinical_hospitD[which(is.na(clinical_hospitD$FIF)),]$FTM - clinical_hospitD[which(is.na(clinical_hospitD$FIF)),]$FIS) + 1
  #check calculation
  #View(clinical_hospitD[,c("FTM","Fecha_Nac","FIS","age","DaysSick")])
  
  #keep only relevant variables
  nonclinic_keepervars = c("code","Study","FTM","age","DaysSick", 
                           "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final","Fecha_Nac")
  namesH_nonLCMS = as.character(c(clinic_varsH_char, nonclinic_keepervars))
  #drop the one guy who is missing "DaysSick" since this variable is too important
  clinical_hospitD_prelim = clinical_hospitD[,namesH_nonLCMS]
  #rename some vars
  names(clinical_hospitD_prelim)[names(clinical_hospitD_prelim)=="ClasificacionPrimerDia"] <- "WHO24h_given"
  names(clinical_hospitD_prelim)[names(clinical_hospitD_prelim)=="Res_Final"] <- "DENV"
  
  ### Re-create initial diagnosis variable ###
  DxDHF_varsH = list("Torniquete","Proteina_tot","Albumina","Plaquetas",
                     "PetequiasEspontaneas","Equimosis","Purpura","Hematoma",
                     "Venopuncion","Hipermenorrea","Vaginal","Subconjuntival","Hemoptisis","Nariz","Encias",
                     "Melena","Hematemesis",
                     "Derrame_","Ascitis","Edema_","Hemoconcentracion","HematocritoElev")
  get_mismatches(colnames(clinical_hospitD_prelim),DxDHF_varsH,"in_data","needed") #we have all needed vars
  DxDSS_varsH = list("Presion_Arterial_Sist","Presion_Arterial_Dias", 
                     "ExtremidadesFrias","Palidez","Sudoracion","Escalofrio","Llenado_Capilar")
  get_mismatches(colnames(clinical_hospitD_prelim),DxDSS_varsH,"in_data","needed") #we have all needed vars
 
  #DHF
  #creates is.torniquete20plus, is.thrombocytopenia, is.Estrechamiento, is.hypotension (among others)
  temp = create_binary_variables(clinical_hospitD_prelim) 
  
  #this is the one categorical "general symptom" variable -- not in cohort data
  temp$is.pulse_rapid = (temp$Pulso=="rapid")
  temp$is.pulse_strong = (temp$Pulso=="strong")
  temp[which(is.na(temp$Pulso)), "is.pulse_rapid"] = NA
  temp[which(is.na(temp$Pulso)), "is.pulse_strong"] = NA
  temp$is.pulse_danger = (temp$Pulso=="R" | temp$Pulso=="N") #danger if rapid or not palpable (#used in DSS classification)
  
  #these are defined according to appearance in definition of DHF -- not in cohort data
  temp$is.Hypoproteinemia = temp$Proteina_tot<4.2  # <4.2 g/dl 
  temp[which(temp$age>=2),"is.Hypoproteinemia"] = temp[which(temp$age>=2),"Proteina_tot"]<6 #looser criteria if age>=2 
  temp$is.Hypoalbuminemia = temp$Albumina<2  # <2 g/dl
  temp[which(temp$age>=1),"is.Hypoalbuminemia"] = temp[which(temp$age>=1),"Albumina"]<3.5  #looser criteria if age>=1
  #DSS definition -- not in cohort data
  temp$is.coldclammy = with(temp, ifelse(ExtremidadesFrias + Palidez + Sudoracion + Escalofrio >0 , 1, 0))
  
  #note: sum will treat NAs as zeros.  So if all vars are NAs (which never happens), then the sum will be zero rather than NA
  temp$Hem1 = apply(temp[,c("PetequiasEspontaneas","Equimosis","Purpura","Hematoma")], 1, function(x) sum(x, na.rm=T)>0)
  temp$Hem2 = apply(temp[,c("Venopuncion","Hipermenorrea","Vaginal","Subconjuntival","Hemoptisis","Nariz","Encias")], 1, function(x) sum(x, na.rm=T)>0) # Menorrhagia, vaginal bleeding, subconjunctival bleeding, hemoptysis, epistaxis, or gingivorrhagia
  temp$Hem3 = apply(temp[,c("Melena","Hematemesis")], 1, function(x) sum(x, na.rm=T)>0) #Melena or hematemesis
  #Pleural Effusion is Derrame_ #check this?? PleuralEffusionLeftRx,  PleuralEffusionRightRx, PulmonaryEdemaRx, PleuralEffusionClinical, PleuralEffusionRightUS or PleuralEffusionLeftUS
  #make sure the "Edema_" variable includes "Edema, periorbital, facial, limbs, hydrocele, generalized and fovea"
  #If (Fever = Yes) AND  ([Tourniquet test] = Yes or [HEM1] = Yes or [HEM2] = Yes or [HEM3] = Yes) AND (Thrombocyto-penia = Yes) 
  #  AND ([Hemoconcentration] = Yes or [ElevHematocrit] = Yes  or [PlerualEffusion] = Yes or [Ascites] = Yes or [Hypoprotein-emia] = Yes or [Hypoalbumin-emia] =Yes or [Edema] = Yes), then [DHF] = Yes
  temp$Leakage = apply(temp[,c("Derrame_","Ascitis","Edema_","is.Hypoproteinemia","is.Hypoalbuminemia","Hemoconcentracion","HematocritoElev")], 1, function(x) sum(x, na.rm=T)>0)
  #I took out is.fever since it can be biphasic (patient could have prior history of fever)
  #when is.thrombocytopenia is NA, can give a 0 but will never give a 1 (gives NA instead if other criteria are satisfied and is.thrombocytopenia is NA -- reasonable)
  temp$DHF = with(temp, ifelse((is.torniquete20plus | Hem1 | Hem2 | Hem3) & is.thrombocytopenia & Leakage, 1, 0 ))
  #View(temp[,c("is.torniquete20plus","Hem1","Hem2","Hem3","is.thrombocytopenia","Leakage","DHF")]) 
  tabletot(temp, "DHF", "WHO24h_given", useNA="always") #poor match -- I have far fewer DHF cases
  
  #DSS (is.Estrechamiento, is.hypotension, is.pulse_danger, is.coldclammy created in "create_binary_variables")
  temp$DSS = with(temp, ifelse( DHF==1 & (is.Estrechamiento | is.hypotension) & (is.pulse_danger | is.coldclammy | Llenado_Capilar), 1, 0 ))
  
  #assume for now all patients in data satisfy the DF (suspected dengue) criteria
  temp$WHO24hr4cat = ifelse(temp$DSS==1,"DSS",ifelse(temp$DHF==1,"DHF","DF")) 
  table(temp$WHO24hr4cat, temp$WHO24h_given, useNA="always") #poor match
  
  ## Additional outcome indicators (see if my DHF indicator matches with Lionel's 12 hr Dx indicator) ##
  outcomes_hospitD = read.delim(paste(clinical_inputsDir,"Additional outcomes_hospital study.txt", sep=""), header=TRUE, nrows=2000)
  outcomes_hospitD$CareLevel12h = substr(outcomes_hospitD$CareLevel12h, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
  outcomes_hospitD$CareLevelFinal = substr(outcomes_hospitD$CareLevelFinal, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  outcomes_hospitD$code = apply(X=outcomes_hospitD, MARGIN = 1, FUN=fcode, var="code")
  #merge with other cleaned hospital data 
  clinical_hospitD_clean = merge(outcomes_hospitD[,c("code","Manejo","CareLevelFinal","WHOFinal4cat","WHORevisedFinal","CareLevel12hr","WHO12hr4cat","WHORevised12hr")], 
                                 temp, by="code", all.y=T)
  tabletot(clinical_hospitD_clean, "ClasificacionFinal", "WHOFinal4cat", useNA="always") #sanity check.  Only one mismatch.
  tabletot(clinical_hospitD_clean, "WHO24hr4cat", "WHO12hr4cat", useNA="always") #not a totally terrible match (but misaligned in both directions)
  tabletot(clinical_hospitD_clean, "WHO24hr4cat", "ClasificacionFinal") #sample size for predicting severe dengue
  clinical_hospitD_clean$WHOFinal4cat = NULL #remove to avoid confusion (will produce it later using ClasificacionFinal and DENV indicator)
  
  # WHO revised classification for first 24 hr data
  #If [Dolor_Abdominal] = Yes or [Vomito]  = Yes or [Ascites, Plamsa leakage or  Edema] = Yes or [Mucosal hemorrhage] = Yes 
  #or [Lethargy] = Yes or [Irritability] = Yes or  [Hepatomegaly] =Yes  or [Decrease in platelets concurrent with increase in hematocrit] = Yes 
  #or [Platelets <100,000 and increased hematocrit] = Yes, then [Dengue with Warning SIgns] = Yes
  


  ###############################################################################
  ################# Data for cohort samples #####################################
  ###############################################################################
  
  #basic_cohortD = read.delim(paste(clinical_inputsDir,"Basic info_cohort study.txt", sep=""), header=TRUE, nrows=100) #old
  basic_cohortD = read.delim(paste(clinical_inputsDir,"Lab_cohort_noninvasive samples.txt", sep=""), header=TRUE, nrows=100) #new (more comprehensive)
  names(basic_cohortD)[names(basic_cohortD)=="ClasificacionPrimerDia"] <- "WHO24hr4cat"
  names(basic_cohortD)[names(basic_cohortD)=="Res_Final"] <- "DENV"
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  #need to take the characters up till the first dot, and then reformat
  formatID = function(x) sprintf("ID%04d",as.numeric(strsplit(as.character(x["Sample"]), "\\.")[[1]][1:1]))
  newIDs_3 = apply(X=basic_cohortD, MARGIN = 1, FUN=formatID)
  basic_cohortD["code"] = newIDs_3
  
  #add an indicator for type of study (will need this later when we do the merge)
  basic_cohortD["Study"] = "Cohort"
  
  
  ## add clinical/lab info
  #clinic_cohortD = read.delim(paste(clinical_inputsDir,"Clinical data_cohort study.txt", sep=""), header=TRUE, nrows=100) #old
  clinic_cohortD = read.delim(paste(clinical_inputsDir,"Clinical_cohort_noninvasive samples.txt", sep=""), header=TRUE, nrows=100) #new (more comprehensive)
  
  newIDs_4 = apply(X=clinic_cohortD, MARGIN = 1, FUN=formatID)
  clinic_cohortD["code"] = newIDs_4
  
  #calculate age of patient as of sampling date (first reformat date variables)
  clinic_cohortD$FTM = as.Date(clinic_cohortD$FECHA_TOMA, format = "%d/%m/%Y")
  clinic_cohortD$Fecha_Nac = as.Date(clinic_cohortD$fechana, format = "%m/%d/%Y")
  clinic_cohortD$FIS = as.Date(clinic_cohortD$FISint, format = "%m/%d/%Y")
  clinic_cohortD$FIF = as.Date(clinic_cohortD$FIFiebre, format = "%m/%d/%Y")
  clinic_cohortD$age = as.numeric(clinic_cohortD$FTM - clinic_cohortD$Fecha_Nac)/365.25
  #calculate days from symptom onset to sampling date (changed from days since symptom onset on 11-25-2014)
  clinic_cohortD$DaysSick = as.numeric(clinic_cohortD$FTM - clinic_cohortD$FIF) + 1
  #check calculation
  #View(clinic_cohortD[,c("FTM","Fecha_Nac","FIS","age","DaysSick", "Edad_toma")])
  
  clinical_cohortD <- merge(basic_cohortD[,c("code","Study", "IR","PCR",
                                             "WHO24hr4cat","ClasificacionFinal", "DENV")], 
                            clinic_cohortD, by=c("code"), all=TRUE)
  
  #keep only relevant variables and 
  #rename clinical/lab variables to match with names in hospital data
  nonclinic_keepervars = c("code","Study","FTM","age","DaysSick", 
                           "WHO24hr4cat","ClasificacionFinal", "DENV","Fecha_Nac")
  oldnames_nonLCMS = as.character(c(clinic_varsC_char, nonclinic_keepervars))
  clinic_varsCnew = clinic_varsD[which(clinic_varsD$in.cohort.data==1 & 
                                      (clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1))
                              ,"variable.name.in.final.data"]
  clinic_varsCnew_char = sapply(clinic_varsCnew, as.character)
  newnames_nonLCMS = as.character(c(clinic_varsCnew_char, nonclinic_keepervars))
  colnames(clinical_cohortD)
  clinical_cohortD_clean = clinical_cohortD[,oldnames_nonLCMS]
  colnames(clinical_cohortD_clean) = newnames_nonLCMS #rename columns
  
  #make the PCR variable a factor variable 
  #ND patients get PCR value of 0 and NAs will be given to indeterminants
  clinical_cohortD_clean[which(clinical_cohortD_clean$DENV=="Negativo"),"temp"]=0
  clinical_cohortD_clean[which(clinical_cohortD_clean$PCR==1),"temp"]=1
  clinical_cohortD_clean[which(clinical_cohortD_clean$PCR==2),"temp"]=2
  clinical_cohortD_clean[which(clinical_cohortD_clean$PCR==3),"temp"]=3
  clinical_cohortD_clean[which(clinical_cohortD_clean$PCR==4),"temp"]=4
  clinical_cohortD_clean$PCR = NULL
  clinical_cohortD_clean$PCR = factor(x = as.character(clinical_cohortD_clean$temp))
  clinical_cohortD_clean$temp = NULL
  table(clinical_cohortD_clean$PCR)
  class(clinical_cohortD_clean$PCR)
  
  #turn Torniquete into factor (rather than numeric) variable
  clinical_cohortD_clean$temp = factor(x = as.character(clinical_cohortD_clean$Torniquete))
  clinical_cohortD_clean["Torniquete"] = NULL
  clinical_cohortD_clean["Torniquete"] = clinical_cohortD_clean$temp
  clinical_cohortD_clean$temp = NULL
  table(clinical_cohortD_clean$Torniquete)
  class(clinical_cohortD_clean$Torniquete)
  
  #these variables should be coded as 0/1 --- replace 2s with NAs 
  #Lionel says 2=unknown in all of the clinical databases at the health center (2-24-2014 email)
  clinical_cohortD_clean[which(clinical_cohortD_clean$Hipermenorrea==2),"Hipermenorrea"]=NA
  clinical_cohortD_clean[which(clinical_cohortD_clean$Hemoconcentracion==2),"Hemoconcentracion"]=NA
  table(clinical_cohortD_clean$Hipermenorrea)
  table(clinical_cohortD_clean$Hemoconcentracion)
  
  temp = create_binary_variables(clinical_cohortD_clean) 
  
  #do we have varibles needed to determine DHF/DSS in first 24 hrs?  No.
  get_mismatches(colnames(clinical_cohortD_clean),DxDHF_varsH,"in_data","needed") #still needed: "Proteina_tot","Albumina","Purpura","Hematoma","Venopuncion","Vaginal","Hemoptisis","Derrame_","Ascitis","Edema_"
  get_mismatches(colnames(clinical_cohortD_clean),DxDSS_varsH,"in_data","needed") #still needed: "Sudoracion", "Escalofrio"
  
  #add additional outcome indicators (new WHO classifications and level of care)
  all_symptomatic = readWorksheetFromFile(paste(clinical_inputsDir,"all_symptomatic_2004_2014_classification_v4.xlsx", sep=""), sheet=1) #609 obs (all the DENV positive patients)
  all_symptomatic["code"] = apply(X=all_symptomatic, MARGIN = 1, FUN=fcode, var="Codigo")
  table(all_symptomatic$OMSTradicional, all_symptomatic$Res_Final) #only 44 patients diagnosed as "DHF/DSS"
  table(all_symptomatic$OMSTradicional, all_symptomatic$Hospitalizado) #all DHF and DSS patients were hospitalized
  table(all_symptomatic[which(all_symptomatic$Hospitalizado=="N"), "Res_Final"], useNA="always") #additional DENV patients in cohort but not hospit
  #level of care
  all_symptomatic$CareLevelFinal = substr(all_symptomatic$nivelcuidado, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
  all_symptomatic[which(all_symptomatic$Hospital=="N"),"CareLevelFinal"] = "1" #Lionel says non-hospitalized patients should all be considered level of care 1
  #revised WHO final diagnosis
  names(all_symptomatic)[names(all_symptomatic)=="OMSNueva"] <- "WHORevisedFinal" #rename variable
  all_symptomatic[which(is.na(all_symptomatic$WHORevisedFinal)),"WHORevisedFinal"] = "Caso sin datos de alarma" #if NA, assume lowest category (should confirm this)
  all_symptomatic[which(all_symptomatic$WHORevisedFinal==""),"WHORevisedFinal"] = "Caso sin datos de alarma" #if empty, assume lowest category
  clinical_cohortD_clean = merge(temp, all_symptomatic[c("code","WHORevisedFinal","CareLevelFinal")], by="code", all.x=T)
  
  #examine potential of cohort data (count of DENV Negativo and Positivo)
  casos = read.delim(paste(clinical_inputsDir,"suspected_dengue_cohorte_2004_2015.txt", sep=""), header=TRUE, nrows=10000) #5372 obs
  casos["code"] = apply(X=casos, MARGIN = 1, FUN=fcode, var="Codigo")
  table(casos$ResFinal.Dengue) #Lionel says to only consider Negativo (3948) and Positivo (610) patients
  
  
  ###############################################################################
  #### Merge hospital study data with cohort data and do additional cleaning ####
  ###############################################################################
  
  #Merge using both personID and study type since person ID is not unique across studies (I verified these are diff people)
  #smartbind allows us to have some variables only present in 1 of the datasets (unlike merge function)
  #(rbind insists on each dataframe containing the same variables)
  clinical_comboD = smartbind(clinical_cohortD_clean, clinical_hospitD_clean)
  
  #give categorical variables more meaningful labels
  clinical_comboD$Pulso <- factor(clinical_comboD$Pulso, levels = c("F","M","R","N"),
                                  labels = c("strong", "moderate", "rapid","not palpable")) 


  #clunky code due to weird stuff with ClasificacionFinal being a factor variable
  #If DENV="negative" then diagnostico="ND"
  clinical_comboD[which(clinical_comboD$DENV=="Negativo"),"WHOFinal4cat"] = "ND"   
  #If DENV="positivo" then diagnostico=value of ClasificacionFinal.  
  clinical_comboD[which(clinical_comboD$DENV=="Positivo" & clinical_comboD$ClasificacionFinal=="DF"),"WHOFinal4cat"] = "DF"
  clinical_comboD[which(clinical_comboD$DENV=="Positivo" & clinical_comboD$ClasificacionFinal=="DHF"),"WHOFinal4cat"] = "DHF"
  clinical_comboD[which(clinical_comboD$DENV=="Positivo" & clinical_comboD$ClasificacionFinal=="DSS"),"WHOFinal4cat"] = "DSS"
  #Else, diagnostico is missing (NA).
  clinical_comboD$WHOFinal4cat_test = factor(clinical_comboD$WHOFinal4cat)
  clinical_comboD$WHOFinal4cat = NULL
  clinical_comboD$WHOFinal4cat = clinical_comboD$WHOFinal4cat_test
  clinical_comboD$WHOFinal4cat_test = NULL
  class(clinical_comboD$WHOFinal4cat)
  table(clinical_comboD$WHOFinal4cat)
  
  #we do not care to separate DHF from DSS in analysis -- combine into one 
  clinical_comboD[which(clinical_comboD$WHOFinal4cat=="DF"), "WHOFinal3cat"] = "DF"
  clinical_comboD[which(clinical_comboD$WHOFinal4cat=="DHF" |  clinical_comboD$WHOFinal4cat=="DSS"), "WHOFinal3cat"] = "DHF_DSS"
  clinical_comboD[which(clinical_comboD$WHOFinal4cat=="ND"), "WHOFinal3cat"] = "ND"
  clinical_comboD$WHOFinal3cat_test = factor(clinical_comboD$WHOFinal3cat)
  clinical_comboD$WHOFinal3cat = NULL
  clinical_comboD$WHOFinal3cat = clinical_comboD$WHOFinal3cat_test
  clinical_comboD$WHOFinal3cat_test = NULL
  
  #All of the cohort samples which have DENV = “Infeccion reciente” also have IR=”primary”.  
  #This seems contradictory to me, as IR should indicate whether the current infection is the first Dengue infection that the patient has had.  
  #(In the hospital data, there are also 10 observations with this same combo.)
  #Lionel says, “the immune response is computed even if the DENV is not positive but it is only meaningful for positives” (email on 2-24-2014)
  #Therefore, I will turn the IR values to NA if DENV is not "Positivo"
  #Also, replace indeterminados with NAs (I don't see a meaningful difference)
  clinical_comboD[which(clinical_comboD$DENV!="Positivo"),"IR"] = NA
  clinical_comboD[which(clinical_comboD$IR=="I"),"IR"] = NA
  
  #create boolean output for ND vs. DEN analysis.  Note that factor gives problem due to inconsistency with level definitions.
  clinical_comboD$is.DEN = (clinical_comboD$WHOFinal3cat=="DF" | clinical_comboD$WHOFinal3cat=="DHF_DSS")
  clinical_comboD[which(is.na(clinical_comboD$WHOFinal3cat)),"is.DEN"] = NA
  
  #create boolean output for NF vs. DHF/DSS analysis with will be NA for ND
  clinical_comboD$is.DHF_DSS = (clinical_comboD$WHOFinal3cat=="DHF_DSS")
  clinical_comboD[which(is.na(clinical_comboD$WHOFinal3cat) | clinical_comboD$WHOFinal3cat=="ND"),"is.DHF_DSS"] = NA

  #convert all yes/no questions into factors (regression algorithms may run suboptimal if coded as numeric)
  YesNo_varsCH_char = c(YesNo_varsH_char, "NAUSEAS","HIPOTERM","CONGNASAL")
  for (name in YesNo_varsCH_char) {
    #note that underlying values will still be 1, 2 (!) and I do not see how to change this
    clinical_comboD[,name] = factor(clinical_comboD[,name], levels=c(0,1), labels=c("No","Yes"))
  } 
  #Sexo is already a factor, but we want to make its underlying values 0/1 rather than 1/2
  clinical_comboD[,"Sexo"] = as.numeric(clinical_comboD$Sexo)-1 # now 0 is female, 1 is male
  clinical_comboD[,"Sexo"] = factor(clinical_comboD[,"Sexo"], levels=c(0,1), labels = c("female","male"))
  
  #Nobody had a "not palpable" pulse and keeping this level in there causes problems with RF, so remove it here
  clinical_comboD$Pulso = factor(as.numeric(clinical_comboD$Pulso), levels=c("1","2","3"), labels=c("strong","moderate","rapid"))

  #Nobody had an "indeterminant" IR and keeping this level in there causes problems with RF, so remove it here
  clinical_comboD$IR = factor(as.numeric(clinical_comboD$IR), levels=c("2","3"), labels=c("primary","secondary"))  
  
  #this seems like an obvious data entry mistake to me
  clinical_comboD[which(clinical_comboD$Engrosamiento_mm>200),"Engrosamiento_mm"]=NA 
  
  #add indicator of presence in LCMS data
  load(paste(outputsDir,"IDs_in_resp_all.RData", sep="")) #loads IDs_in_resp_all
  clinical_comboD = merge(IDs_in_resp_all, clinical_comboD, by=c("code","Study"), all.y=TRUE)
  clinical_comboD$serum.label <- "Indicates LCMS serum batch number"
  
  #see how many DHF/DSS diagnoses there are to predict
  tabletot(clinical_comboD[which(clinical_comboD$Study=="Hospital"),], "WHO24hr4cat", "WHOFinal4cat") # 118 to predict (not including all cohort patients)
  tabletot(clinical_comboD[which(clinical_comboD$Study=="Hospital"),], "WHO12hr4cat", "WHOFinal4cat") # 136 to predict (not including all cohort patients) 
  tabletot(clinical_comboD[which(clinical_comboD$Study=="Cohort"),], "WHO24hr4cat", "WHOFinal4cat")   # 2 to predict among cohort subset

  comment(clinical_comboD) <- "This data contains all clinical/lab/demographic info for data from Nicaragua.  
  Observations with unknown final diagnosis are excluded.  Otherwise, missing values are included and untampered with."
  #note: cat(comment(clinical_full_clean)) will display this data comment.
    
  return(clinical_comboD)
}

