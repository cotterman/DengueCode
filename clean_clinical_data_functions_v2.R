
####################################################################
##################### Clean the clinical data ######################
####################################################################

# function produces clean clinical data
clean_clinical_data = function(clinic_varsD){
  
  ####### List of clinical/lab variables to consider for analysis ########
  
  #names of variables from hospital data to include in either analysis (these are what will become the standardized names)
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
  ####################### Data from Mexico ######################################
  ###############################################################################  
  
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
  #estimate number of patients who end with with DHF but did not have DHF
    #when they appeared to clinic (so days 0 - 2 at most)
  table(temp[which(temp$DxFinal4cat=="DHF"),"DaysSick"])
  
  
  ###############################################################################
  ##################### Data from the hospital study ############################
  ###############################################################################
  # This should contain data from the first batch and from some of the 3rd/4th batches
  
  clinical_hospitD = read.delim(paste(clinical_inputsDir,"Clinical data_hospital study.txt", sep=""), header=TRUE, nrows=2000)
  
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  fcode = function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
  newIDs_4 = apply(X=clinical_hospitD, MARGIN = 1, FUN=fcode)
  clinical_hospitD["code"] = newIDs_4
  
  #code the yes/no variables to take numeric values 1/0 (to be consistent with cohort data)
  #in future, we convert all yes/no variables to logical vectors (booleans)
  YesNo_varsH = clinic_varsD[which((clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1) & 
                                     clinic_varsD$Variable.Output=="Yes, No" & clinic_varsD$in.hospit.data==1),
                             "Variable.Name.in.Hospital.data"]
  YesNo_varsH_char = sapply(YesNo_varsH, as.character)
  
  for (name in YesNo_varsH_char) {
    #print(name)
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
  clinical_hospitD[which(clinical_hospitD$Res_Final=="Negativo"),"PCR"]=0
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
                           "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final")
  namesH_nonLCMS = as.character(c(clinic_varsH_char, nonclinic_keepervars))
  #drop the one guy who is missing "DaysSick" since this variable is too important
  clinical_hospitD_clean = clinical_hospitD[,namesH_nonLCMS]
 
  
  ###############################################################################
  ################# Data for cohort samples #####################################
  ###############################################################################
  
  #basic_cohortD = read.delim(paste(clinical_inputsDir,"Basic info_cohort study.txt", sep=""), header=TRUE, nrows=100) #old
  basic_cohortD = read.delim(paste(clinical_inputsDir,"Lab_cohort_noninvasive samples.txt", sep=""), header=TRUE, nrows=100) #new (more comprehensive)
  
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
                                             "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final")], 
                            clinic_cohortD, by=c("code"), all=TRUE)
  
  #keep only relevant variables and 
  #rename clinical/lab variables to match with names in hospital data
  nonclinic_keepervars = c("code","Study","FTM","age","DaysSick", 
                           "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final")
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
  clinical_cohortD_clean[which(clinical_cohortD_clean$Res_Final=="Negativo"),"temp"]=0
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
  #Lionel says 2=unknonw in all of the clinical databases at the health center (2-24-2014 email)
  clinical_cohortD_clean[which(clinical_cohortD_clean$Hipermenorrea==2),"Hipermenorrea"]=NA
  clinical_cohortD_clean[which(clinical_cohortD_clean$Hemoconcentracion==2),"Hemoconcentracion"]=NA
  table(clinical_cohortD_clean$Hipermenorrea)
  table(clinical_cohortD_clean$Hemoconcentracion)
  
  
  ###############################################################################
  #### Merge hospital study data with cohort data and do additional cleaning ####
  ###############################################################################
  
  #Merge using both personID and study type since person ID is not unique across studies
  #smartbind allows us to have some variables only present in 1 of the datasets (unlike merge function)
  #(rbind insists on each dataframe containing the same variables)
  clinical_comboD = smartbind(clinical_cohortD_clean, clinical_hospitD_clean)
  
  #give categorical variables more meaningful labels
  clinical_comboD$Pulso <- factor(clinical_comboD$Pulso, levels = c("F","M","R","N"),
                                  labels = c("strong", "moderate", "rapid","not palpable")) 
  
  #Create final diagnosis variable 
  #clunky code due to weird stuff with ClasificacionFinal being a factor variable
  #If Res_Final="negative" then diagnostico="ND"
  clinical_comboD[which(clinical_comboD$Res_Final=="Negativo"),"DxFinal4cat"] = "ND"   
  #If Res_Final="positivo" then diagnostico=value of ClasificacionFinal.  
  clinical_comboD[which(clinical_comboD$Res_Final=="Positivo" & clinical_comboD$ClasificacionFinal=="DF"),"DxFinal4cat"] = "DF"
  clinical_comboD[which(clinical_comboD$Res_Final=="Positivo" & clinical_comboD$ClasificacionFinal=="DHF"),"DxFinal4cat"] = "DHF"
  clinical_comboD[which(clinical_comboD$Res_Final=="Positivo" & clinical_comboD$ClasificacionFinal=="DSS"),"DxFinal4cat"] = "DSS"
  #Else, diagnostico is missing (NA).
  clinical_comboD$DxFinal4cat_test = factor(clinical_comboD$DxFinal4cat)
  clinical_comboD$DxFinal4cat = NULL
  clinical_comboD$DxFinal4cat = clinical_comboD$DxFinal4cat_test
  clinical_comboD$DxFinal4cat_test = NULL
  class(clinical_comboD$DxFinal4cat)
  table(clinical_comboD$DxFinal4cat)
  
  #we do not care to separate DHF from DSS in analysis -- combine into one 
  clinical_comboD[which(clinical_comboD$DxFinal4cat=="DF"), "DxFinal3cat"] = "DF"
  clinical_comboD[which(clinical_comboD$DxFinal4cat=="DHF" |  clinical_comboD$DxFinal4cat=="DSS"), "DxFinal3cat"] = "DHF_DSS"
  clinical_comboD[which(clinical_comboD$DxFinal4cat=="ND"), "DxFinal3cat"] = "ND"
  clinical_comboD$DxFinal3cat_test = factor(clinical_comboD$DxFinal3cat)
  clinical_comboD$DxFinal3cat = NULL
  clinical_comboD$DxFinal3cat = clinical_comboD$DxFinal3cat_test
  clinical_comboD$DxFinal3cat_test = NULL
  
  #All of the cohort samples which have Res_Final = “Infeccion reciente” also have IR=”primary”.  
  #This seems contradictory to me, as IR should indicate whether the current infection is the first Dengue infection that the patient has had.  
  #(In the hospital data, there are also 10 observations with this same combo.)
  #Lionel says, “the immune response is computed even if the Res_Final is not positive but it is only meaningful for positives” (email on 2-24-2014)
  #Therefore, I will turn the IR values to NA if Res_Final is not "Positivo"
  #Also, replace indeterminados with NAs (I don't see a meaningful difference)
  clinical_comboD[which(clinical_comboD$Res_Final!="Positivo"),"IR"] = NA
  clinical_comboD[which(clinical_comboD$IR=="I"),"IR"] = NA
  
  #create boolean output for ND vs. DEN analysis.  Note that factor gives problem due to inconsistency with level definitions.
  clinical_comboD$is.DEN = (clinical_comboD$DxFinal3cat=="DF" | clinical_comboD$DxFinal3cat=="DHF_DSS")
  clinical_comboD[which(is.na(clinical_comboD$DxFinal3cat)),"is.DEN"] = NA
  
  #create boolean output for NF vs. DHF/DSS analysis with will be NA for ND
  clinical_comboD$is.DHF_DSS = (clinical_comboD$DxFinal3cat=="DHF_DSS")
  clinical_comboD[which(is.na(clinical_comboD$DxFinal3cat) | clinical_comboD$DxFinal3cat=="ND"),"is.DHF_DSS"] = NA

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
  
  return(clinical_comboD)
}


#function to obtain list of clinical variables based on user-supplied criteria
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

