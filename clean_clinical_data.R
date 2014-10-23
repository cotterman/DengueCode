
####################################################################
##################### Clean the clinical data ######################
####################################################################

# function produces clean clinical data
clean_clinical_data = function(clinic_varsD){
  
  ####### List of clinical/lab variables to consider for analysis ########
  
  #list names of variables from hospital data to include in either analysis (these are what will become the standardized names)
  clinic_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | 
                                      clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),"Variable.Name.in.Hospital.data"]
  clinic_varsH_char = sapply(clinic_varsH, as.character) #convert factor to character
  
  #corresponding names in cohort data (these names will be changed to match the hospital data variable names)
  clinic_varsC = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | 
                                      clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),"Variable.Name.in.Cohort.data"]
  clinic_varsC_char = sapply(clinic_varsC, as.character)
  
  
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
  #in future, we convert all yes/no variables to factors
  YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
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
  
  other_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output!="Yes, No"),
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
  #calculate days from symptom onset to sampling date
  clinical_hospitD$DaysSick = as.numeric(clinical_hospitD$FTM - clinical_hospitD$FIS) + 1
  #check calculation
  #View(clinical_hospitD[,c("FTM","Fecha_Nac","FIS","age","DaysSick")])
  
  #keep only relevant varaibles
  nonclinic_keepervars = c("code","Study","FTM","age","DaysSick", 
                           "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final")
  namesH_nonLCMS = as.character(c(clinic_varsH_char, nonclinic_keepervars))
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
  #calculate days from symptom onset to sampling date
  clinic_cohortD$DaysSick = as.numeric(clinic_cohortD$FTM - clinic_cohortD$FIS) + 1
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
  newnames_nonLCMS = as.character(c(clinic_varsH_char, nonclinic_keepervars))
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
  #smartbind allows us, in the future, to have some variables only present in 1 of the datasets
  #(rbind insists on each dataframe containing the same variables)
  clinical_comboD = rbind(clinical_cohortD_clean, clinical_hospitD_clean)
  
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
  
  #create binary output for ND vs. DEN analysis
  clinical_comboD[which(clinical_comboD$DxFinal3cat=="ND"),"DEN_dum"] = 0  
  clinical_comboD[which(clinical_comboD$DxFinal3cat=="DF" | clinical_comboD$DxFinal3cat=="DHF_DSS"),"DEN_dum"] = 1
  #convert to factor in order to run regression routines that require outcome to be factor
  clinical_comboD$DEN_dum = factor(clinical_comboD$DEN_dum)
  table(clinical_comboD$DxFinal3cat, clinical_comboD$DEN_dum)
  
  #create binary output for NF vs. DHF/DSS analysis with will be NA for ND
  clinical_comboD[which(clinical_comboD$DxFinal3cat=="DF"),"DHF_DSS_dum"] = 0  
  clinical_comboD[which(clinical_comboD$DxFinal3cat=="DHF_DSS"),"DHF_DSS_dum"] = 1
  #convert to factor in order to run regression routines that require outcome to be factor
  clinical_comboD$DHF_DSS_dum = factor(clinical_comboD$DHF_DSS_dum)
  table(clinical_comboD$DxFinal3cat, clinical_comboD$DHF_DSS_dum, useNA="ifany")
  
  #convert all yes/no questions into factors (regression algorithms may run suboptimal if coded as numeric)
  for (name in YesNo_varsH_char) {
    clinical_comboD[,name] = factor(clinical_comboD[,name])
  } 
  
  return(clinical_comboD)
}
