###############################################################################
### Much of this code used to be in "clean_clinical_data_functions_v2.R"
### But I have separated out and present here the code 
### that is exploratory (comparing file contents etc.) 
### so that "clean_clinical_data_functions_v2.R" focuses on building 
### data for analysis without code unused for analysis
###############################################################################

HospitD_cleaning = function(df, clinic_varsD, clinic_varsH_char, YesNo_varsH_char){  
  
  df$temp = NA
  for (name in YesNo_varsH_char) {
    df[which(df[name]=="Si"),"temp"]=1
    df[which(df[name]=="No"),"temp"]=0
    class(df$temp)
    #drop old version and rename new version
    df[name] = NULL
    df[name] = df$temp
    df$temp = NULL
  }
  
  other_varsH = clinic_varsD[which((clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1)
                                   & clinic_varsD$Variable.Output!="Yes, No"),
                             c("Variable.Name.in.Hospital.data","Variable.Output")]
  other_varsH #see what is left
  other_varsH_char = sapply(other_varsH["Variable.Name.in.Hospital.data"], as.character)
  #make sure these variables are compatible between datasets
  #for (name in other_varsH_char) {
  #  print(name)
  #  print(summary(df[name]))
  #  print(summary(clinical_cohort_subset[name]))
  #}
  
  #make the PCR variable a factor variable (else it will mismatch with cohort data and generate error when merged)
  #also, give ND patients a PCR value of 0.  (NAs will be given to indeterminants)
  df[which(df$DENV=="Negativo"),"PCR"]=0
  df$PCR = factor(x = as.character(df$PCR))
  table(df$PCR)
  
  #make IR variable consistent with cohort data
  df[which(df$IR=="Secundario"),"new_IR"]="S"
  df[which(df$IR=="Primario"),"new_IR"]="P"
  df[which(df$IR=="Indeterminado"),"new_IR"]="I"
  df["IR"] = NULL
  df$IR = factor(x = as.character(df$new_IR))
  table(df$IR)
  class(df$IR)
  
  #turn Torniquete into factor (rather than numeric) variable
  df$temp = factor(x = as.character(df$Torniquete))
  df["Torniquete"] = NULL
  df["Torniquete"] = df$temp
  df$temp = NULL
  table(df$Torniquete)
  class(df$Torniquete)
  
  #add an indicator for type of study (will need this later when we merge)
  df["Study"] = "Hospital"
  
  #calculate age of patient as of sampling date (first reformat variables into dates)
  df$FTM = as.Date(df$FTM, format = "%m/%d/%Y")
  df$Fecha_Nac = as.Date(df$Fecha_Nac, format = "%m/%d/%Y")
  df$FIS = as.Date(df$FIS, format = "%m/%d/%Y")
  df$FIF = as.Date(df$FIF, format = "%m/%d/%Y")
  df$age = as.numeric(df$FTM - df$Fecha_Nac)/365.25
  #calculate days from fever onset to sampling date (changed from days since symptom onset on 11-25-2014)
  df$DaysSick = as.numeric(df$FTM - df$FIF) + 1
  #use FIS if FIF is missing
  df[which(is.na(df$FIF)),]$DaysSick = 
    as.numeric(df[which(is.na(df$FIF)),]$FTM - df[which(is.na(df$FIF)),]$FIS) + 1
  #check calculation
  #View(df[,c("FTM","Fecha_Nac","FIS","age","DaysSick")])
  
  #keep only relevant variables
  nonclinic_keepervars = c("code","Study","FTM","age","DaysSick", 
                           "ClasificacionPrimerDia","ClasificacionFinal", "Res_Final","Fecha_Nac")
  namesH_nonLCMS = as.character(c(clinic_varsH_char, nonclinic_keepervars))
  #drop the one guy who is missing "DaysSick" since this variable is too important
  df_prelim = df[,namesH_nonLCMS]
  #rename some vars
  names(df_prelim)[names(df_prelim)=="Res_Final"] <- "DENV"
  
  ### Re-create initial diagnosis variable ###
  DxDHF_varsH = list("Torniquete","Proteina_tot","Albumina","Plaquetas",
                     "PetequiasEspontaneas","Equimosis","Purpura","Hematoma",
                     "Venopuncion","Hipermenorrea","Vaginal","Subconjuntival","Hemoptisis","Nariz","Encias",
                     "Melena","Hematemesis",
                     "Derrame_","Ascitis","Edema_","Hemoconcentracion","HematocritoElev")
  get_mismatches(colnames(df_prelim),DxDHF_varsH,"in_data","needed") #we have all needed vars
  DxDSS_varsH = list("Presion_Arterial_Sist","Presion_Arterial_Dias", 
                     "ExtremidadesFrias","Palidez","Sudoracion","Escalofrio","Llenado_Capilar")
  get_mismatches(colnames(df_prelim),DxDSS_varsH,"in_data","needed") #we have all needed vars
  
  #DHF
  #creates is.torniquete20plus, is.thrombocytopenia, is.Estrechamiento, is.hypotension (among others)
  temp = create_binary_variables(df_prelim) 
  
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
  # Menorrhagia, vaginal bleeding(=Vaginal), subconjunctival bleeding(=Subconjuntival), hemoptysis(=Hemoptisis), epistaxis(=nariz), or gingivorrhagia
  temp$Hem2 = apply(temp[,c("Venopuncion","Hipermenorrea","Vaginal","Subconjuntival","Hemoptisis","Nariz","Encias")], 1, function(x) sum(x, na.rm=T)>0) 
  temp$Hem3 = apply(temp[,c("Melena","Hematemesis")], 1, function(x) sum(x, na.rm=T)>0) #Melena or hematemesis (straight-foward)
  #Pleural Effusion is Derrame_ #check this?? PleuralEffusionLeftRx,  PleuralEffusionRightRx, PulmonaryEdemaRx, PleuralEffusionClinical, PleuralEffusionRightUS or PleuralEffusionLeftUS
  #make sure the "Edema_" variable includes "Edema, periorbital, facial, limbs, hydrocele, generalized and fovea"
  #If (Fever = Yes) AND  ([Tourniquet test] = Yes or [HEM1] = Yes or [HEM2] = Yes or [HEM3] = Yes) AND (Thrombocyto-penia = Yes) 
  #  AND ([Hemoconcentration] = Yes or [ElevHematocrit] = Yes  or [PlerualEffusion] = Yes or [Ascites] = Yes or [Hypoprotein-emia] = Yes or [Hypoalbumin-emia] =Yes or [Edema] = Yes), then [DHF] = Yes
  temp$Leakage = apply(temp[,c("Derrame_","Ascitis","Edema_","is.Hypoproteinemia","is.Hypoalbuminemia","Hemoconcentracion","HematocritoElev")], 1, function(x) sum(x, na.rm=T)>0)
  #I took out is.fever since it can be biphasic (patient could have prior history of fever)
  #when is.thrombocytopenia is NA, can give a 0 but will never give a 1 (gives NA instead if other criteria are satisfied and is.thrombocytopenia is NA -- reasonable)
  temp$DHF = with(temp, ifelse((is.torniquete20plus | Hem1 | Hem2 | Hem3) & is.thrombocytopenia & Leakage, 1, 0 ))
  #View(temp[,c("is.torniquete20plus","Hem1","Hem2","Hem3","is.thrombocytopenia","Leakage","DHF")]) 
  
  #DSS (is.Estrechamiento, is.hypotension, is.pulse_danger, is.coldclammy created in "create_binary_variables")
  temp$DSS = with(temp, ifelse( DHF==1 & (is.Estrechamiento | is.hypotension) & (is.pulse_danger | is.coldclammy | Llenado_Capilar), 1, 0 ))
  
  return(temp)
}



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

#code the yes/no variables to take numeric values 1/0 (to be consistent with cohort data)
#in future, we convert all yes/no variables to logical vectors (booleans)
YesNo_varsH = clinic_varsD[which((clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1) & 
                                   clinic_varsD$Variable.Output=="Yes, No" & clinic_varsD$in.hospit.data==1),
                           "Variable.Name.in.Hospital.data"]
YesNo_varsH_char = sapply(YesNo_varsH, as.character)



###############################################################################
##################### Data from the hospital study ############################
###############################################################################

# This should contain "first 12 hr" data from the hospital
hospit12hr_raw = read.delim(paste(clinical_inputsDir,"Clinical_first_12hrs_hospital_study.txt", sep=""), header=TRUE, nrows=2000) #1771 obs
hospit12hr_raw["code"] = apply(X=hospit12hr_raw, MARGIN = 1, FUN=fcode, var="Code")
hospit12hr_raw["Epigastralgia"] = NA #will prevent error when this variable is expected

# This should contain "first 24 hr" data from the first batch and from some of the 3rd/4th batches
hospit24hr_raw = read.delim(paste(clinical_inputsDir,"Clinical data_hospital study.txt", sep=""), header=TRUE, nrows=2000) #1650 obs
hospit24hr_raw["code2"] = apply(X=hospit24hr_raw, MARGIN = 1, FUN=fcode, var="code")
hospit24hr_raw["code"] = hospit24hr_raw["code2"]

# check out differences in variable content between 12hr and 24hr data
get_mismatches(colnames(hospit12hr_raw),colnames(hospit24hr_raw),"in_12hr","in_24hr") 
#perfect match, except 12hr is missing Age, DaysSick, Epigastralgia.  Age and DaysSick are created later on, so no problem.

#work with 12-hr dta
clinical_hospit12hr_prelim = HospitD_cleaning(hospit12hr_raw, clinic_varsD, clinic_varsH_char, YesNo_varsH_char)
names(clinical_hospit12hr_prelim)[names(clinical_hospit12hr_prelim)=="ClasificacionPrimerDia"] <- "WHO12h_given"
#assume for now all patients in data satisfy the DF (suspected dengue) criteria
clinical_hospit12hr_prelim$WHO12hr4cat = ifelse(clinical_hospit12hr_prelim$DSS==1,"DSS",ifelse(clinical_hospit12hr_prelim$DHF==1,"DHF","DF")) 
table(clinical_hospit12hr_prelim$WHO12hr4cat, clinical_hospit12hr_prelim$WHO12h_given, useNA="always") 
#pretty good -- only 7 that I have as DHF but you dont and 2 vice versa.
#also, I have 5 with NA that you say are DF.  (probably related to missing values?)

#work with 24-hr data
clinical_hospit24hr_prelim = HospitD_cleaning(hospit24hr_raw, clinic_varsD, clinic_varsH_char, YesNo_varsH_char)
names(clinical_hospit24hr_prelim)[names(clinical_hospit24hr_prelim)=="ClasificacionPrimerDia"] <- "WHO24h_given"
#assume for now all patients in data satisfy the DF (suspected dengue) criteria
clinical_hospit24hr_prelim$WHO24hr4cat = ifelse(clinical_hospit24hr_prelim$DSS==1,"DSS",ifelse(clinical_hospit24hr_prelim$DHF==1,"DHF","DF")) 
table(clinical_hospit24hr_prelim$WHO24hr4cat, clinical_hospit24hr_prelim$WHO24h_given, useNA="always") #poor match.  use the one I calculated (WHO24hr4cat)

#merge 12-hr with 24-hr data
clinical_hospitD = merge(clinical_hospit12hr_prelim, clinical_hospit24hr_prelim, by="code", all=T)
#compare the 12-hr and 24-hr data
table(clinical_hospitD$ClasificacionFinal.x, clinical_hospitD$ClasificacionFinal.y, useNA="always") #sample size increase -- 118 added DF and 3 added DHF (final)
table(clinical_hospitD$WHO12hr4cat, clinical_hospitD$WHO24hr4cat, useNA="always")  #12-hr and 24-hr data have identical Dx values
vars_in_both = intersect(colnames(clinical_hospit24hr_prelim), colnames(clinical_hospit12hr_prelim))
vars_in_both_char = sapply(vars_in_both, as.character)
#compare values from 12 vs 24 data
for(var in vars_in_both_char){
  if(var!="code"){
    var.x = as.character(paste(var, ".x",sep=""))
    var.y = as.character(paste(var, ".y",sep=""))
    print(var.x)
    #print(mode(clinical_hospitD[, var.x])) 
    if(mode(clinical_hospitD[, var.x])=="numeric"){
      diff = (clinical_hospitD[, var.y] -  clinical_hospitD[, var.x])
      #positive value means higher in 24hr than 12hr (usually fine), neg values should happen only when small value mean poor health (e.g. platelets)
      print(table(diff))
    }
    if(mode(clinical_hospitD[, var.x])=="logical"){
      table(clinical_hospitD[, var.y], clinical_hospitD[, var.x])
    }
  }
}
#identical values with following exceptions:
#expected direction:
#Higado +2 (n=2), +3 (n=3)
#Rash +1 (n=1)
#unexpected direction: 
#Artralgia -1 (n=1)
#Mialgia -1 (n=1)
#Perdida_Apetito -1 (n=4)
#PetequiasEspontaneas -1 (n=5)
#Esplenomegalia_mm. -111 - -78 (n=4)
#Engrosamiento_mm -2.6 - -1.7 (n=4)
#Hepatomegalia_mm -125 - -98 (n=4)
#Albumina +.6 (n=1)
#unknown direction
#Tos -1 (n=1)
#TGO_ -6 (n=1)

## Additional outcome indicators (see if my DHF indicator matches with Lionel's 12 hr Dx indicator) ##
outcomes_hospitD = read.delim(paste(clinical_inputsDir,"Additional outcomes_hospital study.txt", sep=""), header=TRUE, nrows=2000)
outcomes_hospitD$CareLevel12h = substr(outcomes_hospitD$CareLevel12h, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
outcomes_hospitD$CareLevelFinal = substr(outcomes_hospitD$CareLevelFinal, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
outcomes_hospitD$code = apply(X=outcomes_hospitD, MARGIN = 1, FUN=fcode, var="code")

#merge with other cleaned hospital data 
clinical_hospitD_clean = merge(outcomes_hospitD[,c("code","Manejo","CareLevelFinal","WHOFinal4cat","WHORevisedFinal","CareLevel12hr","WHO12hr4cat","WHORevised12hr")], 
                               clinical_hospit24hr_prelim, by="code", all.y=T)
tabletot(clinical_hospitD_clean, "ClasificacionFinal", "WHOFinal4cat", useNA="always") #sanity check.  Only one mismatch.
tabletot(clinical_hospitD_clean, "WHO24hr4cat", "ClasificacionFinal") #sample size for predicting severe dengue
clinical_hospitD_clean$WHOFinal4cat = NULL #remove to avoid confusion (will produce it later using ClasificacionFinal and DENV indicator)

## Explore initial Dx of the 88 batch 1 samples ## 
hospitD_exper = merge(clinical_hospitD[,c("WHO24h_given", "WHO12h_given", "WHO12hr4cat", "WHO24hr4cat","DENV.x","code")], 
                      outcomes_hospitD[,c("WHOFinal4cat","code")], by="code")
hospitD_exper = merge(hospitD_exper, IDs_in_resp_all, by="code")
hospitD1_exper = hospitD_exper[which(hospitD_exper$serum==1),]
hospitD1_exper[which(hospitD1_exper$DENV.x=="Negativo"),"WHOFinal4cat"]="ND"
with(hospitD1_exper, table(DENV.x, WHOFinal4cat, useNA="always"))
with(hospitD1_exper[which(hospitD1_exper$WHOFinal4cat=="DHF" | hospitD1_exper$WHOFinal4cat=="DSS"),], table(WHO24hr4cat, WHO24h_given), useNA="always")
with(hospitD1_exper[which(hospitD1_exper$WHOFinal4cat=="DHF" | hospitD1_exper$WHOFinal4cat=="DSS"),], table(WHO12hr4cat, WHO12h_given), useNA="always")
with(hospitD1_exper, table(WHOFinal4cat, WHO12h_given), useNA="always") #13 of 30 had DF initially
with(hospitD1_exper, table(WHOFinal4cat, WHO12hr4cat), useNA="always")  #13 of 30 had DF initially
#try using just the outcomes dataset
outhospD_exper = merge(outcomes_hospitD, IDs_in_resp_all, by="code")
outhospD1_exper = outhospD_exper[which(outhospD_exper$serum==1),]
outhospD1_exper[which(outhospD1_exper$Res_Final=="Negativo"),"WHOFinal4cat"]="ND"
with(outhospD1_exper, table(WHOFinal4cat, Res_Final, useNA="always")) 
with(outhospD1_exper, table(WHOFinal4cat, WHO12hr4cat, useNA="always")) #these are consistent with Lionel numbers in email to CO
#12hr Dx across datasets provided by Nica
clinical_hospitD_clean = merge(outcomes_hospitD[,c("code","Manejo","CareLevelFinal","WHOFinal4cat","WHORevisedFinal","CareLevel12hr","WHO12hr4cat","WHORevised12hr")], 
                               clinical_hospit12hr_prelim[,c("WHO12h_given","code")], by="code", all.y=T)
with(clinical_hospitD_clean, table(WHO12hr4cat, WHO12h_given))

# WHO revised classification for first 24 hr data
#If [Dolor_Abdominal] = Yes or [Vomito]  = Yes or [Ascites, Plamsa leakage or  Edema] = Yes or [Mucosal hemorrhage] = Yes 
#or [Lethargy] = Yes or [Irritability] = Yes or  [Hepatomegaly] =Yes  or [Decrease in platelets concurrent with increase in hematocrit] = Yes 
#or [Platelets <100,000 and increased hematocrit] = Yes, then [Dengue with Warning SIgns] = Yes


###############################################################################
### Data for cohort samples (old data -- just for matching with LCMS data) ####
###############################################################################

# I now have data for virtually all cohort patients so no longer need this data subset

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
nonclinic_keepervars = c("code","Sample", "Study","FTM","age","DaysSick", 
                         "WHO24hr4cat","ClasificacionFinal", "DENV","Fecha_Nac")
oldnames_nonLCMS = as.character(c(clinic_varsC_char, nonclinic_keepervars))
clinic_varsCnew = clinic_varsD[which(clinic_varsD$in.cohort.data==1 & 
                                       (clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1))
                               ,"variable.name.in.final.data"]
clinic_varsCnew_char = sapply(clinic_varsCnew, as.character)
newnames_nonLCMS = as.character(c(clinic_varsCnew_char, nonclinic_keepervars))
colnames(clinical_cohortD)
clinical_cohort_subset = clinical_cohortD[,oldnames_nonLCMS]
colnames(clinical_cohort_subset) = newnames_nonLCMS #rename columns

#make the PCR variable a factor variable 
#ND patients get PCR value of 0 and NAs will be given to indeterminants
clinical_cohort_subset[which(clinical_cohort_subset$DENV=="Negativo"),"temp"]=0
clinical_cohort_subset[which(clinical_cohort_subset$PCR==1),"temp"]=1
clinical_cohort_subset[which(clinical_cohort_subset$PCR==2),"temp"]=2
clinical_cohort_subset[which(clinical_cohort_subset$PCR==3),"temp"]=3
clinical_cohort_subset[which(clinical_cohort_subset$PCR==4),"temp"]=4
clinical_cohort_subset$PCR = NULL
clinical_cohort_subset$PCR = factor(x = as.character(clinical_cohort_subset$temp))
clinical_cohort_subset$temp = NULL
table(clinical_cohort_subset$PCR)
class(clinical_cohort_subset$PCR)

#turn Torniquete into factor (rather than numeric) variable
clinical_cohort_subset$temp = factor(x = as.character(clinical_cohort_subset$Torniquete))
clinical_cohort_subset["Torniquete"] = NULL
clinical_cohort_subset["Torniquete"] = clinical_cohort_subset$temp
clinical_cohort_subset$temp = NULL
table(clinical_cohort_subset$Torniquete)
class(clinical_cohort_subset$Torniquete)

#these variables should be coded as 0/1 --- replace 2s with NAs 
#Lionel says 2=unknown in all of the clinical databases at the health center (2-24-2014 email)
clinical_cohort_subset[which(clinical_cohort_subset$Hipermenorrea==2),"Hipermenorrea"]=NA
clinical_cohort_subset[which(clinical_cohort_subset$Hemoconcentracion==2),"Hemoconcentracion"]=NA
table(clinical_cohort_subset$Hipermenorrea)
table(clinical_cohort_subset$Hemoconcentracion)

cohort_subset = create_binary_variables(clinical_cohort_subset) 

#do we have varibles needed to determine DHF/DSS in first 24 hrs?  No.
get_mismatches(colnames(clinical_cohort_subset),DxDHF_varsH,"in_data","needed") #still needed: "Proteina_tot","Albumina","Purpura","Hematoma","Venopuncion","Vaginal","Hemoptisis","Derrame_","Ascitis","Edema_"
get_mismatches(colnames(clinical_cohort_subset),DxDSS_varsH,"in_data","needed") #still needed: "Sudoracion", "Escalofrio"

#add additional outcome indicators (new WHO classifications and level of care)
all_symptomatic = readWorksheetFromFile(paste(clinical_inputsDir,"all_symptomatic_2004_2014_classification_v4.xlsx", sep=""), sheet=1) #609 obs (all the DENV positive patients)
all_symptomatic["code"] = apply(X=all_symptomatic, MARGIN = 1, FUN=fcode, var="Codigo")
table(all_symptomatic$OMSTradicional, all_symptomatic$Res_Final) #only 44 patients diagnosed as "DHF/DSS".  All obs in this data have Res_Final = positivo
table(all_symptomatic$OMSTradicional, all_symptomatic$Hospitalizado) #all DHF and DSS patients were hospitalized
table(all_symptomatic[which(all_symptomatic$Hospitalizado=="N"), "Res_Final"], useNA="always") #additional DENV patients in cohort but not hospit
#level of care
all_symptomatic$CareLevelFinal = substr(all_symptomatic$nivelcuidado, start=11, stop=12) #change "Categoria 1", "Categoria 2", "Categoria 3" to "1", "2", "3"
all_symptomatic[which(all_symptomatic$Hospital=="N"),"CareLevelFinal"] = "1" #Lionel says non-hospitalized patients should all be considered level of care 1
#revised WHO final diagnosis
names(all_symptomatic)[names(all_symptomatic)=="OMSNueva"] <- "WHORevisedFinal" #rename variable
all_symptomatic[which(is.na(all_symptomatic$WHORevisedFinal)),"WHORevisedFinal"] = "Caso sin datos de alarma" #if NA, assume lowest category (should confirm this)
all_symptomatic[which(all_symptomatic$WHORevisedFinal==""),"WHORevisedFinal"] = "Caso sin datos de alarma" #if empty, assume lowest category
clinical_cohort_subset = merge(cohort_subset, all_symptomatic[c("code","WHORevisedFinal","CareLevelFinal")], by="code", all.x=T)

#compare Dx from data subset with full data
#reformat the sample code variable so that the letter is lowercase to be compatable with other data
#need to take the characters up till the first dot, and then reformat
formatID = function(x) sprintf("ID%04d",as.numeric(strsplit(as.character(x["Sample"]), "\\.")[[1]][1:1]))
newIDs_3 = apply(X=basic_cohortD, MARGIN = 1, FUN=formatID)
basic_cohortD["code"] = newIDs_3
compare = merge(clinical_cohort_subset, clinical_cohortD_clean, by.x="Sample", by.y="Cod_Nin") #
