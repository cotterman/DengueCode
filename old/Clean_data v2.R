
################################################################################
############### Prepare and describe data for Dengue Dx project ################
################################################################################

library(lattice)
library(SuperLearner)
library(gtools) #enables smartbind

#### Please see "Guide to Dengue Data.xls" for documentation on data used in this file

#rm(list = ls()) #start with blank slate (clear everything from workspace)

#### Establish directories ####

inputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_inputs/"
outputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_outputs/"
resultsDir = "/home/carolyn/dengue_dx/Dengue_data/R_results/"


####################################################################
#################### Clean the abundance data ######################
####################################################################


####### Data from the first batch of samples ########


##### use data on which Natalia imposed 10 percent filter

abundancesD1 = read.delim(paste(inputsDir,"MS abundances_10percent_first batch.txt", sep=""), header=TRUE, nrows=100)
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
take2decimals = function(x) round(x, digits=2)
MZ_Nums_short = sapply(X=MZ_Nums, take2decimals)
MZ_Names = paste("MZ_", sapply(MZ_Nums, as.character), sep="") #no rounding
#MZ_Names = paste("MZ_", sapply(MZ_Nums_short, as.character), sep="") #with rounding
#MZ_Names #list of MZ names (with MZ prefix)

#drop the RT and MZ columns and relabel the Resp columns with MZ names
RespCols = seq(from=4, to=nColumns, by=3)
respD1_filter10n = abundancesD1[,c(1,RespCols)]
newnames = as.character(c("code", MZ_Names))
#there are two compounds named "MZ_336.32513".  I will relabel them here:
#newnames[newnames %in% "MZ_336.32513"] = c("MZ_336.32513a", "MZ_336.32513b")
colnames(respD1_filter10n) = newnames #rename columns

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
formatID = function(x, start, stop) sprintf("ID%04d",as.numeric(substr(as.character(x["code"]), start=start, stop=stop)))
newIDs_1 = apply(X=respD1_filter10n, MARGIN = 1, FUN=formatID, start=13, stop=17)
respD1_filter10n["code"] = newIDs_1

#add an indicator for type of study (will need this later when we merge)
respD1_filter10n["Study"] = "Hospital"

#add an indicator of LCMS_run (to keep track once merged with other round of LC-MS)
respD1_filter10n["LCMS_run"] = 1

#count  the number of duplicated MZ values
table(table(MZ_Nums_reshaped1)) #without rounding
table(table(MZ_Nums_short)) #with rounding


##### use data on which Natalia imposed 50 percent filter

#this script will run same commands as above, but with filter50n
source("/home/carolyn/dengue_dx/Dengue_code/clean_respD1_filter50n.R")



######## additional data manipulations ########

#merge with clinical info so that we have diagnoses
round1 = merge(clinical_comboD_clean, respD1_filter10n, by=c("code","Study"), all=FALSE)

#impose "50% filter" (consider 3 groups: DF, DHF/DSS, ND even though 1 of the analyses does not use ND observations)
   #(in future, can consider changing filter for NF vs. DHF/DSS analysis)
missingsDFpct = sapply(round1[which(round1$DxFinal3cat=="DF"), ], function(x) 
  sum(is.na(x)) / ( sum(is.na(x)) + sum(!is.na(x)) ))
missingsDHFDSSpct = sapply(round1[which(round1$DxFinal3cat=="DHF_DSS"), ], function(x) 
  sum(is.na(x)) / ( sum(is.na(x)) + sum(!is.na(x)) ))
missingsNDpct = sapply(round1[which(round1$DxFinal3cat=="ND"), ], function(x) 
  sum(is.na(x)) / ( sum(is.na(x)) + sum(!is.na(x)) ))
missingsOverallpct = sapply(round1, function(x) 
  sum(is.na(x)) / ( sum(is.na(x)) + sum(!is.na(x)) ))
#drop variables that do not have at least 50% presence in at least one group
round1_filter50a = round1[, which(missingsDFpct<.50 | missingsDHFDSSpct<.50 | missingsNDpct<.50)]
#drop variables that do not have at least 50% presence in data (regardless of Dx)
round1_filter50b = round1[, which(missingsOverallpct<.50)]

#missingness
missingsD = sapply(round1, function(x) sum(is.na(x)))
test = missingsD[which(missingsD==0)] #examine the compounds that are never missing
names(test) #1159 variables contain no missings
table(round1$PCR)
missingsD["PCR"]
table(round1$IR)
missingsD["IR"]
histogram(missingsD, xlim=c(0,30), xlab="number of missing values within each variable")
summary(missingsD)
quantile(missingsD, c(30, 35, 37, 40, 90, 95, 99)/100) #more percentiles

#do we have any zeros for abundance?
zerosD = sapply(round1, function(x) sum(x==0, na.rm=TRUE))
sum(zerosD!=0, na.rm=TRUE) #no zeros for abundance data

#verify that compounds in Natalia's 10% filtered data are subset of compounds in Natalia's 50% filtered data
commonvars = intersect(colnames(respD1_filter10n),colnames(respD1_filter50n)) #747 in common (good)

#compare compounds after imposing my filter to compounds after Natalia's 50% filter
commonvars = intersect(colnames(round1_filter50a),colnames(respD1_filter50n)) #742 in common (close)
#compounds that Natalia dropped but I did not
vars_in_mine_only = setdiff(colnames(round1_filter50a),colnames(respD1_filter50n)) #3550
vars_in_nats_only = setdiff(colnames(respD1_filter50n),colnames(round1_filter50a)) #5
vars_that_should_be_in_nats = intersect(names(test),vars_in_mine_only)
vars_that_should_be_in_nats

#consider missings to be have zero abundance, add 1, then transform to the log base 2 scale



#### simple analyses ####

#obtain features that demonstrate at least a 2-fold differentiation based on t-test with a p-value<.05.  
  #Take the Xx such features with the lowest p-value.

#Final model is a logistic model using these Xx features.



####### Data from the 3rd and 4th batch of samples ########

abundancesD2 = read.delim(paste(inputsDir,"MS abundances_10percent_batches 3 and 4.txt", sep=""), header=TRUE, nrows=200)
nColumns = dim(abundancesD2)[2] #number of columns in abundance data
#MZ = mass to charge of the compound, it usually is approximately +1 of the mass.  
#This value distinguishes the type of metabolite.  
#RT = retention times.  Can be ignored.
#Resp = response (aka abundance).

#obtained list of compounds (in same order in which they appear in spreadsheet columns)
# take first row of values for variables MZ, MZ.1, MZ.2.,,,,MZ.2232
MZ_Cols = seq(from=7, to=nColumns-2, by=3)
#MZ_Cols #these are the column numbers that we will want to keep
MZ_Nums = abundancesD2[1,MZ_Cols] #take the first row of data containing MZ values
#MZ_Nums #these are the values in the specified columns
MZ_Nums_reshaped2 = t(MZ_Nums) 
take2decimals = function(x) round(x, digits=2)
MZ_Nums_short = sapply(X=MZ_Nums, take2decimals)
MZ_Names = paste("MZ_", sapply(MZ_Nums_short, as.character), sep="")
#MZ_Names #list of MZ names (with MZ prefix)

#drop the RT and MZ columns and relabel the Resp colunmns with MZ names
RespCols = seq(from=9, to=nColumns, by=3)
respD2_filter10n = abundancesD2[,c(1,2,4,5,RespCols)] #select columns to keep
newnames = as.character(c("Study","code","Res_Final2","OMS", MZ_Names))
#there are two compounds named "MZ_336.32513".  I will relabel them here:
#newnames[newnames %in% "MZ_336.32513"] = c("MZ_336.32513a", "MZ_336.32513b")
colnames(respD2_filter10n) = newnames #rename columns

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode =    function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
newIDs_1 = apply(X=respD2_filter10n, MARGIN = 1, FUN=fcode)
respD2_filter10n["code"] = newIDs_1

#add an indicator of LCMS_run (to keep track once merged with other round of LC-MS)
respD2_filter10n["LCMS_run"] = 2

#count  the number of duplicated MZ values
table(table(MZ_Nums_reshaped2)) #without rounding
table(table(MZ_Nums_short)) #with rounding


###### View MZ values in graph ###### 

MZ_Nums_reshaped1 = t(MZ_Nums) #min value is 107, max is 1427
mybreaks = seq(min(MZ_Nums_reshaped1),max(MZ_Nums_reshaped1),10)
histogram(MZ_Nums_reshaped1, nint=100) #overview
histogram(MZ_Nums_reshaped1, xlim=c(300,400), nint=10000) #zoomed
histogram(MZ_Nums_reshaped1, xlim=c(330,380), nint=10000) #zoomed

### kernel density, overview
png("/home/carolyn/dengue_dx/R_results/MZ_density_full.png")
plot(density(MZ_Nums_reshaped1, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
     col="blue", xlab="", main="Kernel Density of MZ Values (full range)")
par(new=TRUE)
plot(density(MZ_Nums_reshaped2, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
     col="red", xlab="", main="") 
legend("topright",legend=c("1st LC-MS run","2nd LC-MS run"), col=c("blue","red"),lty=1)
dev.off()

### kernel density, zoomed in

png("/home/carolyn/dengue_dx/R_results/MZ_density_100.png")
plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
     col="blue", main="Kernel Density of MZ Values", xlab="") #kernal density
lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
      col="red", main="", xlab="") #kernal density
legend("topleft",legend=c("1st LC-MS run","2nd LC-MS run"), col=c("blue","red"),lty=1)
dev.off()

plot_kdensity = function(mz_start, mz_end, ylimit) {
  png(paste("/home/carolyn/dengue_dx/R_results/MZ_density_",mz_start,".png", sep=""))
  plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
       col="blue", main="", xlab="") #kernal density
  lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
        col="red", main="", xlab="") #kernal density
  dev.off()
}
plot_kdensity(mz_start=200, mz_end=300, ylimit=3000)
plot_kdensity(mz_start=300, mz_end=400, ylimit=3000)
plot_kdensity(mz_start=400, mz_end=500, ylimit=3000)
plot_kdensity(mz_start=500, mz_end=600, ylimit=3000)
plot_kdensity(mz_start=600, mz_end=700, ylimit=3000)
plot_kdensity(mz_start=700, mz_end=800, ylimit=3000)
plot_kdensity(mz_start=800, mz_end=900, ylimit=3000)
plot_kdensity(mz_start=900, mz_end=1000, ylimit=3000)
plot_kdensity(mz_start=1000, mz_end=1100, ylimit=1000)
plot_kdensity(mz_start=1100, mz_end=1200, ylimit=1000)
plot_kdensity(mz_start=1200, mz_end=1300, ylimit=1000)
plot_kdensity(mz_start=1300, mz_end=1400, ylimit=1000)
plot_kdensity(mz_start=1400, mz_end=1500, ylimit=1000)


####### Combine the abundance data ########

#list the columns that are in common 
commonvars = intersect(colnames(respD2_filter10n),colnames(respD1_filter10n)) #only 40 compound identifiers in common
resp_comboD = rbind(respD1_filter10n[,commonvars], respD2_filter10n[,commonvars])

#just get list of ID codes (so I can limit clinical data to patients for whom we have LC-MS data)
IDs_in_respData = resp_comboD[,c("code","Study","LCMS_run")]
#IDs_in_respData

##fails
#keepervars1 = c("code","Study","MZ_128.14265")
#keepervars2 = c("code","Study","MZ_899.49805")
#resp_comboD = bindData(respD1_filter10n[,keepervars1], respD2_filter10n[,keepervars2], common=commonvars)
#resp_comboD = smartbind(respD1_filter10n, respD2_filter10n, fill=NA) #fill with NA if not in both datasets
#colnames(resp_comboD)[1:300] #why did smartbind eliminate variables that are common?
#number of variables in common between respD1_filter10n and respD2_filter10n = 85 (by I find 40 above?)
#dim(respD1_filter10n)[2] + dim(respD2_filter10n)[2] - dim(resp_comboD)[2]



####################################################################
##################### Clean the clinical data ######################
####################################################################

####### List of clinical/lab variables to consider for analysis ########

clinic_varsD = read.delim(paste(inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)

#list names of variables from hospital data to include in either analysis
clinic_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 | 
                                    clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),"Variable.Name.in.Hospital.data"]
clinic_varsH_char = sapply(clinic_varsH, as.character) #convert factor to character

#corresponding names in cohort data
clinic_varsC = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),"Variable.Name.in.Cohort.data"]
clinic_varsC_char = sapply(clinic_varsC, as.character)


####### Data from the hospital study ########
# This should contain data from the first batch and from some of the 3rd/4th batches

clinical_hospitD = read.delim(paste(inputsDir,"Clinical data_hospital study.txt", sep=""), header=TRUE, nrows=2000)

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode = function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
newIDs_4 = apply(X=clinical_hospitD, MARGIN = 1, FUN=fcode)
clinical_hospitD["code"] = newIDs_4

#code the yes/no variables to take numeric values 1/0 (to be consistent with cohort data)
   #in future, could consider instead converting cohort data to factors (do algorithms care?)
YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                           "Variable.Name.in.Hospital.data"]
YesNo_varsH_char = sapply(YesNo_varsH, as.character)

for (name in YesNo_varsH_char) {
  print(name)
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


####### Data for cohort samples ########

basic_cohortD = read.delim(paste(inputsDir,"Basic info_cohort study.txt", sep=""), header=TRUE, nrows=100)

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
   #need to take the characters up till the first dot, and then reformat
formatID = function(x) sprintf("ID%04d",as.numeric(strsplit(as.character(x["Sample"]), "\\.")[[1]][1:1]))
newIDs_3 = apply(X=basic_cohortD, MARGIN = 1, FUN=formatID)
basic_cohortD["code"] = newIDs_3

#add an indicator for type of study (will need this later when we do the merge)
basic_cohortD["Study"] = "Cohort"

## add clinical/lab info
clinic_cohortD = read.delim(paste(inputsDir,"Clinical data_cohort study.txt", sep=""), header=TRUE, nrows=100)
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
#View(clinic_cohortD[,c("FTM","Fecha_Nac","FIS","age","DaysSick")])

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


######### Merge hospital study data with cohort data #######

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

#convert all yes/no questions into factors (then R will summarize variables as I want them)
  #TODO: figure out why this does not work
#for (name in YesNo_varsH_char) {
#  clinical_comboD[name] = factor(x = as.character(clinical_comboD[name]))
#}                                
                                 
#keep only observations for which we have LC-MS data
clinical_comboD = merge(IDs_in_respData, clinical_comboD, by=c("code","Study"), all=FALSE)
#drop observations with unknown final dengue dx
clinical_comboD_clean = clinical_comboD[which(!is.na(clinical_comboD$DxFinal4cat)),]

## Write this clinical data to file for easy future access ##
write.csv(x=clinical_comboD_clean, file=paste(outputsDir,"clinical_comboD_clean.txt"), row.names = FALSE)
#clinical_comboD_clean = read.csv(paste(outputsDir,"clinical_comboD_clean.txt"), header=TRUE)


### Summarize lab and clinical data ###

summary(clinical_comboD_clean)

#days between symptom onset and sample collection
table(clinical_comboD_clean$DaysSick) 

#diagnosis variables
table(clinical_comboD_clean$Res_Final) #no missings
table(clinical_comboD_clean$ClasificacionFinal) #no missings
table(clinical_comboD_clean$DxFinal4cat) #no missings

#2-way table
table(clinical_comboD_clean$ClasificacionPrimerDia, clinical_comboD_clean$DxFinal4cat) 
table(clinical_comboD_clean$ClasificacionFinal, clinical_comboD_clean$DxFinal4cat) 

#restrict to just LC-MS round 1
table(clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==1),"ClasificacionPrimerDia"], 
      clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==1),"DxFinal4cat"]) 
table(clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==1),"DaysSick"]) 

#restrict to just LC-MS round 2
table(clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==2),"ClasificacionPrimerDia"], 
      clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==2),"DxFinal4cat"]) 
table(clinical_comboD_clean[which(clinical_comboD_clean$LCMS_run==2),"DaysSick"]) 


##############################################################################
################ Combine abundance data with other data ######################
##############################################################################

comboD <- merge(resp_comboD, clinical_comboD_clean, by=c("code","Study"), all=FALSE)
dim(comboD)[1] #number of rows in resulting data -- a perfect match
#class(comboD$MZ_491.2914[[1]])


##############################################################################
######### Some basic data checks and descriptives of merged data #############
##############################################################################

#diagnosis variables
table(comboD$Res_Final, exclude=NULL) 

#drop observations with unknown final dengue dx
clean_comboD = comboD[which(!is.na(comboD$DxFinal4cat)),]

table(comboD$ClasificacionFinal) 
table(comboD$DxFinal4cat) 

#2-way tables
table(comboD$ClasificacionPrimerDia, comboD$DxFinal4cat) 
table(final_comboD$Study)
table(comboD[which(comboD$Study=="Hospital"),"ClasificacionPrimerDia"], 
      comboD[which(comboD$Study=="Hospital"),"DxFinal4cat"]) 

#days between symptom onset and sample collection
table(clean_comboD$DaysSick) 

#histogram of abundance values
#histogram(comboD$MZ_522.3552)
#histogram(c(comboD$MZ_522.3552 , comboD$MZ_287.63376))

#missingness
missingsD = sapply(clean_comboD, function(x) sum(is.na(x)))
missingsD[which(missingsD==0)] #examine the compounds that are never missing
histogram(missingsD, xlim=c(0,30), xlab="number of missing values within each observation")
summary(missingsD)
quantile(missingsD, c(30, 35, 37, 40, 90, 95, 99)/100) #more percentiles

#do we have any zeros for abundance?
zerosD = sapply(clean_comboD, function(x) sum(x==0, na.rm=TRUE))
sum(zerosD!=0, na.rm=TRUE) #no zeros

#create a version of the data which replaces values for missing compounds with zeros
  #other options: fill with min value observed for that compound,
      #fill with 1/2 of the min value, input value
#Note: I do not have any missings for non-MZ variables so this works, but beware in future
wzeros_comboD = clean_comboD
wzeros_comboD[is.na(clean_comboD)] <- 0


###############################################################################
####### Predict diagnosis using abundance data with Super Learner #############
###############################################################################


############################################################### 
#### dengue (DF/DHF/DSS) vs. non-dengue (ND) prediction ####### 
############################################################### 

#create dummy for non-dengue (0) versus dengue (1)
nrows_DF = dim(wzeros_comboD)[1] #number of rows
dummyDF = vector(mode="numeric", length=nrows_DF)
dummyDF[which(wzeros_comboD["DxFinal4cat"]=="ND")] = 0
dummyDF[which(wzeros_comboD["DxFinal4cat"]!="ND")] = 1
wzeros_comboD$dummyDF = dummyDF
class(wzeros_comboD$dummyDF) 

#obtain list of explanatory variables for super learner
myvars = colnames(wzeros_comboD)
explanatory_vars = myvars[grepl("MZ_", myvars)]
#explanatory_vars

#get the V-fold cross-validated risk estimate for super learner.
set.seed(999)
#SL.gam, SL.glm produced errors (singularities and non-convergence)
mySL.library = c("SL.mean", "SL.knn", "SL.randomForest")
SLresultsDF = CV.SuperLearner(Y=dummyDF, X=wzeros_comboD[explanatory_vars], 
                               V=10, family=binomial(), SL.library=mySL.library,
                               method = "method.NNLS", verbose=FALSE)
SLresultsDF
summary(SLresultsDF)
coef(SLresultsDF)


#Plot the CV risk to sample size to see if it’s starting to plateau.


############################################################
###### severe dengue (DSS) vs. non-severe dengue (DHF) #####
############################################################

#drop the ND and DF folks
DSS_and_DHF = comboD[which(comboD["diagnosis"]=="DSS" | comboD["diagnosis"]=="DHF"),]
nrows_DSS = dim(DSS_and_DHF)[1] #number of rows
dummyDSS = vector(mode="numeric", length=nrows_DSS)
dummyDSS[which(DSS_and_DHF["diagnosis"]=="DHS")] = 0
dummyDSS[which(DSS_and_DHF["diagnosis"]=="DSS")] = 1
DSS_and_DHF$dummyDSS = dummyDSS
DSS_and_DHF[, c("diagnosis","dummyDSS")]
class(DSS_and_DHF$dummyDSS) 

#drop variables with one or more missing values (SL cannot deal with missings)
#better to input missing values? (when/why do missings occur with LC-MS?)
dim(DSS_and_DHF)[2]  #number of columns before dropping those with missing values
missingsDSSDHF = sapply(DSS_and_DHF, function(x) sum(is.na(x)))
DSS_and_DHF = DSS_and_DHF[,which(missingsDSSDHF==0)]
dim(DSS_and_DHF)[2]  #number of columns after dropping those with missing values
#covariates drop from 744 to 448

#obtain list of explanatory variables for super learner
myvars = colnames(DSS_and_DHF)
explanatory_vars = myvars[grepl("MZ_", myvars)]
#explanatory_vars

#get the V-fold cross-validated risk estimate for super learner.
set.seed(999)
mySL.library = c("SL.glm", "SL.knn", "SL.randomForest")
SLresultsDSS = CV.SuperLearner(Y=dummyDSS, X=DSS_and_DHF[explanatory_vars], 
                               V=10, family=binomial(), SL.library=mySL.library,
                               method = "method.NNLS", verbose=FALSE)
SLresultsDSS
summary(SLresultsDSS)
coef(SLresultsDSS)

