################################################################################
############### Prepare and describe data for Dengue Dx project ################
################################################################################

library(lattice)
library(SuperLearner)

#### Please see "Notes on data.txt" for documentation on data used in this file

#### Establish directories ####

processedDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/"


####################################################################
##################### Clean the clinical data ######################
####################################################################

####### Data from the hosptial study ########
# I think this should contain data from first batch and also from some of the 3rd/4th batches

clinical_hospitD = read.delim(paste(processedDir,"Clinical data_hospital study.txt", sep=""), header=TRUE, nrows=2000)

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode = function(x) sprintf("ID%04d",as.numeric(as.character(x["Code"])))
newIDs_4 = apply(X=clinical_hospitD, MARGIN = 1, FUN=fcode)
#newIDs_4
clinical_hospitD["code"] = newIDs_4

#days between symptom onset and sample collection
table(clinical_hospitD$Day.of.illness) 

#values taken by the diagnosis variable
table(clinical_hospitD$Res_Final) #no missings
table(clinical_hospitD$ClasificacionFinal) #no missings

#2-way table
table(clinical_hospitD$ClasificacionPrimerDia, clinical_hospitD$ClasificacionFinal) 
table(clinical_hospitD$Res_Final, clinical_hospitD$ClasificacionFinal)


####### Data from the first batch of samples ########

clinicalD1 = read.delim(paste(processedDir,"Clinical data_batch 1.txt", sep=""), header=TRUE, nrows=100)

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode = function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
newIDs_3 = apply(X=clinicalD1, MARGIN = 1, FUN=fcode)
newIDs_3
clinicalD1["code"] = newIDs_3

#days between symptom onset and sample collection
table(clinicalD1$DiaEnfermedadBD) 

#values taken by the diagnosis variable
table(clinicalD1$Res_Final) #no missings
table(clinicalD1$ClasificacionFinal) #no missings

#2-way table
table(clinicalD1$ClasificacionPrimerDia, clinicalD1$ClasificacionFinal) 
table(clinicalD1$Res_Final, clinicalD1$ClasificacionFinal)

#Is first batch completely contained within hospital study data?  YES
testc = merge(clinicalD1, clinical_hospitD, by="code", all=TRUE)


####################################################################
#################### Clean the abundance data ######################
####################################################################

####### Data from the 3rd and 4th batch of samples ########

hospitD = read.delim(paste(processedDir,"codes_in_hospital study.txt", sep=""), header=TRUE, nrows=100)
fcode = function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
newIDs_5 = apply(X=hospitD, MARGIN = 1, FUN=fcode)
hospitD["code"] = newIDs_5
combo_hospitD <- merge(hospitD, clinical_hospitD, by="code", all=FALSE)


#values taken by the diagnosis variable
table(combo_hospitD$Res_Final) #no missings
table(combo_hospitD$ClasificacionFinal) #no missings

#2-way table
table(combo_hospitD$ClasificacionPrimerDia, combo_hospitD$ClasificacionFinal) 
table(combo_hospitD$Res_Final, combo_hospitD$ClasificacionFinal)


####### Data from the fist batch of samples ########

abundancesD = read.delim(paste(processedDir,"MS abundances_10percent_first batch.txt", sep=""), header=TRUE, nrows=100)
nColumns = dim(abundancesD)[2] #number of columns in abundance data
#MZ = mass to charge of the compound, it usually is approximately +1 of the mass.  
  #This value distinguishes the type of metabolite.  
#RT = retention times.  Can be ignored.
#Resp = response (aka abundance).

#obtained list of compounds (in same order in which they appear in spreadsheet columns)
  # take first row of values for variables MZ, MZ.1, MZ.2.,,,,MZ.2232
MZ_Cols = seq(from=2, to=nColumns-2, by=3)
#MZ_Cols #these are the column numbers that we will want to keep
MZ_Nums = abundancesD[1,MZ_Cols] #take the first row of data containing MZ values
#MZ_Nums #these are the values in the specified columns
MZ_Names = paste("MZ_", sapply(MZ_Nums, as.character), sep="")
#MZ_Names #list of MZ names (with MZ prefix)

#drop the RT and MZ columns and relabel the Resp colunmns with MZ names
RespCols = seq(from=4, to=nColumns, by=3)
respD = abundancesD[,c(1,RespCols)]
newnames = as.character(c("code", MZ_Names))
#there are two compounds named "MZ_336.32513".  I will relabel them here:
#newnames[newnames %in% "MZ_336.32513"] = c("MZ_336.32513a", "MZ_336.32513b")
colnames(respD) = newnames #rename columns

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
formatID = function(x, start, stop) sprintf("ID%04d",as.numeric(substr(as.character(x["code"]), start=start, stop=stop)))
newIDs_1 = apply(X=respD, MARGIN = 1, FUN=formatID, start=13, stop=17)
newIDs_1
respD["code"] = newIDs_1


###############################################################################
############## Clean the data containing classification and dates #############
###############################################################################

#This data needs to be updated (it is from excel file that contained mistakes)

class_and_datesD = read.delim(paste(processedDir,"classifications and dates_first batch.txt", sep=""), header=TRUE, nrows=100)

#list the variables in this data
colnames(class_and_datesD)

#check out the variable types -- make sure they are proper
class(class_and_datesD["symptom_days"][[1]])
sapply(class_and_datesD, class) #shortcut to get type of all variables

#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fID = function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
newIDs_2 = apply(X=class_and_datesD, MARGIN = 1, FUN=fID)
newIDs_2
class_and_datesD["code"] = newIDs_2

#keep only variables we care about
dx_and_daysD = class_and_datesD[c("code","diagnosis","symptom_days")]


##############################################################################
#### Combine abundance data with data containing classification and dates ####
##############################################################################

comboD <- merge(respD, clinicalD, by="code", all=TRUE)
dim(abundancesD)[1] #number of rows in resulting data -- a perfect match
class(comboD$MZ_491.2914[[1]])


##############################################################################
######### Some basic data checks and descriptives of merged data #############
##############################################################################

#values taken by the diagnosis variable
table(comboD$diagnosis) #no missings

#days between symptom onset and sample collection
table(comboD$symptom_days) #1 missing value

#2-way table
table(comboD$symptom_days, comboD$diagnosis) 

#histogram of abundance values
histogram(comboD$MZ_522.3552)
histogram(c(comboD$MZ_522.3552 , comboD$MZ_287.63376))

#missingness
missingsD = sapply(comboD, function(x) sum(is.na(x)))
histogram(missingsD, xlim=c(0,30), xlab="number of missing values within each observation")
summary(missingsD)
quantile(missingsD, c(30, 35, 37, 40, 90, 95, 99)/100) #more percentiles
#so we loose almost 40% of SMBs if they cannot be missing

#do we have any zeros for abundance?
zerosD = sapply(comboD, function(x) sum(x==0, na.rm=TRUE))
sum(zerosD!=0, na.rm=TRUE) #no zeros


###############################################################################
####### Predict diagnosis using abundance data with Super Learner #############
###############################################################################


########################################################
#### dengue (DF) vs. non-dengue (NEG) prediction ####### 
########################################################

#drop the DSS and DHS folks
DF_and_NEG = comboD[which(comboD["diagnosis"]=="NEG" | comboD["diagnosis"]=="DF"),]
nrows_DF = dim(DF_and_NEG)[1] #number of rows
dummyDF = vector(mode="numeric", length=nrows_DF)
dummyDF[which(DF_and_NEG["diagnosis"]=="NEG")] = 0
dummyDF[which(DF_and_NEG["diagnosis"]=="DF")] = 1
DF_and_NEG$dummyDF = dummyDF
DF_and_NEG[, c("diagnosis","dummyDF")]
class(DF_and_NEG$dummyDF) 

#drop variables with one or more missing values (SL cannot deal with missings)
  #better to input missing values? (when/why do missings occur with LC-MS?)
dim(DF_and_NEG)[2]  #number of columns before dropping those with missing values
missingsDFNEG = sapply(DF_and_NEG, function(x) sum(is.na(x)))
DF_and_NEG = DF_and_NEG[,which(missingsDFNEG==0)]
dim(DF_and_NEG)[2]  #number of columns after dropping those with missing values
#covariates drop from 744 to 367

#obtain list of explanatory variables for super learner
myvars = colnames(DF_and_NEG)
explanatory_vars = myvars[grepl("MZ_", myvars)]
#explanatory_vars

#get the V-fold cross-validated risk estimate for super learner.
set.seed(999)
mySL.library = c("SL.glmnet", "SL.gam", "SL.mean", "SL.glm", "SL.knn", "SL.randomForest")
SLresultsDF = CV.SuperLearner(Y=dummyDF, X=DF_and_NEG[explanatory_vars], 
                               V=10, family=binomial(), SL.library=mySL.library,
                               method = "method.NNLS", verbose=FALSE)
SLresultsDF
summary(SLresultsDF)
coef(SLresultsDF)


#########################################################
#### severe dengue (DSS) vs. non-severe dengue (DHS) ####
#########################################################

#drop the DF and NEG folks
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

