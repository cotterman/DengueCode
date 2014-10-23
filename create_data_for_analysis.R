###############################################################################
####################### Create data to be used in analysis ####################
###############################################################################
# Refer to main_code.R to establish directories
# See "Guide to Dengue Data.xls" for documentation on data used in this file


#### Esablish functions to be run #### 
source(paste(codeDir, "clean_data_functions.R", sep=""))


###############################################################################
#################### Clean the abundance data #################################
###############################################################################

#for Nicaragua Serum samples
respD1_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_first batch.txt", lcms_run=1, printme=TRUE) #88 samples
respD1_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_first batch.txt", lcms_run=1, printme=TRUE)
respD2_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples
respD2_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples

#for non-invasive Nicaraguan samples (urine and saliva)
respD3_filter50n = clean_LCMS(infile="LCMS_saliva_Nica_50percent.txt", lcms_run=3, printme=TRUE) #86 (as in LCMS excel)
temp = clean_LCMS(infile="LCMS_urine_Nica_50percent.txt", lcms_run=5, printme=TRUE) #91 (as in LCMS excel)
#get rid of duplicated code/study ID (note: this wil keep one of the observation in each duplicate group -- will remove later)
respD5_filter50n = temp[!duplicated(temp[,c("code","Study")]),] #86 (and still more to remove)

# Combine the abundance data (runs 1 and 2) 
#list the columns (mz values) that are in common between serum runs 1 and 2
commonvars = intersect(colnames(respD2_filter50n),colnames(respD1_filter50n)) #only 40 compound identifiers in common (10% data)
resp_D1D2 = rbind(respD1_filter50n[,commonvars], respD2_filter50n[,commonvars])


###############################################################################
########################### Clean the clinical data ###########################
###############################################################################

source(paste(codeDir, "clean_clinical_data.R", sep="")) #contains clean_clinical_data function
clinic_varsD = read.delim(paste(clinical_inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
clinical_full_clean = clean_clinical_data(clinic_varsD) #1726 observations (no restriction based on LCMS presence)
save(clinical_full_clean, file=paste(outputsDir,"clinical_full_clean.RData", sep=""))

## restrict clinical data sample to observations in LCMS data that have known diagnosis ##
#drop observations with unknown final dengue dx
clinical_prelim = clinical_full_clean[which(!is.na(clinical_full_clean$DxFinal4cat)),]
#Get list of patient ID codes that are present in LCMS data (Nicaragua serum R1 and R2, urine, and saliva)
IDs_in_resp_all = get_IDs_in_common(resp_D1D2, respD3_filter50n, respD5_filter50n)
#keep only clinical observations for which we have LC-MS data
clinical_restricted_clean = merge(IDs_in_resp_all, clinical_prelim, by=c("code","Study"), all=FALSE)

## Write this clinical data to file for easy future access ##
save(clinical_restricted_clean, file=paste(outputsDir,"clinical_restricted_clean.RData", sep=""))
#now load file in future programs like this (will load dataframe object named "clinical_restricted_clean.RData")
#load(paste(outputsDir,"clinical_restricted_clean.RData", sep="")) 

## Summarize lab and clinical data ##
summarize_clinical(clinical_restricted_clean)


###############################################################################
################ Combine abundance with clinical data #########################
###############################################################################

#Data that includes serum runs 1 and 2
comboD1D2 <- merge(clinical_restricted_clean, resp_D1D2, by=c("code","Study"), all=F)
dim(comboD1D2)[1] #number of rows in resulting data -- a perfect match

#Data for serum run 1 only
comboD1_filter50n <- merge(clinical_restricted_clean, respD1_filter50n, by=c("code","Study"), all=F) #88 obs
save(comboD1_filter50n, file=paste(outputsDir,"comboD1_filter50n.RData", sep=""))
#Data for serum run 2 only
comboD2_filter50n <- merge(clinical_restricted_clean, respD2_filter50n, by=c("code","Study"), all=F) #75 obs

#Data for Nicaragua saliva
temp <- merge(clinical_restricted_clean, respD3_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n = temp[which(temp$code!="ID1408"),]  #85
table(comboD3_filter50n$DxFinal4cat) #cannot to DF vs DHF/DSS analysis with only 7 DHF/DSS (!)

#Data for Nicaragua urine
temp <- merge(clinical_restricted_clean, respD5_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID and those which were duplicated in this data
comboD5_filter50n = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                                 temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
table(comboD5_filter50n$DxFinal4cat) #cannot to DF vs DHF/DSS analysis with only 7 DHF/DSS (!)

