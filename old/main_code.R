
################################################################################
############### Prepare and describe data for Dengue Dx project ################
################################################################################

rm(list = ls()) #start with blank slate (clear everything from workspace)

library(lattice)
library(SuperLearner)
library(gtools) #enables smartbind
#library(leaps) #for best subset selection
library(bestglm) #for best subset selection
library(MESS) #for area under the curve
library(randomForest)
library(glmnet) #lasso
#library(xlsReadWrite) #need newer version of R for this package

#### Please see "Guide to Dengue Data.xls" for documentation on data used in this file


#### Establish directories ####

inputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_inputs/"
outputsDir = "/home/carolyn/dengue_dx/Dengue_data/Processed/R_outputs/"
resultsDir = "/home/carolyn/dengue_dx/R_results/"


#### Esablish functions to be run #### 
source("/home/carolyn/dengue_dx/Dengue_code/clean_data_functions.R")


###############################################################################
######################### Clean the abundance data ############################
###############################################################################


#notes on function call:
#if roundme parameter is not specified, then it is assumed to be false and resulting dataset will contain original MZ values

#for Nicaragua Serum samples
respD1_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_first batch.txt", lcms_run=1, printme=TRUE)
respD1_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_first batch.txt", lcms_run=1, printme=TRUE)
respD2_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #only 80 samples
respD2_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples

#for non-invasive Nicaraguan samples (urine and saliva)
respD3_filter50n = clean_LCMS(infile="LCMS_saliva_Nica_50percent.txt", lcms_run=3, printme=TRUE) #86 (as in LCMS excel)
temp = clean_LCMS(infile="LCMS_urine_Nica_50percent.txt", lcms_run=5, printme=TRUE) #91 (as in LCMS excel)
#get rid of duplicated code/study ID (note: this wil keep one of the observation in each duplicate group -- will remove later)
respD5_filter50n = temp[!duplicated(temp[,c("code","Study")]),] #86 (and still more to remove)

####### Compare frequences of MZ values across datasets (histograms) ########

graph_MZ_frequencies(respD1_filter50n, respD2_filter50n, "Round1_filter50", "Round2_filter50", "f50_")
graph_MZ_frequencies(respD1_filter10n, respD2_filter10n, "Round1_filter10", "Round2_filter10", "f10_")


####### Combine the abundance data (runs 1 and 2) ########

#list the columns that are in common between serum runs 1 and 2
commonvars = intersect(colnames(respD2_filter10n),colnames(respD1_filter10n)) #only 40 compound identifiers in common (10% data)
resp_comboD = rbind(respD1_filter10n[,commonvars], respD2_filter10n[,commonvars])

#just get list of ID codes (so I can limit clinical data to patients for whom we have LC-MS data)
IDs_in_resp_D1_D2 = resp_comboD[,c("code","Study","LCMS_run")]
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
#are any saliva and serum samples from same people? just 1 (might be mistake -- ID1408 in cohort, which was already problematic)
table(IDs_in_resp_all$serum, IDs_in_resp_all$saliva, useNA="ifany")


######### Other checks and descriptives of abundance data #############

#missingness (exclude "LCMS_run","Study", and "code" columns)
missingsD = sapply(respD1_filter50n[grep("MZ_",colnames(respD1_filter50n))], function(x) sum(is.na(x))) 
#missingsD[which(missingsD==0)] #examine the compounds that are never missing
fileprefix = "R1_filter50"
png(paste("/home/carolyn/dengue_dx/R_results/",fileprefix,"_missingMZ.png", sep=""))
histogram(missingsD, xlim=c(0,30), xlab="number of observations with missing values within each MZ value",
          main="LC-MS Run 1 with 50% filter")
dev.off()
summary(missingsD)
table(missingsD)
sum(missingsD==0)/dim(respD1_filter50n)[2] #fraction of compounds that have no missings
sum(missingsD<2)/dim(respD1_filter50n)[2] #fraction of compounds that have 0 or 1 missings
quantile(missingsD, c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99)/100) #more percentiles
#do we have any zeros for abundance?
zerosD = sapply(resp_comboD, function(x) sum(x==0, na.rm=TRUE))
sum(zerosD!=0, na.rm=TRUE) #no zeros


####################################################################
##################### Clean the clinical data ######################
####################################################################

#Process text files
#Create clinical_comboD 
clinic_varsD = read.delim(paste(inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
source("/home/carolyn/dengue_dx/Dengue_code/clean_clinical_data.R") #produces clinical_comboD


#drop observations with unknown final dengue dx
clinical_comboD_prelim = clinical_comboD[which(!is.na(clinical_comboD$DxFinal4cat)),]
#keep only observations for which we have LC-MS data
clinical_comboD_clean = merge(IDs_in_resp_all, clinical_comboD_prelim, by=c("code","Study"), all=FALSE)

## Write this clinical data to file for easy future access ##
write.csv(x=clinical_comboD_clean, file=paste(outputsDir,"clinical_comboD_clean.txt"), row.names = FALSE)
#clinical_comboD_clean = read.csv(paste(outputsDir,"clinical_comboD_clean.txt"), header=TRUE)
#WARNING: after reading in the csv, must reconvert variables to factors


####### Summarize lab and clinical data ####### 

summarize_clinical(clinical_comboD_clean)


##############################################################################
################ additional investigation of abundance data ##################
##############################################################################

#source("/home/carolyn/dengue_dx/Dengue_code/investigate_abundance_filtering.R")
#CONCLUSION: the filtering that Natalia is doing differs from what we thought



##############################################################################
################ Combine abundance with clinical data ########################
##############################################################################

#Data that includes serum runs 1 and 2
comboD <- merge(clinical_comboD_clean, resp_comboD, by=c("code","Study"), all=F)
dim(comboD)[1] #number of rows in resulting data -- a perfect match

#Data for serum run 1 only
comboD1_filter50n <- merge(clinical_comboD_clean, respD1_filter50n, by=c("code","Study"), all=F) #88 obs
#Data for serum run 2 only
comboD2_filter50n <- merge(clinical_comboD_clean, respD2_filter50n, by=c("code","Study"), all=F) #75 obs

#Data for Nicaragua saliva
temp <- merge(clinical_comboD_clean, respD3_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n = temp[which(temp$code!="ID1408"),]  #85
table(comboD3_filter50n$DxFinal4cat) #cannot to DF vs DHF/DSS analysis with such few obs

#Data for Nicaragua urine
temp <- merge(clinical_comboD_clean, respD5_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID and those which were duplicated in this data
comboD5_filter50n = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                    temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
table(comboD5_filter50n$DxFinal4cat) #cannot to DF vs DHF/DSS analysis with such few obs


###############################################################################
#################### Variable Importance Analysis without CV ##################
###############################################################################

#### Esablish functions to be run #### 
source("/home/carolyn/dengue_dx/Dengue_code/prediction_functions.R")


D1_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, dim_reduce_covar=F, "D1_covarF")
D1_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, dim_reduce_covar=T, "D1_covarT")
D1_VIM = merge(D1_VIM_covarT, D1_VIM_covarF, by="row.names", all=T)

D2_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, dim_reduce_covar=F, "D2_covarF")
D2_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, dim_reduce_covar=T, "D2_covarT")
D2_VIM = merge(D2_VIM_covarT, D2_VIM_covarF, by="row.names", all=T)

D3_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, dim_reduce_covar=F, "D3_covarF")
D3_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, dim_reduce_covar=T, "D3_covarT")
D3_VIM = merge(D3_VIM_covarT, D3_VIM_covarF, by="row.names", all=T)

D5_VIM_covarF = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, dim_reduce_covar=F, "D5_covarF")
D5_VIM_covarT = get_VIM_RF(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, dim_reduce_covar=T, "D5_covarT")
D5_VIM = merge(D5_VIM_covarT, D5_VIM_covarF, by="row.names", all=T)

## combine all variable importance lists ##
serum_VIM = merge(D1_VIM, D2_VIM, by="Row.names", all=T)
noninvasive_VIM = merge(D3_VIM, D5_VIM, by="Row.names", all=T)
VIM_RF = merge(serum_VIM, noninvasive_VIM, by="Row.names", all=T)
write.csv(x=VIM_RF, file=paste(resultsDir,"VIM_RF_DEN.txt"), row.names = FALSE)


###############################################################################
############################ Predict diagnosis ################################
###############################################################################


#### Specify location to place output ####
sink(paste(resultsDir,"prediction_output_DEN_D3.txt",sep=""), append=FALSE, split=TRUE)


#### Esablish functions to be run #### 
source("/home/carolyn/dengue_dx/Dengue_code/prediction_functions.R")


### Non-invasive samples from Nicaragua -- ND vs DEN only ###

resultsD3_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD3_filter50n, "saliva, D3", reduce_in_CV=T)
#selected_D1_MFs_DEN = resultsD1_DEN[[2]]
resultsD3_DEN[[1]]
write.csv(x=resultsD3_DEN[[1]], file=paste(resultsDir,"resultsD3_DEN_v7.txt"), row.names = FALSE)

resultsD5_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD5_filter50n, "urine, D5", reduce_in_CV=T)
#selected_D1_MFs_DEN = resultsD1_DEN[[2]]
resultsD5_DEN[[1]]
write.csv(x=resultsD5_DEN[[1]], file=paste(resultsDir,"resultsD5_DEN_v7.txt"), row.names = FALSE)


### Serum Samples from Nicaragua ###

resultsD1_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD1_filter50n, "serum, D1", reduce_in_CV=T)
#selected_D1_MFs_DEN = resultsD1_DEN[[2]]
resultsD1_DEN[[1]]
write.csv(x=resultsD1_DEN[[1]], file=paste(resultsDir,"resultsD1_DEN_v6a.txt"), row.names = FALSE)

resultsD1_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD1_filter50n, "serum, D1", reduce_in_CV=T)
#selected_D1_MFs_SDEN = resultsD1_SDEN[[2]]
write.csv(x=resultsD1_SDEN[[1]], file=paste(resultsDir,"resultsD1_SDEN_v6a.txt"), row.names = FALSE)

resultsD2_DEN = run_predictions(clinic_varsD, "ND.vs.DEN", comboD2_filter50n, "serum, D2", reduce_in_CV=T)
resultsD2_DEN[[1]]
selected_D2_MFs_DEN = resultsD2_DEN[[2]]
write.csv(x=resultsD2_DEN[[1]], file=paste(resultsDir,"resultsD2_DEN_v6a.txt"), row.names = FALSE)

resultsD2_SDEN = run_predictions(clinic_varsD, "DF.vs.DHF.DSS", comboD2_filter50n, "serum, D2", reduce_in_CV=T)
selected_D2_MFs_SDEN = resultsD2_SDEN[[2]]
write.csv(x=resultsD2_SDEN[[1]], file=paste(resultsDir,"resultsD2_SDEN_v6a.txt"), row.names = FALSE)

## combine lists of MFs selected by dimension reduction methods
flist_DEN = merge(selected_D1_MFs_DEN, selected_D2_MFs_DEN,  by="Var1", all=TRUE)
flist_DEN
write.csv(x=flist_DEN, file=paste(resultsDir,"MFs_selected_by_RF_DEN_v6a.txt"), row.names = FALSE)
flist_SDEN = merge(selected_D1_MFs_SDEN, selected_D2_MFs_SDEN,  by="Var1", all=TRUE)
flist_SDEN
write.csv(x=flist_SDEN, file=paste(resultsDir,"MFs_selected_by_RF_SDEN_v6a.txt"), row.names = FALSE)
