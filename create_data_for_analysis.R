###############################################################################
####################### Create data to be used in analysis ####################
###############################################################################
# Refer to main_code.R to establish directories and functions
# See "Guide to Dengue Data.xls" for documentation on data used in this file


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
#get rid of duplicated code/study ID (note: this will keep one of the observation in each duplicate group -- will remove later)
respD5_filter50n = temp[!duplicated(temp[,c("code","Study")]),] #86 (and still more to remove)

# Combine the abundance data (runs 1 and 2) 
#list the columns (mz values) that are in common between serum runs 1 and 2
commonvars = intersect(colnames(respD2_filter50n),colnames(respD1_filter50n)) #only 40 compound identifiers in common (10% data)
resp_D1D2 = rbind(respD1_filter50n[,commonvars], respD2_filter50n[,commonvars])

#Get list of patient ID codes that are present in LCMS data (Nicaragua serum R1 and R2, urine, and saliva)
IDs_in_resp_all = get_IDs_in_common(resp_D1D2, respD3_filter50n, respD5_filter50n)
save(IDs_in_resp_all, file=paste(outputsDir,"IDs_in_resp_all.RData", sep="")) #for future access


###############################################################################
########################### Clean the clinical data ###########################
###############################################################################

# this function was established in "clean_clinical_data_functions_v2.R"
clin24_full = clean_clin24_data(clinic_varsD, IDs_in_resp_all) #1726 observations (no restriction based on LCMS presence)

# drop observations with unknown final dengue dx
clin24_full_clean = clin24_full[which(!is.na(clin24_full$WHOFinal4cat)),]

#keep only clinical observations for which we have LC-MS data
load(paste(outputsDir,"IDs_in_resp_all.RData", sep="")) #loads IDs_in_resp_all
clin24_restricted_clean = merge(IDs_in_resp_all, clin24_full_clean, by=c("code","Study"), all=FALSE)

#keep only clinical observations for batch 1 LCMS data
clin24_D1_clean = merge(respD1_filter50n[,c("code","Study")], clin24_full_clean, by=c("code","Study"), all=FALSE)

## Write this clinical data to file for easy future access ##
save(clin24_full_clean, file=paste(outputsDir,"clin24_full_clean.RData", sep=""))
save(clin24_restricted_clean, file=paste(outputsDir,"clin24_restricted_clean.RData", sep=""))
save(clin24_D1_clean, file=paste(outputsDir,"clin24_D1_clean.RData", sep=""))

#now load file in future programs like this (will load dataframe object named "clin24_restricted_clean.RData")
#load(paste(outputsDir,"clin24_restricted_clean.RData", sep="")) 


###############################################################################
################ Additional (optional) data manipulations #####################
###############################################################################

#get list of clinical variables for purposes of imputing missing values
  #note: right now I'm using vars that are part of the DHF/DSS definition
    #letting outcome="either" here is ok for ND vs DEN since IR and PCR have missings
      #and I only use fully non-missing vars for imputations (but beware if change procedure)
varlist = get_clinic_var_list(clinic_varsD, outcome="either", 
                              eliminate_vars_with_missings=F, eliminate_constant_vars=T, 
                              eliminate_vars_with_minXnomiss=50,
                              XD=clin24_full_clean, restrict_to_cohort_vars=F, 
                              restrict_to_hospit_vars=T,UltraX=T, BloodLab=T)
#impute missing values and add indicator for imputation (todo: add SL as method option)
clin24_full_wImputedRF1 = impute_missings(XD=clin24_full_clean, vars_to_impute=varlist, 
                                        predictor_vars_prelim=varlist, exclude_miss_vars=T, method="RF")
#write to file
save(clin24_full_wImputedRF1, file=paste(outputsDir,"clin24_full_wImputedRF1.RData", sep=""))

#verify missings have been imputed
  #sum(is.na(clin_full_wImputedRF1[varlist])) #should be zero
  #compare values before/after
  #clin24_full_clean[c(1:25),"Eosi"]
  #clin_full_wImputedRF1[c(1:25),"Eosi"]
  #verify categorical variables are still categorical (contents function?)
  #contents(clin_full_wImputedRF1[varlist])


###############################################################################
################ Combine abundance with clinical data #########################
###############################################################################

########## Using clinical data WITHOUT imputations ##########

#Data that includes serum runs 1 and 2
comboD1D2 <- merge(clin24_restricted_clean, resp_D1D2, by=c("code","Study"), all=F)
dim(comboD1D2)[1] #number of rows in resulting data -- a perfect match

#Data for serum run 1 only
comboD1_filter50n <- merge(clin24_restricted_clean, respD1_filter50n, by=c("code","Study"), all=F) #88 obs
save(comboD1_filter50n, file=paste(outputsDir,"comboD1_filter50n.RData", sep=""))
#Data for serum run 2 only
comboD2_filter50n <- merge(clin24_restricted_clean, respD2_filter50n, by=c("code","Study"), all=F) #75 obs

#Data for Nicaragua saliva
temp <- merge(clin24_restricted_clean, respD3_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n = temp[which(temp$code!="ID1408"),]  #85
table(comboD3_filter50n$DxFinal4cat) #cannot do DF vs DHF/DSS analysis with only 7 DHF/DSS (!)
save(comboD3_filter50n, file=paste(outputsDir,"comboD3_filter50n.RData", sep=""))
comment(comboD3_filter50n) <- "Data contains Nicaragua saliva LCMS (Mass Hunter data from Natalia) and corresponding clinical data.  
      Nothing has been imputed (contains missings)."

#Data for Nicaragua urine
temp <- merge(clin24_restricted_clean, respD5_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID and those which were duplicated in this data
comboD5_filter50n = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                                 temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
table(comboD5_filter50n$DxFinal4cat) #cannot do DF vs DHF/DSS analysis with only 7 DHF/DSS (!)
save(comboD5_filter50n, file=paste(outputsDir,"comboD5_filter50n.RData", sep=""))
comment(comboD5_filter50n) <- "Data contains Nicaragua urine LCMS (Mass Hunter data from Natalia) and corresponding clinical data.  
      Nothing has been imputed (contains missings)."


########## Using clinical data WITH imputations ##########

#Data that includes serum runs 1 and 2
comboD1D2_wImpRF1 <- merge(clin_full_wImputedRF1, resp_D1D2, by=c("code","Study"), all=F)

#Data for serum run 1 only
comboD1_filter50n_wImpRF1 <- merge(clin24_full_wImputedRF1, respD1_filter50n, by=c("code","Study"), all=F) #88 obs
save(comboD1_filter50n_wImpRF1, file=paste(outputsDir,"comboD1_filter50n_wImpRF1.RData", sep=""))
#Data for serum run 2 only
comboD2_filter50n_wImpRF1 <- merge(clin_full_wImputedRF1, respD2_filter50n, by=c("code","Study"), all=F) #75 obs

#Data for Nicaragua saliva
temp <- merge(clin_full_wImputedRF1, respD3_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n_wImpRF1 = temp[which(temp$code!="ID1408"),]  #85
save(comboD3_filter50n_wImpRF1, file=paste(outputsDir,"comboD3_filter50n_wImpRF1.RData", sep=""))
#Data for Nicaragua urine
temp <- merge(clin_full_wImputedRF1, respD5_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID and those which were duplicated in this data
comboD5_filter50n_wImpRF1 = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                                 temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
save(comboD5_filter50n_wImpRF1, file=paste(outputsDir,"comboD5_filter50n_wImpRF1.RData", sep=""))

