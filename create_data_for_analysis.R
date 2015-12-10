###############################################################################
####################### Create data to be used in analysis ####################
###############################################################################
# Refer to main_code.R to establish directories and functions
# See "Guide to Dengue Data.xls" for documentation on data used in this file


###############################################################################
#################### Clean the abundance data #################################
###############################################################################

#for Nicaragua Serum samples.  Note that the first batch data seems to be of higher quality.  
  #Might not want to publish results from batches 3 or 4 ("run 2") 
respD1_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_first batch.txt", lcms_run=1, printme=TRUE) #88 samples
respD1_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_first batch.txt", lcms_run=1, printme=TRUE)
respD2_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples
respD2_filter10n = clean_LCMS(infile="LCMS_serum_Nica_10percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples

#for non-invasive Nicaraguan samples (urine and saliva)
temp = clean_LCMS(infile="LCMS_saliva_Nica_50percent.txt", lcms_run=3, printme=TRUE) #86 (as in LCMS excel)
#remove observation that probably has incorrect ID
respD3_filter50n = temp[which(temp$code!="ID1408"),]  #85
write.table(respD3_filter50n, paste(outputsDir,"respD3_filter50n.txt", sep=""), sep="\t")

temp = clean_LCMS(infile="LCMS_urine_Nica_50percent.txt", lcms_run=5, printme=TRUE) #91 (as in LCMS excel)
#remove observation that probably has incorrect ID and those which were duplicated in this data
respD5_filter50n = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                                         temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
write.table(respD5_filter50n, paste(outputsDir,"respD5_filter50n.txt", sep=""), sep="\t")


# Combine the abundance data (runs 1 and 2) 
#list the columns (mz values) that are in common between serum runs 1 and 2
commonvars = intersect(colnames(respD2_filter50n),colnames(respD1_filter50n)) #only 40 compound identifiers in common (10% data)
resp_D1D2 = rbind(respD1_filter50n[,commonvars], respD2_filter50n[,commonvars])

#Get list of patient ID codes that are present in LCMS data (Nicaragua serum R1 and R2, urine, and saliva)
IDs_in_resp_all = get_IDs_in_common(resp_D1D2, respD3_filter50n, respD5_filter50n)
save(IDs_in_resp_all, file=paste(outputsDir,"IDs_in_resp_all.RData", sep="")) #for future access
load(paste(outputsDir,"IDs_in_resp_all.RData", sep="")) #loads IDs_in_resp_all

###############################################################################
########################### Clean the clinical data ###########################
###############################################################################

# this function was established in "clean_clinical_data_functions_v2.R"
clin24_full_clean = clean_clin_initial_data(clinic_varsD, IDs_in_resp_all, time_period=24) 
clin12_full_clean = clean_clin_initial_data(clinic_varsD, IDs_in_resp_all, time_period=12) #5876 observations (no restriction based on LCMS presence)

#keep only clinical observations for which we have LC-MS data
clin24_restricted_clean = merge(IDs_in_resp_all, clin24_full_clean, by=c("code","Study"), all=FALSE)

#keep only clinical observations for batch 1 LCMS data
clin24_D1_clean = merge(respD1_filter50n[,c("code","Study")], clin24_full_clean, by=c("code","Study"), all=FALSE)
clin12_D1_clean = merge(respD1_filter50n[,c("code","Study")], clin12_full_clean, by=c("code","Study"), all=FALSE)
clin12_D3_clean = merge(respD3_filter50n[,c("code","Study","Cod_Nin")], clin12_full_clean, by=c("code","Study","Cod_Nin"), all=FALSE)
clin12_D5_clean = merge(respD5_filter50n[,c("code","Study","Cod_Nin")], clin12_full_clean, by=c("code","Study","Cod_Nin"), all=FALSE)


## Write this clinical data to file for easy future access ##
save(clin24_full_clean, file=paste(outputsDir,"clin24_full_clean.RData", sep=""))
save(clin12_full_clean, file=paste(outputsDir,"clin12_full_clean.RData", sep=""))
write.table(clin24_full_clean, paste(outputsDir,"clin24_full_clean.txt", sep=""), sep="\t")
write.table(clin12_full_clean, paste(outputsDir,"clin12_full_clean.txt", sep=""), sep="\t")
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
                              XD=clin12_full_clean, restrict_to_cohort_vars=F, 
                              restrict_to_hospit_vars=T,UltraX=T, BloodLab=T)
#note: imputation uses only vars that are never missing -- for cohort + hospital this leaves 10 vars
  # currently I use cohort data only as a test set
  # for hospit data, imput missings without knowledge of cohort data values (this is fair and proper)
  # for cohort data, imput missings with the help of the hospit data 
      #(predictions for cohort are conditional on hospit data being available so I think this is not unreasonable)
      #we could impute cohort data without using hospit data but then all non-cohort vars will be missing (not terrible)
#impute missing values and add indicator for imputation - all data (we will take subset for cohort)
prelim_cohort = impute_missings(XD=clin12_full_clean, vars_to_impute=varlist, 
                                                 predictor_vars_prelim=varlist, exclude_miss_vars=T, method="RF")
prelim_cohort = prelim_cohort[which(prelim_cohort$Study=="Cohort"),]
prelim_hospit = impute_missings(XD=clin12_full_clean[which(clin12_full_clean$Study=="Hospital"),], 
                        vars_to_impute=varlist, predictor_vars_prelim=varlist, exclude_miss_vars=T, method="RF")
#drop the imputation dummies that are in cohort but not hospit data (model to to hospit will not include these anyhow)
c_h = colnames(prelim_cohort)[colnames(prelim_cohort) %in% colnames(prelim_hospit)==TRUE] #list of vars to include
prelim = rbind(prelim_cohort[c_h], prelim_hospit)
#note: since missings were imputed for categoricals but not binary equivalents,
#must create these binary vars now rather than using the ones already in data
prelim = create_binary_variables(prelim)
clin12_full_wImputedRF1 = create_hospit_binary_vars(prelim)
#write to file
save(clin12_full_wImputedRF1, file=paste(outputsDir,"clin12_full_wImputedRF1.RData", sep=""))
write.table(clin12_full_wImputedRF1, paste(outputsDir,"clin12_full_wImputedRF1.txt", sep=""), sep="\t")


#now do the imputations using just hospital data (to use all 40-something non-missing vars)
clin12_hospit_clean = clin12_full_clean[which(clin12_full_clean$Study=="Hospital"),]
varlist = get_clinic_var_list(clinic_varsD, outcome="either", 
                              eliminate_vars_with_missings=F, eliminate_constant_vars=T, 
                              eliminate_vars_with_minXnomiss=50,
                              XD=clin12_hospit_clean, restrict_to_cohort_vars=F, 
                              restrict_to_hospit_vars=T,UltraX=T, BloodLab=T)
#impute missing values and add indicator for imputation 
clin12_hospit_clean$PCR = droplevels(clin12_hospit_clean$PCR) #get rid of empty factor levels (will cause error)
clin12_hospit_wImputedRF1_prelim = impute_missings(XD=clin12_hospit_clean, vars_to_impute=varlist, 
                                                 predictor_vars_prelim=varlist, exclude_miss_vars=T, method="RF")
#note: since missings were imputed for categoricals but not binary equivalents,
#must create these variables now rather than using the ones already in data
clin12_hospit_wImputedRF1 = create_binary_variables(clin12_hospit_wImputedRF1_prelim)
clin12_hospit_wImputedRF1 = create_hospit_binary_vars(clin12_hospit_wImputedRF1)
#write to file
save(clin12_hospit_wImputedRF1, file=paste(outputsDir,"clin12_hospit_wImputedRF1.RData", sep=""))
write.table(clin12_hospit_wImputedRF1, paste(outputsDir,"clin12_hospit_wImputedRF1.txt", sep=""), sep="\t")



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

#Data for Nicaragua saliva -- might need to redo with Cod_Nin taken into account
temp <- merge(clin24_restricted_clean, respD3_filter50n, by=c("code","Study"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n = temp[which(temp$code!="ID1408"),]  #85
table(comboD3_filter50n$DxFinal4cat) #cannot do DF vs DHF/DSS analysis with only 7 DHF/DSS (!)
save(comboD3_filter50n, file=paste(outputsDir,"comboD3_filter50n.RData", sep=""))
comment(comboD3_filter50n) <- "Data contains Nicaragua saliva LCMS (Mass Hunter data from Natalia) and corresponding clinical data.  
      Nothing has been imputed (contains missings)."

#Data for Nicaragua urine -- might need to redo with Cod_Nin taken into account
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
comboD1_filter50n_wImpRF1 <- merge(clin12_full_wImputedRF1, respD1_filter50n, by=c("code","Study"), all=F) #88 obs
save(comboD1_filter50n_wImpRF1, file=paste(outputsDir,"comboD1_filter50n_wImpRF1.RData", sep=""))
write.table(comboD1_filter50n_wImpRF1, paste(outputsDir,"comboD1_filter50n_wImpRF1.txt", sep=""), sep="\t")

#Data for serum run 2 only
comboD2_filter50n_wImpRF1 <- merge(clin_full_wImputedRF1, respD2_filter50n, by=c("code","Study"), all=F) #75 obs

#Data for Nicaragua saliva
temp <- merge(clin12_full_wImputedRF1, respD3_filter50n, by=c("code","Study","Cod_Nin"), all=F) 
#remove observation that probably has incorrect ID
comboD3_filter50n_wImpRF1 = temp[which(temp$code!="ID1408"),]  #85
save(comboD3_filter50n_wImpRF1, file=paste(outputsDir,"comboD3_filter50n_wImpRF1.RData", sep=""))
write.table(comboD3_filter50n_wImpRF1, paste(outputsDir,"comboD3_filter50n_wImpRF1.txt", sep=""), sep="\t")

#Data for Nicaragua urine
temp <- merge(clin12_full_wImputedRF1, respD5_filter50n, by=c("code","Study","Cod_Nin"), all=F) 
#remove observation that probably has incorrect ID and those which were duplicated in this data
comboD5_filter50n_wImpRF1 = temp[which(temp$code!="ID1408" & temp$code!="ID1304" & temp$code!="ID1315" & 
                                 temp$code!="ID1342" & temp$code!="ID1348" & temp$code!="ID1380"),]  #80
save(comboD5_filter50n_wImpRF1, file=paste(outputsDir,"comboD5_filter50n_wImpRF1.RData", sep=""))
write.table(comboD5_filter50n_wImpRF1, paste(outputsDir,"comboD5_filter50n_wImpRF1.txt", sep=""), sep="\t")

### for Kristofs batch 1 patients ###
RP_batch1 = read.csv(paste(clinical_inputsDir,"RP_batch1_sample_merge_info.csv", sep=""), header=TRUE, nrows=1000)
names(RP_batch1)[names(RP_batch1)==c("Code")] = c("code") #rename column
RP_batch1$code = apply(X=RP_batch1, MARGIN = 1, FUN=fcode, var="code")
RP_batch1$Cod_Nin = tolower(RP_batch1$Sample)
RP_batch1[which(RP_batch1$Study=="Hospital"),c("Cod_Nin")] = NA
clinical_RPb1_clean <- merge(clin12_full_wImputedRF1, RP_batch1, by=c("code","Study","Cod_Nin"), all=F)
save(clinical_RPb1_clean, file=paste(outputsDir,"clinical_RPb1_clean.RData", sep=""))
write.table(clinical_RPb1_clean, paste(outputsDir,"clinical_RPb1_clean.txt", sep=""), sep="\t")

