
##############################################################################
################ additional investigation of abundance data ##################
##############################################################################


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