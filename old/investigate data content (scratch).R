
####### Data checks (pre main_code) #########



############### Basic cohort info #####################

B1 = read.delim(paste(inputsDir,"Basic info_cohort study.txt", sep=""), header=TRUE, nrows=100)
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
#need to take the characters up till the first dot, and then reformat
formatID = function(x) sprintf("ID%04d",as.numeric(strsplit(as.character(x["Sample"]), "\\.")[[1]][1:1]))
newIDs_3 = apply(X=B1, MARGIN = 1, FUN=formatID)
B1["code"] = newIDs_3

B2 = read.delim(paste(inputsDir,"Lab_cohort_noninvasive samples.txt", sep=""), header=TRUE, nrows=100)
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
#need to take the characters up till the first dot, and then reformat
newIDs_3 = apply(X=B2, MARGIN = 1, FUN=formatID)
B2["code"] = newIDs_3

get_mismatches(B1[,"code"], B2[,"code"], "B1","B2") #B2 has every ID that B1 has except for 1584 (so I added this one by hand)



############### Clinical cohort info #####################

C1 = read.delim(paste(inputsDir,"Clinical data_cohort study.txt", sep=""), header=TRUE, nrows=100)
formatID = function(x) sprintf("ID%04d",as.numeric(strsplit(as.character(x["Sample"]), "\\.")[[1]][1:1]))
newIDs_4 = apply(X=C1, MARGIN = 1, FUN=formatID)
C1["code"] = newIDs_4

C2 = read.delim(paste(inputsDir,"Clinical_cohort_noninvasive samples.txt", sep=""), header=TRUE, nrows=100)
newIDs_4 = apply(X=C2, MARGIN = 1, FUN=formatID)
C2["code"] = newIDs_4

get_mismatches(C1[,"code"], C2[,"code"], "C1","C2")


############# compare clinical to lab ################

get_mismatches(B1[,"code"], C1[,"code"], "B1","C1")
get_mismatches(B2[,"code"], C2[,"code"], "B2","C2")


############# compare clinical and lab to LCMS ###############

get_mismatches(B2[,"code"], respD3_filter50n[,"code"], "B2","resp_saliva")
get_mismatches(B2[,"code"], respD5_filter50n[,"code"], "B2","resp_urine")


