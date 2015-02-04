
###################### Apply simple learners to Dengue data ###################

###############################################################################
######################## Obtain list of predictors  ###########################
###############################################################################

#list of clinical variables
DENvarlist_all = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                                    eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                    XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
DENvarlist_CohortRestrict = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                                               eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                               XD=clinical_full_clean, restrict_to_cohort_vars=T, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
DENvarlist_noUltraX = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                                         eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                         XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=F, BloodLab=T)
DENvarlist_genOnly = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                                        eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                        XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=F, BloodLab=F)
DENvarlist_noMiss = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=T, 
                                       eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                       XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=F, UltraX=T, BloodLab=T)
DENbinary_genOnly = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                                        eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50, restrict_to_binary=T,
                                        XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=F, BloodLab=F)

#list of MFs 
cindices = grep("MZ_",colnames(dataf)) #column indices in dataset
allMFs = colnames(dataf[cindices]) #names
MFcount = length(allMFs)


###############################################################################
############################ Ordinary CART ####################################
###############################################################################

# [probably don't care to present this since ctree is better --- this doesn't use data to determine parameters.]
XD = clin_full_wImputedRF1
fit.rpart=rpart(is.DEN ~ ., data=XD[,c("is.DEN",DENvarlist_genOnly)], method="class"
                control = rpart.control(cp = 0.001, minsplit = 40, xval = 20, minbucket = 20, maxdepth=10), 
                method = "class") #, cost = rep(0.000001,xp)
fit.rpart #view list of splits
plot.new()
plot(fit.rpart) #pretty ugly plot!!
text(fit.rpart, use.n=T)
summary(fit.rpart)


###############################################################################
############################ CART Alternative #################################
###############################################################################

fit.ctree <- ctree(is.DEN~., data=clin_full_wImputedRF1[,c("is.DEN",DENvarlist_genOnly)])
fit.ctree <- ctree(is.DEN~., data=clin_full_wImputedRF1[,c("is.DEN","Sexo","age","Temperatura")])
png(paste(resultsDir,"ctree_clinOnly_genOnly.png", sep=""), width=800, height=500)
plot(fit.ctree) #why does fit look perfect and splits unhelpful???
dev.off()

###############################################################################
#################### Oridinary Logic regression ###############################
###############################################################################

XD_binary = create_binary_variables(clin_full_wImputedRF1)
#list of newly crated variables that are transformations of continuous and categorical "general symptoms/hem"
#use these in addition to already binary "general symptoms/hem" vars in logic regression
vars_binary_genOnly = c("is.LiverEnlarged", "is.hypertension","is.hypotension","is.braquicardia","is.fever","is.resp_normal",
                        "is.torniquete10plus","is.torniquete20plus","is.pulse_rapid","is.pulse_strong",
                        "is.age_under1","is.age_1or2","is.age_3or4","is.age_5to8","is.age_9to12","is.age_13to15",
                        "is.DaysSick_under4","is.DaysSick_4to6","is.DaysSick_7plus")
#convert factors into logicals
for (name in DENbinary_genOnly) {
  #print(name)
  XD_binary$temp = (XD_binary[name]=="Yes" | XD_binary[name]=="female")
  #drop old version and rename new version
  XD_binary[name] = NULL
  XD_binary[name] = XD_binary["temp"]
  XD_binary["temp"] = NULL
}
vars_binary = c(vars_binary_genOnly, DENbinary_genOnly)
vars_binary
fit.logreg = logreg(XD_binary$is.DEN, XD_binary[,c(vars_binary)], select=1)
fit.logreg #see the sentence
plot.new()
par(mfrow=c(1,1))
png(paste(resultsDir,"logreg_clinOnly_genOnlyBinary.png", sep=""), width=960, height=600)
plot(fit.logreg) #sentence in graph form
dev.off()

