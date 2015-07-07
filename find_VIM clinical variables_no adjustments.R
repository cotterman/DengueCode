
library(multiPIM)
library(tmle) #works only for binary treatments
library(tmle.npvi) #works for continuous treatments but must provide reference value for X (or assumed 0)
    #and I think it is estimating a linear fit for VIM, which is not cool
  #I could force continuous variables to be on scale from 0 to 1....should find out more about how reference value of X is used in tmle.
  #I could also turn to binary by giving value of 0/1 if less than/greater than 75th percentile (or something)
  #For now, I could just do the variable importance stuff for general signs & symptoms, which I have already dichotomized
    #or, forget about tMLE altogether

#data to use
myhospit = clin12_full_wImputedRF1[which(clin12_full_wImputedRF1$Study=="Hospital"),]
#obtain list of clinical variables to include
clinic_vars = get_clinic_var_list(clinic_varsD, outcome="either", eliminate_vars_with_missings=F, 
                                  eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                                  XD=myhospit, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
print(clinic_vars)
print(length(clinic_vars))
#Here is an opportunity to add imputation dummies (include_imp_dum="all")
  #this function is created in "prediction_functions_v2.R"
  #It will also eliminate IR, PCR for ND.vs.DEN
myvars = get_clinical_predictors("ND.vs.DEN", clinic_varsD, clinic_vars, myhospit, include_imp_dums="all")
#make outcome a factor so RF will do classification
myhospit[which(myhospit$is.DEN==TRUE), "DEN_dum"] = "1"
myhospit[which(myhospit$is.DEN==FALSE), "DEN_dum"] = "0"
myhospit$DEN_dum = factor(myhospit$DEN_dum)
rf = randomForest(x=myhospit[,c(clinic_vars)], y=myhospit[,"DEN_dum"], ntree=500, nPerm=1, importance=TRUE)
VIM_OOB = importance(rf, type=1) #"mean decrease accuracy" from permuting variable in oob data, divided by SD
oob = as.data.frame(cbind(as.numeric(VIM_OOB), rownames(VIM_OOB) ))
VIM_Gini = importance(rf, type=2) #mean decrease Gini 
gini = as.data.frame(cbind(as.numeric(VIM_Gini), rownames(VIM_Gini) ))
VIM_rf = merge(oob, gini, by="V2")
colnames(VIM_rf) = c("varname","RF_OOB","RF_Gini")
#sort -- grr.....values are interpreted as chars not numerics for sorting.  strange.
VIM_rf = VIM_rf[order(VIM_rf$RF_OOB),]
plot(VIM_rf$RF_OOB, VIM_rf$RF_Gini)
#now export so we can read it into python and make graphs with other VIMs
write.table(VIM_rf, paste(outputsDir,"VIM_rf.txt", sep=""), sep="\t")




###############################################################################
############################## old analyses ###################################
###############################################################################

get_VIM = function(XD, method, suffix, clinic_vars, rf_type){
  
  #impute missing values using RF proximity measure (see RandomForest documentation)
  set.seed(101)
  XD_missimput = rfImpute(x=XD[,c(clinic_vars)], y=XD[,"DEN_dum"], ntree=500)
  colnames(XD_missimput) = c("DEN_dum",clinic_vars)
  
  #RF variable importance
  if(method=="rfVIM"){
    set.seed(102)
    rf = randomForest(x=XD_missimput[,c(clinic_vars)], y=XD_missimput[,"DEN_dum"], ntree=500, localImp=FALSE)
    if(rf_type=="oob"){
      rf_VIPS = importance(rf, type=1) #type=1 uses the permutation method.  
    }
    if(rf_type=="gini"){
      rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  
    }
    rfVIP = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
    myVIM = as.data.frame(cbind(names(rfVIP), rfVIP))
    colnames(myVIM) = c("VarName",paste("rfVIM",suffix,sep="_"))
  }
  #convert variables into numerics
  YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                             "Variable.Name.in.Hospital.data"]
  YesNo_varsH_char = sapply(YesNo_varsH, as.character)
  XD_tmle = convert_factor_to_numeric(XD_missimput, c("is.DEN","Torniquete","Sexo",YesNo_varsH_char))
  
  #TMLE variable importance
  set.seed(102)
  #obtain marginal treatment effects of each dichotomized point treatment
  #tmleOut = tmle(Y=XD_tmle[,"DEN_dum"], A=XD_tmle[,"Sexo"], W=XD_tmle[,c("Cefalea","Dolor_ocular","Rash")], 
  #               Q.SL.library = c("SL.mean", "SL.knn", "SL.randomForest"), family="binomial")
  #summary(tmleOut)
  #tmleOut$estimates$RR$psi #estimate of relative risk (could alternatively request OR or ATE)
  #tmleOut$estimates$RR$pvalue #corresponding pvalue
  
  #Lasso variable importance
  if(method=="lasso"){
    set.seed(102)
    cvLasso = cv.glmnet(x=as.matrix(XD_tmle[,c(clinic_vars)]),
                        y=as.matrix(XD_tmle[,"DEN_dum"]),
                        family='binomial', alpha=1, nfolds=10, type.measure="auc")
    cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
    MFpreds = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs
    myVIM = as.data.frame(cbind(rownames(cvLcoefs), cvLcoefs))
    colnames(myVIM) = c("VarName",paste("Lasso",suffix,sep="_"))
  }
  
  return(myVIM)
}


#find most important variables
rfVIM_n88 = get_VIM(XD=comboD1_filter50n, method="rfVIM", suffix="n88", clinic_vars)
rfVIM_n1624 = get_VIM(XD=clinical_prelim, method="rfVIM", suffix="n1624", clinic_vars)
lassoVIM_n88 = get_VIM(XD=comboD1_filter50n, method="lasso", suffix="n88", clinic_vars)
lassoVIM_n1624 = get_VIM(XD=clinical_prelim, method="lasso", suffix="n1624", clinic_vars)

#combine variable importance measures into 1 chart
VIM_merged = merge(rfVIM_n88, rfVIM_n1624,  by="VarName", all=TRUE)
VIM_merged = merge(VIM_merged, lassoVIM_n88,  by="VarName", all=TRUE)
VIM_merged = merge(VIM_merged, lassoVIM_n1624,  by="VarName", all=TRUE)
VIM_merged
write.csv(x=VIM_merged, file=paste(resultsDir,"VIM_merged.txt"), row.names = FALSE)
