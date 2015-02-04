
library(multiPIM)
library(tmle) #works only for binary treatments
library(tmle.npvi) #works for continuous treatments but must provide 

get_VIM = function(XD, method, suffix, clinic_vars){

  #impute missing values using RF proximity measure (see RandomForest documentation)
  set.seed(101)
  XD_missimput = rfImpute(x=XD[,c(clinic_vars)], y=XD[,"DEN_dum"], ntree=500)
  colnames(XD_missimput) = c("DEN_dum",clinic_vars)
  
  #RF variable importance
  if(method=="rfVIM"){
    set.seed(102)
    rf = randomForest(x=XD_missimput[,c(clinic_vars)], y=XD_missimput[,"DEN_dum"], ntree=500, localImp=FALSE)
    rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
    rfVIP = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
    myVIM = as.data.frame(cbind(names(rfVIP), rfVIP))
    colnames(myVIM) = c("VarName",paste("rfVIM",suffix,sep="_"))
  }
  #convert variables into numerics
  YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                             "Variable.Name.in.Hospital.data"]
  YesNo_varsH_char = sapply(YesNo_varsH, as.character)
  XD_tmle = convert_factor_to_numeric(XD_missimput, c("DEN_dum","Torniquete","Sexo",YesNo_varsH_char))
  
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

#obtain list of clinical variables to include (will impute when missing) 
clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                          "Variable.Name.in.Hospital.data"]),"age","DaysSick") #43+2
print(clinic_vars)
clinicCount = length(clinic_vars)
print(clinicCount)

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
