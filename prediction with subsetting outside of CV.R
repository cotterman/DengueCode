if(reduce_in_CV==F){
  
  ############# choose variable subset using all data (not within CV) #########
  
  #obtain features that demonstrate at least a 2-fold differentiation based on t-test with a p-value<.05.
  #Take the 10-15 such features with the lowest p-values.
  ttest_top5MFs_noC = row.names(ttest_select(XD, topX=5, use_clinical=FALSE))
  ttest_top10MFs_noC = row.names(ttest_select(XD, topX=10, use_clinical=FALSE))
  
  #[optional] Add 1 and then transform to the log base 2 scale.
  #XD[cindices] = log2(XD[cindices]+1)
  
  #Choose a subset of the Xx features using “best subsets selection” 
  #to identify the logistic regression model with the smallest AIC 
  #(bestglm will also implement CV rather than AIC, which is a superior method). 
  #The argument Xy should be a data-frame containing the design matrix in the first columns
  #and the outcome in the last column
  #note: this regression cannot be run when all clinical variables are included.
  #infact, seemed to have trouble when using more than 10 MFs in input
  set.seed(922) #this will matter if we do CV to find best subset
  modelD = cbind(XD[,c(paste(ttest_top10MFs_noC, sep=","))],
                 XD[,"is.DEN"])
  mybestglm = bestglm(Xy=modelD, family=binomial, IC="AIC", 
                      TopModels = 1, method = "exhaustive", intercept = TRUE, 
                      nvmax = "default")
  print("Logistic regression using best subset selection:")
  print(mybestglm)
  #mybestglm$Subsets #for each subset size, view which variables are selected
  #mybestglm$BestModel
  bestLogit_MFs = colnames(mybestglm$BestModel$model)[-1]
  bestLogitCount = length(bestLogit_MFs)
  
  #Use random forest variable importance to determine subset of Xx features
  #cutoff seems to be the parameter to toggle in order to develop ROC:  
  #This parameter dictates the fraction of votes that win for classification.  
  #Ex: if set to c(.25, .75) then observation will only be classified as diseased if it gets 75% of the tree votes
  #(or, I can just take the predicted probabilities from rf output and do my own ROC calc as usual)
  #sampsize seems to give us the option of sampling within each strata (to gaurantee we get both 0s and 1s in each bootstrap sample)
  #should I set this?
  #Importance: set to TRUE if variable importance calculation using Gini coefficients is desired
  #localImp: set to TRUE if permutation-based variable importance calculation is desired
  set.seed(102)
  rf = randomForest(x=XD[,c(allMFs)], y=as.factor(XD[,"is.DEN"]), ntree=501, localImp=FALSE)
  rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
  print("Top 5 VIPS from random forest")
  rfVIP5 = names(rf_VIPS[order(rf_VIPS, decreasing=TRUE)[1:5],])
  print(rfVIP5) #take a look at what rf chose
  rfVIP4 = names(rf_VIPS[order(rf_VIPS, decreasing=TRUE)[1:4],])
  rfVIP3 = names(rf_VIPS[order(rf_VIPS, decreasing=TRUE)[1:3],])
  rfVIP2 = names(rf_VIPS[order(rf_VIPS, decreasing=TRUE)[1:2],])
  rfVIP1 = rfVIP5[1]
  
  #Use the lasso to determine subset of Xx features
  #cv.glmnet chooses the optimal tuning parameter based on CV.
  #type.measure could also be set to "auc", which would find the lambda that minimizes the area under the ROC
  #I do not see documentation which allows penality.factor to be used with cv.glmnet, 
  #so I might need to find the optimal lamba myself using glmnet with CV
  set.seed(456)
  cvLasso = cv.glmnet(x=as.matrix(XD[cindices]),
                      y=as.matrix(XD[,"is.DEN"]),
                      family='binomial', alpha=1, nfolds=10, type.measure="deviance")
  plot(cvLasso) #shows deviance as a function of the tuning parameter lambda
  #lambda.1se if the largest value of lambda such that error is within 1 std err of the minimum
  #lambda.min is the value of lambda that gives the minimum  mean CV error
  cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
  LassoVIP = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs (MFs with non-zero coef)
  LassoCount = length(LassoVIP)
  print("MFs with non-zero coefficients from lasso")
  print(LassoVIP)
  
  
  ## Methods from the literature ##
  
  #The "Mahapatra" method
  #note: must first take log transform of data, as she did
  p_bestLogit = run_CV_and_report(outD, XD, "logit", bestLogit_MFs, "MFs", "t_AIC", bestLogitCount, F) 
  
  #The "Denery" method 
  outD = run_CV_and_report_bound(outD, "RF", ttest_top5MFs_noC, "MFs", "t-test",5,F)
  outD = run_CV_and_report_bound(outD, "RF", c(ttest_top5MFs_noC, cvars_noMiss), "all", "t-test",5,F)
  
  ## Cotterman methods ##
  
  #No subsetting --- using all info with random forest
  outD = run_CV_and_report_bound(outD, "RF", cvars_noMiss, "non-MFs", "none",clinicCount,F)
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "none",MFcount,F)
  outD = run_CV_and_report_bound(outD, "RF", c(allMFs, cvars_noMiss), "all", "none",allCount,F)
  
  #Subset using RF variable importance measure, then random forest on selected variables
  outD = run_CV_and_report_bound(outD, "RF", rfVIP5, "MFs", "RF",5,F)
  outD = run_CV_and_report_bound(outD, "RF", c(rfVIP5, cvars_noMiss), "all", "RF",5,F)
  
  #Subset using lasso, then random forest on selected variables
  outD = run_CV_and_report_bound(outD, "RF", LassoVIP, "MFs","lasso",LassoCount,F)
  outD = run_CV_and_report_bound(outD, "RF", c(LassoVIP, cvars_noMiss), "all","lasso",LassoCount,F)
  
  #No subsetting --- using all info with super learner
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "none",MFcount,F)
  outD = run_CV_and_report_bound(outD, "SL", c(allMFs, cvars_noMiss), "all", "none",allCount,F)
  
  #Subset using RF variable importance measure, then super learner on selected variables
  outD = run_CV_and_report_bound(outD, "SL", rfVIP5, "non-MFs", "RF",5,F) 
  outD = run_CV_and_report_bound(outD, "SL", c(rfVIP5, cvars_noMiss), "all", "RF",5,F)
  outD = run_CV_and_report_bound(outD, "SL", c(rfVIP4, cvars_noMiss), "all", "RF",4,F)
  outD = run_CV_and_report_bound(outD, "SL", c(rfVIP3, cvars_noMiss), "all", "RF",3,F)
  outD = run_CV_and_report_bound(outD, "SL", c(rfVIP2, cvars_noMiss), "all", "RF",2,F)
  outD = run_CV_and_report_bound(outD, "SL", c(rfVIP1, cvars_noMiss), "all", "RF",1,F)
  
  #Subset using lasso, then super learner on selected variables
  outD = run_CV_and_report(outD,data_numeric, "SL", c(LassoVIP), "MFs", "lasso",LassoCount,F)
  outD = run_CV_and_report(outD,data_numeric, "SL", c(LassoVIP, cvars_noMiss), "all", "lasso",LassoCount,F)
  
}