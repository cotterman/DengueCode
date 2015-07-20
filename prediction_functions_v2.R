#v2 takes out code that selects variables outside of CV (do not have reason to run this again)


#obtain features that demonstrate at least a 2-fold differentiation based on t-test with a p-value<.05.  
#note: t-test default is to set var.equal=FALSE, which means using the Welch test (good)
get_ttest_results <- function(mydata){
  myttest <- function(x){
    ttest = t.test(x[which(mydata$is.DEN==T)],
                   x[which(mydata$is.DEN==F)])
    fold.diff <- ttest$estimate[[1]]/ttest$estimate[[2]]
    return(c(fold.diff, ttest$p.value))
  }
  MZindices = grep("MZ_",colnames(mydata))
  ttest_results = data.frame(t(sapply(mydata[MZindices], FUN=myttest)))
  colnames(ttest_results) = c("fold_diff", "pval")
  head(ttest_results)
  tsig = ttest_results[which( (ttest_results[,"fold_diff"]>2 | ttest_results[,"fold_diff"]<.5) & 
                                ttest_results[,"pval"]<.05),]
  return(tsig)
}


#Take the 10-15 such features with the lowest p-values.
ttest_select = function(tdata, topX=5, use_clinical) {
  tsigs = get_ttest_results(tdata) 
  #the compounds with the topX lowest p-value, provided there is at least a 2 fold change
  sigX = tsigs[order(tsigs$pval, tsigs$fold_diff)[1:topX],]
  #print(paste("Selected", topX, "compounds based on t-test:"))
  #print(row.names(sigX))
  return(sigX)
}


#calulate the area under the curve determined by list of x,y coordinates
   #uses the trapezoid rule
myauc = function(x,y){
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2) 
}


##Obtain sensitivity, specificity, overall error rate and area under the ROC
report_performance = function(pdata) {
  
  measure_performance = function(mydata, threshold){
    #change the probability threshold for 0 vs. 1 classifications in order to trace ROC curve?
    mydata$is.DENpred = NA
    mydata$is.DENpred[mydata$DEN_prob<threshold] = F 
    mydata$is.DENpred[mydata$DEN_prob>=threshold] = T 
    mydata[,c("is.DEN","DEN_prob","is.DENpred")] #view actual versus predicted
    #Sensitivity (also called the true positive rate) measures the proportion of actual positives 
    #which are correctly identified as such (e.g. the percentage of sick people who are correctly identified as having the condition). 
    sensitivity = round(sum(mydata$is.DEN==T & mydata$is.DENpred==T)/sum(mydata$is.DEN==T),3)
    sensitivity
    #Specificity measures the proportion of negatives which are correctly identified as such 
    #(e.g. the percentage of healthy people who are correctly identified as not having the condition, 
    #sometimes called the true negative rate).
    specificity = round(sum(mydata$is.DEN==F & mydata$is.DENpred==F)/sum(mydata$is.DEN==F),3)
    specificity
    #Overall error rate
    pred_err = round(sum(mydata$is.DEN != mydata$is.DENpred)/sum(!is.na(mydata$is.DEN)),3)
    pred_err
    #ROC curves plot the false positive rate (1-specificity) against the true positive rate (sensitivity)
    return(c(sensitivity, specificity, pred_err))
  }
  performance = data.frame(t(sapply(X=seq(1, 0, -.05), FUN=measure_performance, mydata=pdata)))
  colnames(performance) = c("sensitivity", "specificity", "pred_err")
  #write.csv(x=performance, file="/srv/scratch/carolyn/Results/data_for_ROC.txt", row.names = FALSE) #data to plot for ROC curve
  #ROC plot (uncomment if desired)
  #plot(x=(1-performance$specificity), y=performance$sensitivity, 
  #     col="grey", lwd=2, type="l",
  #     ylab="true positive rate", xlab="false positive rate",main="ROC curve with 10-fold cross validation")
  #specificity, sensitivity, and error rates when threshold is chosen so as to minimize error rate
  report = performance[which(performance$pred_err==min(performance$pred_err)),][1,] #chose the first when there is a tie
  #calculate area under the ROC curve (does not work ever since installing the glmnet package --- asked MESS package auc function)
  report$ROC_area = round(myauc(x=(1-performance$specificity), y=performance$sensitivity),3)
  #calculate mean squared error (using actual outcome and predicted probability)
  report$MSE = round(mean( (as.vector(pdata$DEN_prob) - as.numeric(pdata$is.DEN))^2 ),3)
  report$SD_SqEr = round(sd( (as.vector(pdata$DEN_prob) - as.numeric(pdata$is.DEN))^2 ),3)
  return(report)
}

#same as built-in SL.knn, but convert all variables into numerics
my.SL.knn <- function (Y, X, newX, family, k = 10, ...) {
  if (family$family != "binomial") {
    stop("SL.knn only available for family = binomial()")
  }
  # convert all covariates to numerics (knn will not accept factors)
  Xnumeric = as.data.frame(sapply(X, FUN=function(x) as.numeric(x)))
  newXnumeric = as.data.frame(sapply(newX, FUN=function(x) as.numeric(x)))
  # now run knn as usual
  fit.knn <- knn(train = Xnumeric, test = newXnumeric,
                 k = k, cl = Y, prob = TRUE)
  pred <- (as.numeric(fit.knn) - 1) * attr(fit.knn, "prob") +
    (1 - (as.numeric(fit.knn) - 1)) * (1 - attr(fit.knn, "prob"))
  fit <- list(k = k)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.knn")
  return(out)
}
#add my wrapper to SuperLearner namespace
#environment(my.SL.knn) <-asNamespace("SuperLearner")

#same as built-in SL.caret, but run with method=party rather than randomForest
my.SL.caret = function (Y, X, newX, family, obsWeights, method = "ctree", tuneLength = 3, 
                        trControl = trainControl(method = "cv", number = 10, verboseIter = TRUE), 
                        metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) {
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "prob")[,2] #make sure this is giving predicted probabilities (to clear error)
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}
#add my wrapper to SuperLearner namespace
#environment(my.SL.caret) <-asNamespace("SuperLearner")

#run model and obtain CV predicted probabilities
run_CV_and_report = function(output_desc, outD, mydata, pred_method, predictors, predictor_cats, dim_reduce_method, 
                             dim_reduce_num, dim_reduce_covar, cVars, LCMS_vars, clinic_varsD, V, clusternum) {
  #parameters:
    #output_desc: this is a string label to help keep track of output
    #outD: function will return this dataframe with added row of output
    #mydata: data to use in analysis
    #pred_method: options are "logit","RF","SL","ctree". Employed after dimension reduction (if any).
    #predictors: when no dimension reduction is used, these are the predictors used in the regression.
    #predictor_cats: options are "all", "MFs", "non-MFs".  Determines which variables are used in regression.  
    #dim_reduce_method: options are "RF","t-test","lasso".  Determines how subset of predictors should be selected.
    #dim_reduce_num: number of elements in subset selected by dim reduction method.
    #[TODO - dim_reduce_target: options are "all_simult", "all_sep", "MFs", "non-MFs". 
        #Determines on which predictors we should apply subset routine.]
    #dim_reduce_covar: T if we should include non-MF covariates when running dim reduction on MFs. 
      #F if we should ignore covariates when selecting MF subset
      #currently operable for "RF" and "lasso" dimension reduction methods; "t-test" will no consider non-MFs.
    #cVars: clinical variables for consideration (=cvars_noMiss by default)
    #LCMS_vars: MFs for consideration (=allMFs by default)
    #V: number of folds to use in CV
    #clusternum: number of cores (to run multithreading with SL)
  
  cat(paste("\n\n** Running with", predictor_cats, "using", dim_reduce_method, "dimension method",
            "with", dim_reduce_num, "covariates using", pred_method,"prediction method**\n\n"))
  
  #this list will contain variables selected in the subsetting process
  selected_vars = list()
  #this list will contain counts of number of selected variables for the lasso
  dim_reduce_nums = list()
  
  #partition data into 10 equally sized peices -- with stratification according to outcome (10/22/2014 modification)
  set.seed(123)
  #Y0 will be a list of 10 elements, each containing index numbers of a random 1/10 of observations for which Y=0 
  Y0 <- split(sample(which(mydata$is.DEN==F)), rep(1:V, length=length(which(mydata$is.DEN==F)))) 
  #Y1 will be a list of 10 elements, each containing index numbers of a random 1/10 of observations for which Y=1 
  Y1 <- split(sample(which(mydata$is.DEN==T)), rep(1:V, length=length(which(mydata$is.DEN==T))))
  #folds will be a list of 10 elements, each containing 1 of Y0's elements combined with 1 of Y1's elements
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
  
  ### cross-validation step: fit model and predict for left out fold ###
  for(v in seq(from=1, to=V)) {
    
    trainD = mydata[-folds[[v]],] #training set for first run 
    testD = mydata[folds[[v]],]   #test set for first run
    
    ### determine which variables to include in model ###
    if(dim_reduce_method!="none"){
      
      if(dim_reduce_method=="RF"){
        if(dim_reduce_covar==F){
          rf = randomForest(x=trainD[,c(LCMS_vars)], y=as.factor(trainD[,"is.DEN"]), ntree=501, localImp=FALSE)
        }else if(dim_reduce_covar==T){
          rf = randomForest(x=trainD[,c(cVars, LCMS_vars)], y=as.factor(trainD[,"is.DEN"]), ntree=501, localImp=FALSE)
        }
        rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
        rf_VIPS_ordered = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
        rf_VIP_MFs = names(rf_VIPS_ordered[grep("MZ_",names(rf_VIPS_ordered))]) #exclude non-MF variables from list
        MFpreds = rf_VIP_MFs[1:dim_reduce_num] #take top MFs
        
      } else if(dim_reduce_method=="t-test"){
        MFpreds = row.names(ttest_select(trainD, topX=dim_reduce_num, use_clinical=FALSE))
        
      } else if(dim_reduce_method=="lasso"){
        if(dim_reduce_covar==F){
          cvLasso = cv.glmnet(x=as.matrix(trainD[,c(LCMS_vars)]),
                              y=as.matrix(trainD[,"is.DEN"]),
                              family='binomial', alpha=1, nfolds=10, type.measure="deviance")
          cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
          MFpreds = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs (MFs with non-zero coef)
        }else if(dim_reduce_covar==T){
          cvLasso = cv.glmnet(x=as.matrix(trainD[,c(cVars, LCMS_vars)]),
                              y=as.matrix(trainD[,"is.DEN"]),
                              family='binomial', alpha=1, nfolds=10, type.measure="deviance")
          cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
          #put in descending order the non-zero coefficients, excluding intercept
          nonzeros = cvLcoefs[which(cvLcoefs>0),][-1]
          LasSigX = nonzeros[order(nonzeros)]
          MFpreds = names(LasSigX[grep("MZ_",names(LasSigX))]) #exclude the non-MFs
        }
        dim_reduce_nums = append(dim_reduce_nums, length(MFpreds)) #will later take average over all folds
      }
      selected_vars = append(selected_vars, MFpreds) #compile list of selected MFss from all runs
      if(predictor_cats=="all"){
        predictors = c(MFpreds, cVars)
      } else if(predictor_cats=="MFs"){
        predictors = MFpreds
      }
      #print(paste("Number of predictors:", length(predictors)))
    }
    myformula = as.formula(paste("is.DEN ~ ", paste(predictors, collapse= "+")))
    
    ### fit prediction algorithm to training set and predict on test set ###
    if(pred_method=="logit"){
      #logit will return same output whether Y is boolean or a 0/1 numeric factor
      mylogit <- glm(myformula, data=trainD,, family = "binomial")
      #print(summary(mylogit))
      #obtain predicted probabilities for left out chunk using regression fit
      mydata[folds[[v]],"DEN_prob"] <- predict(mylogit, newdata=testD, type="response")
    }
    if(pred_method=="ctree"){
      cfit = ctree(myformula, data=trainD)
      mydata[folds[[v]],"DEN_prob"] = predict(cfit, testD[,c(predictors)], type="response")
    }
    if(pred_method=="RF"){
      #if Y is numeric (even if it is a boolean), RF will run a regression
      #if Y is a factor, RF will do classification, so I therefore convert outcome to factor
        #unlike SL, RF is smart about factors
      rf = randomForest(x=trainD[,c(predictors)], y=as.factor(trainD[,"is.DEN"]), ntree=501)
      mydata[folds[[v]],"DEN_prob"] = predict(rf, testD[,c(predictors)], type="prob")[,2]
    }
    if(pred_method=="SL"){

      #outcome variable must be a numeric vector (SL does not accept other types)
      #knn (for super learner) also has problems with X vars that are factors with 2 levels so must convert everything to numeric
      #family="gaussian" in SL means when SL calls RandomForest, it will first convert Y into a factor so so RF will do classification
      #YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
      #                           "Variable.Name.in.Hospital.data"]
      #YesNo_varsH_char = sapply(YesNo_varsH, as.character)
      #trainD_num = convert_factor_to_numeric(trainD, c("is.DEN","Torniquete","Sexo",YesNo_varsH_char))
      #testD_num = convert_factor_to_numeric(testD, c("is.DEN","Torniquete","Sexo",YesNo_varsH_char))
      
      #Super Learner
      #note: not all algorithms are estimable given that I have more predictors than observations 
      #gam, glm will not run when all clinical varaibles are included
      #glmnet produces warnings when using just the clinical data, though seems fine with all data
        #also, it will crash if fed only 1 predictor
      #default for SL.glmnet is to use lambda that gives minimum cvm
      #caret is a general template that can run a variety of methods, including party
      #
      if(length(predictors)>1){
        #Error in trimLogit(Z, trim = control$trimLogit) : argument "control" is missing, with no default
        mySL.library = c("SL.mean","SL.randomForest","SL.glmnet") #"my.SL.caret", "my.SL.knn", ,"SL.polymars"
      }else{
        #I might want to add glm for when predictor set is small, but won't do it now since I don't want to 
          #present results for which prediction error decreases with fewer predictors
        mySL.library = c("SL.mean", "SL.randomForest")  #removing knn for now (11-28-2014)
      }
      cl <- makeCluster(clusternum, type = "PSOCK") # can use different types here
      clusterSetRNGStream(cl, iseed = 2343)
      clusterExport(cl, c("my.SL.knn","my.SL.caret")) #copy my created SL functions to all clusters
      SuperLearner.CV.control = list(10L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
      #SLresults = SuperLearner(Y=as.numeric(trainD[, "is.DEN"]), 
      #                             X=as.data.frame(trainD[,c(predictors)]), newX=as.data.frame(testD[,c(predictors)]),
      #                             family=binomial(), SL.library=mySL.library, cvControl=SuperLearner.CV.control,
      #                             method = "method.NNloglik", verbose=FALSE) #9-28-14: changed method from NNLS to NNloglik (Alan's suggestion)
      SLresults = snowSuperLearner(cluster=cl, Y=as.numeric(trainD[, "is.DEN"]), 
                               X=as.data.frame(trainD[,c(predictors)]), newX=as.data.frame(testD[,c(predictors)]),
                               family=binomial(), SL.library=mySL.library, cvControl=SuperLearner.CV.control,
                               method = "method.NNloglik", verbose=FALSE) #9-28-14: changed method from NNLS to NNloglik (Alan's suggestion)
      mydata[folds[[v]],"DEN_prob"] = SLresults$SL.predict #predicts for obs in newX, as we want
      if(v==1){
        SLcoefs = as.data.frame(t(c(v,SLresults$coef)))
        #print(SLcoefs)
      }else{
        newrow = as.data.frame(t(c(v,SLresults$coef)))
        #print(newrow)
        SLcoefs = smartbind(SLcoefs, newrow, fill=NA)
      }
      stopCluster(cl)
    }
    
  } #end of CV step 
  
  #present coefficients for algorithms chosen by Super Learner
  if(pred_method=="SL"){
    print("Coefficients for algorithms chosen by Super Learner (for each CV loop)")
    print(SLcoefs)
  }
  
  #compile list of selected predictors from across cross-validation fits
  if(dim_reduce_method!="none"){
    MF_freqs = subset(as.data.frame(table(t(subset(as.data.frame(table(selected_vars)), select=-c(Freq)))[,1])))
    MF_freqs_ordered = MF_freqs[order(MF_freqs$Freq, decreasing=TRUE),]
    print("Frequencies at which MFs were chosen (across CV loop)")
    print(MF_freqs_ordered)
    if(dim_reduce_method=="lasso"){
      print("Number of chosen MFs during lasso selection (across CV loop)")
      print(as.data.frame(dim_reduce_nums))
      dim_reduce_num = rowMeans(as.data.frame(dim_reduce_nums))
    }
  } else {
    MF_freqs_ordered = list() #create this empty list so there is no complaining
  }
  
  cvAUCreport <- ci.cvAUC(predictions=mydata$DEN_prob, labels=mydata$is.DEN, folds=folds, confidence=0.95)
  print(cvAUCreport)
  myreport = report_performance(mydata[,c("is.DEN","DEN_prob")])
  print(myreport)
  output = list(output_desc, predictor_cats, dim_reduce_covar, dim_reduce_method, dim_reduce_num, pred_method, 
                myreport$sensitivity, myreport$specificity, myreport$pred_err, 
                myreport$ROC_area, myreport$MSE, myreport$SD_SqEr, 
                cvAUCreport$cvAUC, cvAUCreport$se, cvAUCreport$ci[[1]], cvAUCreport$ci[[2]], cvAUCreport$confidence)
  outD = rbind(outD, output)
  #return(outD)
  return(list(outD, MF_freqs_ordered)) 
}

#merge lists of selected MFs
merge_lists = function (list1, name1, list2, name2, nameC){
  colnames(list1)[2] = name1
  colnames(list2)[2] = name2
  m = merge(list1, list2, by="Var1", all=TRUE)
  FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
  m[,name1]= sapply(m[,name1], FUN=FillMissingsWithZeros)
  m[,name2]= sapply(m[,name2], FUN=FillMissingsWithZeros)
  m[,nameC] = m[,name1] + m[,name2]
  selected_MFs = m[order(m[,nameC], decreasing=TRUE),]
  print("Significant MFs across runs")
  print(selected_MFs)
  return(selected_MFs)
}


#function takes preliminary set of clinical predictors and cleans it up according to fixed and user-supplied criteria
get_clinical_predictors = function(compare_grps, clinic_varsD, covarlist, XD, include_imp_dums){
  
  if(compare_grps=="OFI.v.DEN" | compare_grps=="OFIDF.v.DHFDSS"){
    #obtain list of candidate clinical variables to include in ND vs. DEN prediction
    #drop the variables in covarlist that are not applicable to ND vs DEN prediction (e.g., IR and PCR)
    clinic_vars = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% covarlist==T 
                                     & clinic_varsD$Use.in.ND.vs.DEN.prediction==1),"variable.name.in.final.data"]
  }
  if(compare_grps=="DF.v.DHFDSS"){
    #obtain list of candidate clinical variables to include in DF vs. DHF/DSS prediction
    #drop the variables in covarlist that are not applicable to DF vs. DHF/DSS prediction (e.g., vars that define DHF/DSS)
      #Note: we have reconsidered analysis and are now using even the vars that are part of the DHF/DSS definition to predict DHF/DSS
        #however, we are now (in python) dropping patients who already qualify as DHF/DSS in first consultation
    clinic_vars = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% covarlist==T 
                                     & clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),"variable.name.in.final.data"]
  }
  
  #obtain list of clinical variables that have no missing values in data
  clinic_vars_missCount = sapply(XD[,c(as.character(clinic_vars))], function(x) sum(is.na(x)))
  cvars_noMiss = names(which(clinic_vars_missCount==0))
  cat("\n clinical variables that shall be included in regressions \n")
  print(cvars_noMiss) #view 
  cat("\n variables excluded due to missing values \n")
  print(names(which(clinic_vars_missCount>0))) #view 
  
  #add the dummy variables indicating imputed values to list of clinical variables 
  if(include_imp_dums=="all" | include_imp_dums=="Only"){
    allvars = colnames(XD)
    imp_found = allvars[grepl("_imputed", allvars)] #all imputed indicators in data
    imp_imagine = paste(cvars_noMiss, "_imputed", sep="") 
    imp_matched = imp_found[imp_found %in% imp_imagine==TRUE] 
    cat("\n preliminary list of imputation indicators to include regressions \n")
    print(imp_matched)
    XD$is.cohort = (XD$Study=="Cohort")
    is.diff_fun = function(var, XD) (sum(XD$is.cohort!=XD[,var]))
    diffs = sapply(imp_matched, FUN=is.diff_fun, XD)
    #get list of vars that do not perfectly correspond to is.cohort
    unique_imp_matched = names(diffs[diffs!=0])
    #add in one of the is.cohort variables if there is one
    if(length(names(diffs[diffs==0]))>0) {
      unique_imp_matched = c(unique_imp_matched, names(diffs[diffs==0])[1])
    }
    cat("\n imputation indicators to include in regressions (excludes redundancies) \n")
    print(unique_imp_matched)
    if (include_imp_dums=="Only"){
      cvars_noMiss = unique_imp_matched #include only missing indicators (experiment)
    }
    else{
      cvars_noMiss = c(cvars_noMiss, unique_imp_matched)
    }
  }
  #add study (cohort v hospital) indicator
  if(include_imp_dums=="study"){
    cat("\n adding indicator of study type (cohort v hospital) to clinical var list \n")
    cvars_noMiss = c(cvars_noMiss, "is.cohort")
  }
  return(cvars_noMiss)
}

#Run a variety of prediction algorithms and output results
run_predictions = function(clinic_varsD, covarlist, compare_grps, XD, output_desc, include_imp_dums="all"){
  #clinic_varsD is the data containing information on each clinical/lab/demographic variable in data
  #covarlist is the list of clinical/lab/demogrpahic info to include in the analysis
    #only manipulation to this list will be to eliminate vars inappropriate for predicting outcome
    #and also to eliminate vars that have missing values
  #compare_grps: what we want to predict -- either OFI.v.DEN" or "DF.v.DHFDSS" or "OFIDF.v.DHFDSS"
  #XD is the data from which we should find covariate and outcome values
  #output_desc is a string used in output to describe the prediction run
  #include_imp_dums: 
      #if "all" then prediction will be run with all clinical variables found in the data that end with "imputed" and 
         #that match with a clinical variable selected for use in the prediction algorithms (these are indicators of missingness in data)
      #if "study" then include a dummy indicating hospital study, but no other imputation indicators
      #if "Only" then prediction will be run with only the imputation indicators (and not the actual clinical variable values)
      #otherwise, do not include
  
  #for testing purposes: 
  #compare_grps = "OFI.v.DEN"
  #XD = clinical_full_clean
 
  cat("\n\n",paste("***************************************************************************"),"\n")
  cat(paste("************ Running",compare_grps,"analysis with", output_desc, " *************"),"\n")
  cat(paste("***************************************************************************"),"\n")
  
  ###############################################################################
  ############ develop formulas (indicating covariates) to try ##################
  
  #todo: to get more honest inference, should I imput missing values within the CV step? 
  
  if(compare_grps=="DF.v.DHFDSS"){
    #drop OFI observations 
    XD = XD[which(XD$is.DEN==TRUE),]
    #for now, rename DF vs. DHF/DSS outcome to the ND vs. DEN outcome (avoids having to change code)
    XD$is.DEN = XD$is.DHF_DSS
  }
  #add study (cohort v hospital) indicator
  if(include_imp_dums=="study"){
    XD$is.cohort = (XD$Study=="Cohort")
  }
  
  #obtain list of clinical predictors
  cvars_noMiss = get_clinical_predictors(compare_grps, clinic_varsD, covarlist, XD, include_imp_dums)
  clinicCount = length(cvars_noMiss)
  cat(paste("Total number of clinical variables to include",clinicCount),"\n")
  
  #obtain list of MFs 
  cindices = grep("MZ_",colnames(XD)) #column indices in dataset
  allMFs = colnames(XD[cindices]) #names
  MFcount = length(allMFs)
  
  #count of predictors when using MFs and clinical vars
  allCount = clinicCount + MFcount

  #initialize dataframe to hold results
  outD=data.frame(x0=character(0),x1=character(0), x2=logical(0), x3=character(0), x4=integer(0), x5=character(0),
                  x6=integer(0), x7=integer(0), x8=integer(0), x9=integer(0), x10=integer(0), 
                  x11=integer(0), x12=integer(0), x13=integer(0), x14=integer(0), x15=integer(0), 
                  x16=integer(0), stringsAsFactors = FALSE)
  colnames(outD) = c("output_desc","predictors", "dim_reduce_covar", "reduce_method", "r_num", "pred_method", 
                     "sensitivity", "specificity", "pred_err","ROC_area","MSE", "SD_SqEr",
                     "cvAUC","cvAUC_se", "cvAUC_ci_lower","cvAUC_ci_upper", "cvAUC_confidence")
  blanks = list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) 
  outD[1,] = blanks #seems ridiculous, but adding this blank row is only way I could get it working
  #this will hold info on which MFs were selected, if such an analysis is desired
  selected_MFs = "NA"
  
  ## Rashika and Natalia supplied the following list for me to use in prediction (3-12-2014 email)
     # (I did a fuzzy match to obtain exact mz values found in Nicaragua serum runs of LC-MS)
  if(output_desc=="serum, D1"){
    selectMFs = c("MZ_522.3552", "MZ_496.34033", "MZ_496.34","MZ_782.5712","MZ_784.5826",
                "MZ_772.549","MZ_279.231","MZ_329.2458","MZ_209.15569",
                "MZ_369.35104","MZ_369.35147","MZ_417.33545","MZ_417.33557","MZ_116.07061")
  }else if(output_desc=="serum, D2"){
    selectMFs = c("MZ_522.355", "MZ_496.3399", "MZ_209.15659", "MZ_369.3519", "MZ_417.3369")
  }else{
    selectMFs = c()
  }
  selectCountMF = length(selectMFs)
  selectCountAll = selectCountMF + clinicCount
  
  #avoid having to pass in certain variables each time (mydata, cVars, LCMS_vars, clinic_varsD)
  run_CV_and_report_bound = function(
    outD, pred_method, predictors, predictor_cats, dim_reduce_method, 
    dim_reduce_num, dim_reduce_covar) {
    return(run_CV_and_report(
      output_desc=output_desc, outD, mydata=XD, pred_method, predictors, predictor_cats, dim_reduce_method, 
      dim_reduce_num, dim_reduce_covar,
      cVars=cvars_noMiss, LCMS_vars=allMFs, clinic_varsD=clinic_varsD, V=10, clusternum=16))
  }
  
  ###############################################################################
  ######## Obtain predictions using 10-fold CV and report performance ###########  
  
  outD = run_CV_and_report_bound(outD, "SL", cvars_noMiss, "non-MFs", "none", clinicCount, T)[[1]]
  #outD = run_CV_and_report_bound(outD, "RF", c(allMFs, cvars_noMiss), "all", "none", allCount, T)[[1]]
  #outD = run_CV_and_report_bound(outD, "ctree", cvars_noMiss, "non-MFs", "none", clinicCount, T)[[1]]
  
  #### Random Forests Runs ####
  if(F==T){
 
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "none", MFcount, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", c(allMFs, cvars_noMiss), "all", "none", allCount, T)[[1]]
 
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "t-test", 5, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "t-test", 4, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "t-test", 3, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "t-test", 2, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "t-test", 1, F)[[1]]
  #outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "lasso", 99, F)[[1]]
  
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "RF", 5, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "RF", 4, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "RF", 3, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "RF", 2, F)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "RF", 1, F)[[1]] #works!
  
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "t-test", 5, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "t-test", 4, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "t-test", 3, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "t-test", 2, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "t-test", 1, T)[[1]]
  #outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "lasso", 99, F)[[1]] #for lasso, do not use covars (might not choose any MFs)
  
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 5, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 4, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 3, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 2, T)[[1]]
  outD = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 1, T)[[1]]
}
  
  #### Super Learner Runs ####
if(F==T){
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "none", MFcount, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", cvars_noMiss, "non-MFs", "none", clinicCount, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", c(allMFs, cvars_noMiss), "all", "none", allCount, T)[[1]]
  
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 5, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 4, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 3, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 2, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 1, F)[[1]] #gives error
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "lasso", 99, F)[[1]]

  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 5, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 4, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 3, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 2, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 1, F)[[1]] #gives error
  
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 5, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 4, F)[[1]] #trimLogit err, DEN, all clinical, D1
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 3, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 2, F)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 1, F)[[1]]
  #outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "lasso", 99, F)[[1]]
  
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 5, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 4, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 3, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 2, T)[[1]]
  outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 1, T)[[1]]
  }
  ## Obtain results for SL run with MFs selected by Rashika and Natalia
  #outD = run_CV_and_report_bound(outD, "SL", selectMFs, "MFs", "none", selectCountMF,T)[[1]]
  #outD = run_CV_and_report_bound(outD, "SL", c(selectMFs, cvars_noMiss), "all", "none",selectCountAll,T)[[1]]

  #### Obtain predictors chosen by random forest over the course of the CV step ####
  #outA = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 10, T)
  #outD = outA[[1]]
  #selectMFs_RF10 = outA[[2]] 
  #outA = run_CV_and_report_bound(outD, "RF", allMFs, "all", "RF", 5, T)
  #outD = outA[[1]]
  #selectMFs_RF5 = outA[[2]]
  #selected_MFs = merge_lists(list1=selectMFs_RF10, name1="RF10",list2=selectMFs_RF5,  name2="RF5", nameC="total")
  #selected_MFs

  main_out = outD[-1,] #drop the first row (it's blank)
  #main_out = outD
  #return(main_out) 
  return(list(main_out,selected_MFs)) 
}

  