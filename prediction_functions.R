
#create a version of the data which replaces values for missing compounds with zeros
#other options: fill with min value observed for that compound,
#fill with 1/2 of the min value, input value
process_data = function(mydata){
  cindices = grep("MZ_",colnames(mydata))
  processedD = mydata
  #note: must have a return statement b/c assignment operator will not return anything (and default is 0)
  FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
  processedD[cindices] = lapply(mydata[cindices],FUN=FillMissingsWithZeros)
  return(processedD)
}


#obtain features that demonstrate at least a 2-fold differentiation based on t-test with a p-value<.05.  
#note: t-test default is to set var.equal=FALSE, which means using the Welch test (good)
get_ttest_results <- function(mydata){
  myttest <- function(x){
    ttest = t.test(x[which(mydata$DEN_dum==1)],
                   x[which(mydata$DEN_dum==0)])
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
    mydata$DEN_dumP = NA
    mydata$DEN_dumP[mydata$DEN_prob<threshold] = 0 
    mydata$DEN_dumP[mydata$DEN_prob>=threshold] = 1 
    mydata[,c("DEN_dum","DEN_prob","DEN_dumP")] #view actual versus predicted
    #Sensitivity (also called the true positive rate) measures the proportion of actual positives 
    #which are correctly identified as such (e.g. the percentage of sick people who are correctly identified as having the condition). 
    sensitivity = round(sum(mydata$DEN_dum==1 & mydata$DEN_dumP==1)/sum(mydata$DEN_dum==1),3)
    sensitivity
    #Specificity measures the proportion of negatives which are correctly identified as such 
    #(e.g. the percentage of healthy people who are correctly identified as not having the condition, 
    #sometimes called the true negative rate).
    specificity = round(sum(mydata$DEN_dum==0 & mydata$DEN_dumP==0)/sum(mydata$DEN_dum==0),3)
    specificity
    #Overall error rate
    pred_err = round(sum(mydata$DEN_dum != mydata$DEN_dumP)/sum(!is.na(mydata$DEN_dum)),3)
    pred_err
    #ROC curves plot the false positive rate (1-specificity) against the true positive rate (sensitivity)
    return(c(sensitivity, specificity, pred_err))
  }
  performance = data.frame(t(sapply(X=seq(1, 0, -.05), FUN=measure_performance, mydata=pdata)))
  colnames(performance) = c("sensitivity", "specificity", "pred_err")
  #ROC plot (uncomment if desired)
  #plot(x=(1-performance$specificity), y=performance$sensitivity, 
  #     col="grey", lwd=2, type="l",
  #     ylab="true positive rate", xlab="false positive rate",main="ROC curve with 10-fold cross validation")
  #specificity, sensitivity, and error rates when threshold is chosen so as to minimize error rate
  report = performance[which(performance$pred_err==min(performance$pred_err)),][1,] #chose the first when there is a tie
  #calculate area under the ROC curve (does not work ever since installing the glmnet package --- asked MESS package auc function)
  report$ROC_area = round(myauc(x=(1-performance$specificity), y=performance$sensitivity),3)
  #calculate mean squared error (using actual outcome and predicted probability)
  report$MSE = round(mean( (as.vector(pdata$DEN_prob) - as.numeric(as.character(pdata$DEN_dum)))^2 ),3)
  report$SD_SqEr = round(sd( (as.vector(pdata$DEN_prob) - as.numeric(as.character(pdata$DEN_dum)))^2 ),3)
  return(report)
}


#run model and obtain CV predicted probabilities
run_CV_and_report = function(outD, mydata, pred_method, predictors, predictor_cats, dim_reduce_method, 
                             dim_reduce_num, dim_reduce_covar, reduce_in_CV, cVars, LCMS_vars, clinic_varsD, V) {
  #parameters:
    #outD: function will return this dataframe with added row of output
    #mydata: data to use in analysis
    #pred_method: options are "logit","RF","SL". Employed after dimension reduction (if any).
    #predictors: when no dimension reduction is used, these are the predictors used in the regression.
    #predictor_cats: options are "all", "MFs", "non-MFs".  Determines which variables are used in regression.  
    #dim_reduce_method: options are "RF","t-test","lasso".  Determines how subset of predictors should be selected.
    #dim_reduce_num: number of elements in subset selected by dim reduction method.
    #[TODO - dim_reduce_target: options are "all_simult", "all_sep", "MFs", "non-MFs". Determines on which predictors we should apply subset routine.]
    #dim_reduce_covar: T if we should include non-MF covariates when running dim reduction on MFs. 
      #F if we should ignore covariates when selecting MF subset
      #currently only operable for "RF" and "lasso" dimension reduction methods
    #reduce_in_CV: T if dimension reduction should be done in CV step.  F is dim reduc should be done before CV step.
    #cVars: clinical variables for consideration (=cvars_noMiss by default)
    #LCMS_vars: MFs for consideration (=allMFs by default)
    #V: number of folds to use in CV
  cat(paste("\n\n**** running for", predictor_cats, dim_reduce_method, dim_reduce_num, pred_method,"****\n\n"))
  
  #this list will contain variables selected in the subsetting process
  selected_vars = list()
  #this list will contain counts of number of selected variables for the lasso
  dim_reduce_nums = list()
  
  #partition data into 10 equally sized peices
  #set.seed(123)
  #nobs = nrow(mydata)
  #partition_nobs <- ceiling(0.1 * nobs) #size of each peice
  #nvec = rep(seq(from=1, to=10),partition_nobs)[1:nobs] #vector of values 1 - 10 of length nobs with [approximately] equal numbers of each value
  #rvec = permute(nvec) #random permutation of vector
  
  #partition data into 10 equally sized peices -- with stratification according to outcome (10/22/2014 modification)
  set.seed(123)
  #Y0 will be a list of 10 elements, each containing index numbers of a random 1/10 of observations for which Y=0 
  Y0 <- split(sample(which(mydata$DEN_dum==0)), rep(1:V, length=length(which(mydata$DEN_dum==0)))) 
  #Y1 will be a list of 10 elements, each containing index numbers of a random 1/10 of observations for which Y=1 
  Y1 <- split(sample(which(mydata$DEN_dum==1)), rep(1:V, length=length(which(mydata$DEN_dum==1))))
  #folds will be a list of 10 elements, each containing 1 of Y0's elements combined with 1 of Y1's elements
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
  
  ### cross-validation step: fit model and predict for left out fold ###
  for(v in seq(from=1, to=V)) {
    
    trainD = mydata[-folds[[v]],] #training set for first run 
    testD = mydata[folds[[v]],]   #test set for first run
    
    ### determine which variables to include in model ###
    if(reduce_in_CV==T & dim_reduce_method!="none"){
      if(dim_reduce_method=="RF"){
        if(dim_reduce_covar==F){
          rf = randomForest(x=trainD[,c(LCMS_vars)], y=trainD[,"DEN_dum"], ntree=500, localImp=FALSE)
        }else if(dim_reduce_covar==T){
          rf = randomForest(x=trainD[,c(cVars, LCMS_vars)], y=trainD[,"DEN_dum"], ntree=500, localImp=FALSE)
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
                              y=as.matrix(trainD[,"DEN_dum"]),
                              family='binomial', alpha=1, nfolds=10, type.measure="deviance")
          cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
          MFpreds = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs (MFs with non-zero coef)
        }else if(dim_reduce_covar==T){
          cvLasso = cv.glmnet(x=as.matrix(trainD[,c(cVars, LCMS_vars)]),
                              y=as.matrix(trainD[,"DEN_dum"]),
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
    myformula = as.formula(paste("DEN_dum ~ ", paste(predictors, collapse= "+")))
    #print(myformula)
    
    ### fit prediction algorithm to training set and predict on test set ###
    if(pred_method=="logit"){
      mylogit <- glm(myformula, data = trainD, family = "binomial")
      #print(summary(mylogit))
      #obtain predicted probabilities for left out chunk using regression fit
      mydata[folds[[v]],"DEN_prob"] <- predict(mylogit, newdata=testD, type="response")
    }
    if(pred_method=="RF"){
      rf = randomForest(myformula, data=trainD, ntree=500)
      mydata[folds[[v]],"DEN_prob"] = predict(rf, testD, type="prob")[,2]
    }
    if(pred_method=="SL"){

      #outcome variable must be a numeric vector
      #knn (for super learner) also has problems with X vars that are factors with 2 levels 
      #so must convert everything to numeric
      YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                                 "Variable.Name.in.Hospital.data"]
      YesNo_varsH_char = sapply(YesNo_varsH, as.character)
      trainD_num = convert_factor_to_numeric(trainD, c("DEN_dum","Torniquete","Sexo",YesNo_varsH_char))
      testD_num = convert_factor_to_numeric(testD, c("DEN_dum","Torniquete","Sexo",YesNo_varsH_char))
      
      #Super Learner
      #note: not all algorithms are estimable given that I have more predictors than observations 
      #gam, glm will not run when all clinical varaibles are included
      #glmnet produces warnings when using just the clinical data, though seems fine with all data
        #also, it will crash if fed only 1 predictor
      #default for SL.glmnet is to use lambda that gives minimum cvm
      if(length(predictors)>1){
        mySL.library = c("SL.mean", "SL.knn", "SL.randomForest", "SL.glmnet")   
      }else{
        #I might want to add glm for when predictor set is small, but won't do it now since I don't want to 
          #present results for which prediction error decreases with fewer predictors
        mySL.library = c("SL.mean", "SL.knn", "SL.randomForest")
      }
      #print(colnames(as.data.frame(trainD_num[,c(predictors)])))
      #print(colnames(as.data.frame(testD_umn[,c(predictors)])))
      #print(mySL.library)
      SLresults = SuperLearner(Y=trainD_num[, "DEN_dum"], 
                               X=as.data.frame(trainD_num[,c(predictors)]), newX=as.data.frame(testD_num[,c(predictors)]),
                               family=binomial(), SL.library=mySL.library,
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
    }
    
  } #end of CV step 
  
  #present coefficients for algorithms chosen by Super Learner
  if(pred_method=="SL"){
    print("Coefficients for algorithms chosen by Super Learner (for each CV loop)")
    print(SLcoefs)
  }
  
  #compile list of selected predictors from across cross-validation fits
  if(reduce_in_CV==T & dim_reduce_method!="none"){
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
  
  cvAUCreport <- ci.cvAUC(predictions=mydata$DEN_prob, labels=mydata$DEN_dum, folds=folds, confidence=0.95)
  print(cvAUCreport)
  myreport = report_performance(mydata[,c("DEN_dum","DEN_prob")])
  print(myreport)
  output = list(predictor_cats, dim_reduce_covar, dim_reduce_method, dim_reduce_num, pred_method, 
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


#Run a variety of prediction algorithms and output results
run_predictions = function(clinic_varsD, outcome, mydata, sample_name, reduce_in_CV=T){
  
  #for testing purposes: 
  #outcome = "ND.vs.DEN"
  #mydata = comboD1_filter50n
 
  cat("\n\n",paste("************ Running",outcome,"analysis with", sample_name, " *************"),"\n")
  
  #create a version of the data which replaces values for missing compounds with zeros
  XD = process_data(mydata)
  
  ###############################################################################
  ############ develop formulas (indicating covariates) to try ##################
  
  if(outcome=="ND.vs.DEN"){
    #obtain list of candidate clinical variables to include in ND vs. DEN prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                              "Variable.Name.in.Hospital.data"]),"age","DaysSick") #43+2
  }
  if(outcome=="DF.vs.DHF.DSS"){
    #obtain list of candidate clinical variables to include in DF vs. DHF/DSS prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),
                                              "Variable.Name.in.Hospital.data"]),"age","DaysSick") #30
    #drop ND observations 
    XD = XD[!is.na(XD$DHF_DSS_dum),]
    #for now, rename DF vs. DHF/DSS outcome to the ND vs. DEN outcome (avoids having to change code)
    XD$DEN_dum = XD$DHF_DSS_dum
  }
  
  #obtain list of MFs 
  cindices = grep("MZ_",colnames(XD)) #column indices in dataset
  allMFs = colnames(XD[cindices]) #names
  MFcount = length(allMFs)
  
  #obtain list of clinical variables that have no missing values in data
  clinic_vars_missCount = sapply(XD[clinic_vars], function(x) sum(is.na(x)))
  cvars_noMiss = clinic_vars[clinic_vars_missCount==0] #40 
  clinicCount = length(cvars_noMiss)
  cat("\n clinical variables that shall be included in regressions \n")
  print(cvars_noMiss) #view 
  cat("\n variables excluded due to missing values \n")
  print(clinic_vars[clinic_vars_missCount>0]) #view 
  allCount = clinicCount + MFcount
  
  #initialize dataframe to hold results
  outD=data.frame(x1=character(0), x2=logical(0), x3=character(0), x4=integer(0), x5=character(0),
                  x6=integer(0), x7=integer(0), x8=integer(0), x9=integer(0), x10=integer(0), 
                  x11=integer(0), x12=integer(0), x13=integer(0), x14=integer(0), x15=integer(0), 
                  x16=integer(0), stringsAsFactors = FALSE)
  colnames(outD) = c("predictors", "dim_reduce_covar", "reduce_method", "r_num", "pred_method", 
                     "sensitivity", "specificity", "pred_err","ROC_area","MSE", "SD_SqEr",
                     "cvAUC","cvAUC_se", "cvAUC_ci_lower","cvAUC_ci_upper", "cvAUC_confidence")
  blanks = list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) 
  outD[1,] = blanks #seems ridiculous, but adding this blank row is only way I could get it working
  #this will hold info on which MFs were selected, if such an analysis is desired
  selected_MFs = "NA"
  
  ## Rashika and Natalia supplied the following list for me to use in prediction (3-12-2014 email)
     # (I did a fuzzy match to obtain exact mz values found in Nicaragua serum runs of LC-MS)
  if(sample_name=="serum, D1"){
    selectMFs = c("MZ_522.3552", "MZ_496.34033", "MZ_496.34","MZ_782.5712","MZ_784.5826",
                "MZ_772.549","MZ_279.231","MZ_329.2458","MZ_209.15569",
                "MZ_369.35104","MZ_369.35147","MZ_417.33545","MZ_417.33557","MZ_116.07061")
  }else if(sample_name=="serum, D2"){
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
      outD, mydata=XD, pred_method, predictors, predictor_cats, dim_reduce_method, 
      dim_reduce_num, dim_reduce_covar, reduce_in_CV=reduce_in_CV, 
      cVars=cvars_noMiss, LCMS_vars=allMFs, clinic_varsD=clinic_varsD, V=10))
  }
  
  ###############################################################################
  ######## Obtain predictions using 10-fold CV and report performance ###########  
  
  
  if(reduce_in_CV==T){
    
    #### Random Forests Runs ####
    
    outD = run_CV_and_report_bound(outD, "RF", cvars_noMiss, "non-MFs", "none",clinicCount, T)[[1]]
    outD = run_CV_and_report_bound(outD, "RF", allMFs, "MFs", "none", MFcount, T)[[1]]
    outD = run_CV_and_report_bound(outD, "RF", c(allMFs, cvars_noMiss), "all", "none", allCount, T)[[1]]
    if(F==T){
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
    
    outD = run_CV_and_report_bound(outD, "SL", cvars_noMiss, "non-MFs", "none",clinicCount, T)[[1]]
    if(F==T){
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "none", MFcount, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", c(allMFs, cvars_noMiss), "all", "none", allCount, T)[[1]]
    
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 5, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 4, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 3, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 2, F)[[1]]
    #outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "t-test", 1, F)[[1]] #gives error
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "lasso", 99, F)[[1]]

    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 5, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 4, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 3, F)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 2, F)[[1]]
    #outD = run_CV_and_report_bound(outD, "SL", allMFs, "MFs", "RF", 1, F)[[1]] #gives error
    
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 5, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 4, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 3, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 2, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "t-test", 1, T)[[1]]
    outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "lasso", 99, F)[[1]]
    
    #outD = run_CV_and_report_bound(outD, "SL", allMFs, "all", "RF", 5, F)[[1]]
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
    
  }
  
  ###############################################################################
  ######## Obtain predictions using 10-fold CV and report performance ###########
  
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
                   XD[,"DEN_dum"])
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
    rf = randomForest(x=XD[,c(allMFs)], y=XD[,"DEN_dum"], ntree=500, localImp=FALSE)
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
                        y=as.matrix(XD[,"DEN_dum"]),
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

  main_out = outD[-1,] #drop the first row (it's blank)
  #main_out = outD
  #return(main_out) 
  return(list(main_out,selected_MFs)) 
}


#Get variable importance indicators from Random Forest (with no CV)
get_VIM_RF = function(clinic_varsD, outcome, mydata, dim_reduce_covar, mycolname){
  set.seed(133)
  
  #for testing purposes: 
  #outcome = "ND.vs.DEN"
  #mydata = comboD3_filter50n
  #dim_reduce_covar = T
  
  #create a version of the data which replaces values for missing compounds with zeros
  XD = process_data(mydata)
  
  if(outcome=="ND.vs.DEN"){
    #obtain list of candidate clinical variables to include in ND vs. DEN prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                              "Variable.Name.in.Hospital.data"]),"age","DaysSick") #43
  }
  if(outcome=="DF.vs.DHF.DSS"){
    #obtain list of candidate clinical variables to include in DF vs. DHF/DSS prediction
    clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.DF.vs.DHF.DSS.prediction==1),
                                              "Variable.Name.in.Hospital.data"]),"age","DaysSick") #30
    #drop ND observations 
    XD = XD[!is.na(XD$DHF_DSS_dum),]
    #for now, rename DF vs. DHF/DSS outcome to the ND vs. DEN outcome (avoids having to change code)
    XD$DEN_dum = XD$DHF_DSS_dum
  }
  
  #obtain list of MFs 
  cindices = grep("MZ_",colnames(XD)) #column indices in dataset
  allMFs = colnames(XD[cindices]) #names
  
  if(dim_reduce_covar==F){
    rf = randomForest(x=XD[,c(allMFs)], y=XD[,"DEN_dum"], ntree=500, localImp=FALSE)
  }else if(dim_reduce_covar==T){
    #obtain list of clinical variables that have no missing values in data
    clinic_vars_missCount = sapply(XD[clinic_vars], function(x) sum(is.na(x)))
    cvars_noMiss = clinic_vars[clinic_vars_missCount==0]  
    rf = randomForest(x=XD[,c(cvars_noMiss, allMFs)], y=XD[,"DEN_dum"], ntree=500, localImp=FALSE)
  }
  rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
  rf_VIPS_ordered = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
  rf_VIP_MFs = rf_VIPS_ordered[grep("MZ_",names(rf_VIPS_ordered))] #exclude non-MF variables from list
  rf_VIP_MFs_formated = as.data.frame(t(t(rf_VIP_MFs)))
  colnames(rf_VIP_MFs_formated)=mycolname
  
  return(rf_VIP_MFs_formated)
}
  
  