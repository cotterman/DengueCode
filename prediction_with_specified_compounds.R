


#partition data into 10 equally sized peices
just_run_SL = function(mydata, MFpreds, predictor_cats){
  
  set.seed(123)
  nobs = nrow(mydata)
  partition_nobs <- ceiling(0.1 * nobs) #size of each peice
  nvec = rep(seq(from=1, to=10),partition_nobs)[1:nobs] #vector of values 1 - 10 of length nobs with [approximately] equal numbers of each value
  rvec = permute(nvec) #random permutation of vector

  for(i in seq(from=1, to=10)) {
    trainD = mydata[which(rvec!=i),] #training set for first run 
    testD = mydata[which(rvec==i),]  #test set for first run
    
      if(predictor_cats=="all"){
        predictors = c(MFpreds, cVars)
      } else if(predictor_cats=="MFs"){
        predictors = MFpreds
      }

      mySL.library = c("SL.mean", "SL.knn", "SL.randomForest", "SL.glmnet")   
      SLresults = SuperLearner(Y=trainD[, "DEN_dum"], 
                               X=trainD[,c(predictors)], newX=testD[,c(predictors)],
                               family=binomial(), SL.library=mySL.library,
                               method = "method.NNLS", verbose=FALSE)
      mydata[which(rvec==i),"DEN_prob"] = SLresults$SL.predict #predicts for obs in newX, as we want
    
  }
  
  myreport = report_performance(mydata[,c("DEN_dum","DEN_prob")])
  print(myreport)
  output = list(predictor_cats, dim_reduce_method, dim_reduce_num, pred_method, 
                myreport$sensitivity, myreport$specificity, myreport$pred_err, myreport$ROC_area)
  return(output)
}

processedD1 = process_data(comboD1_filter50n)

just_run_SL(mydata=processedD1, MFpreds="MFs", predictor_cats="all")

522.355
496.34
782.571
784.583
772.549
279.231
329.246
209.156
369.35
417.335
116.071


