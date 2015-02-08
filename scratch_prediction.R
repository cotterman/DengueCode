
###############################################################################
###################### Experiment with prediction algorithms ##################
###############################################################################


################## impute values using SL ########################


#obtain list of clinical variables that have no missing values in data
clinic_vars_missCount = sapply(clinical_full_clean[varlist], function(x) sum(is.na(x)))
cvars_noMiss = varlist[clinic_vars_missCount==0] 
cvars_noMiss
#get list of vars with missing values in data
to_impute_missCount = sapply(clinical_full_clean[,varlist], function(x) sum(is.na(x)))
vars_to_impute = varlist[to_impute_missCount>0]
vars_to_impute
mySL.library = c("SL.mean", "SL.glmnet", "SL.knn", "SL.randomForest") 
mySL.library = c("SL.randomForest","SL.glmnet") 

#fit
clinical_full_clean$is.Rash = clinical_full_clean$Rash %in% c(1)
clinical_full_clean[which(is.na(clinical_full_clean$Rash)),"is.Rash"] = NA
clinical_full_clean$is.Ascitis = clinical_full_clean$Ascitis %in% c(1)
clinical_full_clean[which(is.na(clinical_full_clean$Ascitis)),"is.Ascitis"] = NA
mean(as.numeric(clinical_full_clean$is.Ascitis),na.rm=T) #correct
var = "is.Ascitis"
XD = clinical_full_clean #restart

#RF 
  #RF will assume classification if Y is factor and will otherwise assume regression
    #even logical (boolean) vectors it will run regression on, so much convert to a number and then a factor
set.seed(1)
rf = randomForest(x=XD[which(!is.na(XD[var])),c(cvars_noMiss)], y=as.factor(as.numeric(XD[which(!is.na(XD[var])),var])),  
                  ntree=501, na.action=na.fail) 
#fill var with predicted value when var is missing
tis_r_as = predict(rf, XD[which(is.na(XD[var])),cvars_noMiss], type="response") 
XD = clinical_full_clean #restart
tis_p_as = predict(rf, XD[which(is.na(XD[var])),cvars_noMiss], type="prob") #not meaningful for regression



#SL
             
XD = clinical_full_clean #restart
#SL.mean gives a number between 1 and 2 when I use Rash as a factor, but between 0 and 1 when I use is.Rash as a logical
  #safer to use logical variables
var="is.DEN" #Pulso is 1,2, and 3 factor levels

#knn did not work for Rash(factor) nor for is.Rash(logical) -- 
    #todo: solve problem so we can use knn here too (see stuff on SL wrappers).  


###############################################################################
###############################################################################

cvars = get_clinic_var_list(clinic_varsD, outcome="ND.vs.DEN", eliminate_vars_with_missings=F, 
                            eliminate_constant_vars=T, eliminate_vars_with_minXnomiss=50,
                            XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T, UltraX=T, BloodLab=T)
XD = clin_full_wImputedRF1
#mySL.library = c("SL.mean","SL.randomForest","SL.glmnet","SL.polymars","my.SL.knn")  #,  did not work
mySL.library = c("my.SL.caret","SL.mean","my.SL.knn")
XD$DEN_prob = NULL
#cl <- makeCluster(16, type = "FORK") # can use different types here
#clusterSetRNGStream(cl, iseed = 2343)
SLresults = SuperLearner(Y=as.numeric(XD[,"is.DEN"]),
                         X=subset(XD, select=c(cvars)), 
                         family=binomial(), SL.library=mySL.library,
                         method = "method.NNLS", verbose=FALSE)
#stopCluster(cl)
SLresults
XD[,"DEN_prob"] = predict(SLresults, subset(XD, select=c(cvars)))
XD[c(1:20),c("is.DEN","DEN_prob")] #check results

###############################################################################
###############################################################################

#trees
predictors = DENvarlist_genOnly
XD = clin_full_wImputedRF1
cfit = ctree(is.DEN~., data=XD[,c(predictors,"is.DEN")])
predicted = predict(cfit, XD[,c(predictors)], type="response")
is.vector(predicted)
predicted[c(1:20),]
plot(cfit)

mySL.library = c("SL.mean","SL.glmnet")
SuperLearner.CV.control = list(10L, stratifyCV = FALSE, shuffle = TRUE,validRows = NULL)
var="is.DEN"
mySL = SuperLearner(Y=as.numeric(XD[which(!is.na(XD[var])),var]), X=XD[which(!is.na(XD[var])),c(predictors)], 
                    SL.library=mySL.library, family=binomial(), cvControl=SuperLearner.CV.control, method = "method.NNloglik")
#XD[which(is.na(XD[var])),var] = predict(mySL, newdata=XD[which(is.na(XD[var])),c(cvars_noMiss)], onlySL=T)
#for binary variables, turn predicted probabilities into assignments
prelim = predict(mySL, newdata=XD[c(1:10),c(predictors)],onlySL=T)$pred  #pred will give us predicted values.  otherwise, prelim contains fits for all libraries also.
prelim[c(1:20),]
XD[c(1:10),"DEN_dum"]
#XD[c(1:10),"is.Rash_pred3"] = prelim
#is.logical(XD$is.Rash_pred3) #passed the check

SLresults = SuperLearner(Y=as.numeric(trainD[, "DEN_dum"]), 
                         X=as.data.frame(trainD[,c(predictors)]), newX=as.data.frame(testD[,c(predictors)]),
                         family=binomial(), SL.library=mySL.library, cvControl=SuperLearner.CV.control,
                         method = "method.NNloglik", verbose=FALSE) #9-28-14: changed method from NNLS to NNloglik (Alan's suggestion)
mydata[folds[[v]],"DEN_prob"] = SLresults$SL.predict #predicts for obs in newX, as we want


##########################################
#todo:  
  #convert all binary variables into logical vectors (including outcomes)
  #add something in above code that says if is.logical=false and factor=false then family=gaussian, otherwise family=binomial
  #add the above SL code into the imputation function
  #run it to fill in values for missings!
##########################################  
#todo: explore multicore use with mcSuperLearner


#obtain list of clinical variables that have no missing values in data
clinic_vars_missCount = sapply(XD[,c(as.character(clinic_vars))], function(x) sum(is.na(x)))
cvars_noMiss = names(which(clinic_vars_missCount==0))
clinicCount = length(cvars_noMiss)
cat("\n clinical variables that shall be included in regressions \n")
print(cvars_noMiss) #view 
cat("\n variables excluded due to missing values \n")
print(names(which(clinic_vars_missCount>0))) #view 
allCount = clinicCount + MFcount

#obtain list of clinical variables that includes the dummy variables for when using imputed values


#predict XD[which(is.na(XD[var])),var] = 
prelim = predict(mySL, newdata=XD[which(is.na(XD[var])),c(cvars_noMiss)], 
                 X=XD[which(!is.na(XD[var])),c(cvars_noMiss)], 
                 Y=as.numeric(XD[which(!is.na(XD[var])),var]), onlySL=T)[[1]]

################## set up data #############################

#create a version of the data which replaces values for missing compounds with zeros
XD = process_data(comboD1_filter50n)
#obtain list of candidate clinical variables to include in ND vs. DEN prediction
clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                          "Variable.Name.in.Hospital.data"]),"age","DaysSick") #43

#obtain list of MFs 
cindices = grep("MZ_",colnames(XD)) #column indices in dataset
allMFs = colnames(XD[cindices]) #names
MFcount = length(allMFs)

#obtain list of clinical variables that have no missing values in data
clinic_vars_missCount = sapply(XD[clinic_vars], function(x) sum(is.na(x)))
cvars_noMiss = clinic_vars[clinic_vars_missCount==0] #40 
clinicCount = length(cvars_noMiss)
print("clinical variables that shall be included in regressions")
print(cvars_noMiss) #view 
print("variables excluded due to missing values")
print(clinic_vars[clinic_vars_missCount>0]) #view 
allCount = clinicCount + MFcount

#knn (for super learner) also has problems with factors variables with only 2 levels 
#so must convert to numeric
YesNo_varsH = clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1 & clinic_varsD$Variable.Output=="Yes, No"),
                           "Variable.Name.in.Hospital.data"]
YesNo_varsH_char = sapply(YesNo_varsH, as.character)
data_numeric = convert_factor_to_numeric(XD, c("DEN_dum","Torniquete","Sexo",YesNo_varsH_char))

###############################################################################

set.seed(100)
rf = randomForest(x=XD[,c(cvars_noMiss, allMFs)], y=XD[,"DEN_dum"], ntree=500, localImp=FALSE)
rf
set.seed(100)
rf = randomForest(x=data_numeric[,c(cvars_noMiss, allMFs)], y=XD[,"DEN_dum"], ntree=500, localImp=FALSE)
rf
set.seed(100)
rf = randomForest(x=data_numeric[,c(cvars_noMiss, allMFs)], y=data_numeric[,"DEN_dum"], ntree=500, localImp=FALSE)
rf

#conclusions: rf does not give error so long as the outcome is a factor.  
#however, rf does give different answer if predictors are numeric as opposed to factors

###############################################################################


#Super Learner
#note: not all algorithms are estimable given that I have more predictors than observations 
#gam, glm will not run when all clinical variables are included
mySL.library = c("SL.mean", "SL.randomForest", "SL.glmnet", "SL.knn")
mySL.library = c("SL.mean", "SL.randomForest", "SL.knn")

mynumeric = data_numeric[,c(cvars_noMiss, "DEN_dum")]
mydata$DEN_prob = NULL
SLresults = SuperLearner(Y=mynumeric[, "DEN_dum"],
                         X=subset(mynumeric, select=-c(DEN_dum)), 
                         family=binomial(), SL.library=mySL.library,
                         method = "method.NNLS", verbose=FALSE)
SLresults
mydata[,"DEN_prob"] = predict(SLresults, subset(mynumeric, select=-c(DEN_dum)))
mydata$DEN_prob

VIPs = c("MZ_114.09146","MZ_315.1347","MZ_412.2911","MZ_129.0546","MZ_734.64233")
VIPs = c("MZ_114.09146") #this produces error with glmnet
mynumeric = data_numeric[,c(VIPs, "DEN_dum")]
SLresults = SuperLearner(Y=mynumeric[, "DEN_dum"],
                         X=subset(mynumeric, select=-c(DEN_dum)), newX=subset(mynumeric, select=-c(DEN_dum)),
                         family=binomial(), SL.library=mySL.library,
                         method = "method.NNLS", verbose=FALSE)
mynumeric$DEN_prob = SLresults$SL.predict
#calculate mean squared error
mynumeric$sq_err = (mynumeric$DEN_prob - mynumeric$DEN_dum)^2
MSE = mean(mynumeric$sq_err)
MSE
MSE = mean( (as.vector(mynumeric$DEN_prob) - as.vector(mynumeric$DEN_dum)^2 )
MSE
t1 = as.data.frame(t(c(1,SLresults$coef)))
t1
t2 = as.data.frame(t(c(2,SLresults$coef)))
t2
t4 = smartbind(t1,t2,by="V1", all=TRUE)
t4

predict(SLresults, subset(mynumeric, select=-c(DEN_dum)))

shortD = mydata
train <- subset(shortD[seq(1, dim(shortD)[1], 2), ], select=-c(DEN_dum))
test <- subset(shortD[seq(2, dim(shortD)[1], 2), ], , select=-c(DEN_dum)) 
cl = as.numeric(subset(shortD[seq(1, dim(shortD)[1], 2), ], select=c(DEN_dum))$DEN_dum)
knn(train, test, cl, k = 3, prob=TRUE)



###############################################################################

#Super Learner
#note: not all algorithms are estimable given that I have more predictors than observations 
#gam, glm will not run when all clinical varaibles are included
mySL.library = c("SL.mean", "SL.knn", "SL.randomForest", "SL.glmnet")
data_numeric = convert_factor_to_numeric(processedD, c("DEN_dum","Torniquete","Sexo",YesNo_varsH_char))
mydata= data_numeric[,c(cvars_noMiss, "DEN_dum")]

train <- mydata[seq(1, dim(shortD)[1], 2), ]
test <-  mydata[seq(2, dim(shortD)[1], 2), ]

SLresults = CV.SuperLearner(Y=mydata[, "DEN_dum"], 
                            X=subset(mydata, select=-c(DEN_dum)), newX=subset(test, select=-c(DEN_dum)),
                            family=binomial(), SL.library=mySL.library, V=10,
                            method = "method.NNLS", verbose=FALSE)
SLresults$SL.predict
mydata[seq(2, dim(shortD)[1], 2),"DEN_prob"] = SLresults$SL.predict



###############################################################################

#a matrix with one row for each input observation and one column for each class, 
#giving the fraction or number of (OOB) ‘votes’ from rf

#processedD = process_data(clinic_varsD)
processedD = process_data(comboD1_filter50n)
cindices = grep("MZ_",colnames(processedD)) #column indices in dataset
allMFs = colnames(processedD[cindices]) #names

clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                            "Variable.Name.in.Hospital.data"]),"age","DaysSick") #43

#obtain list of clinical variables that have no missing values in data
clinic_vars_missCount = sapply(processedD[clinic_vars], function(x) sum(is.na(x)))
cvars_noMiss = clinic_vars[clinic_vars_missCount==0] #40 
clinicCount = length(cvars_noMiss)
print("clinical variables that shall be included in regressions")
print(cvars_noMiss) #view 
print("varibles excluded due to missing values")
print(clinic_vars[clinic_vars_missCount>0]) #view 


rf = randomForest(DEN_dum ~ MZ_322.30972 + MZ_352.35504 + MZ_350.34042 + MZ_324.32462 + 
                    MZ_170.15321, data=processedD, ntree=500)
votes(rf) 
processedD$DEN_prob = predict(rf, processedD, type="prob")[,2]
processedD[,c("DEN_dum","DEN_prob")] #check predictions

###############################################################################

#Put clinical features and MFs in RF model and select the top Xx non-clinical features from the list of all top features
#partition data into 10 equally sized peices
mydata=processedD
nobs = nrow(mydata)

set.seed(183)
partition_nobs <- ceiling(0.1 * nobs) #size of each peice
nvec = rep(seq(from=1, to=10),partition_nobs)[1:nobs] #vector of values 1 - 10 of length nobs with [approximately] equal numbers of each value
rvec = permute(nvec) #random permutation of vector

full_list = list()
for(i in seq(from=1, to=10)) {
  trainD = mydata[which(rvec!=i),] #training set for first run 
  testD = mydata[which(rvec==i),]  #test set for first run
  
  rf = randomForest(x=processedD[,c(cvars_noMiss, allMFs)], y=processedD[,"DEN_dum"], ntree=500, localImp=FALSE)
  rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
  rf_VIPS_ordered = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
  rf_VIP_MFs = names(rf_VIPS_ordered[grep("MZ_",names(rf_VIPS_ordered))])
  MFpreds = rf_VIP_MFs[1:1]
  MFpreds
  MFpreds = names(rf_VIPS[order(rf_VIPS, decreasing=TRUE)[1],])
  full_list = append(full_list, MFpreds)
}
full_list
list_freqs = subset(as.data.frame(table(t(subset(as.data.frame(table(full_list)), select=-c(Freq)))[,1])))
list_freqs
list_freqs_ordered2 = list_freqs[order(list_freqs$Freq, decreasing=TRUE),]

FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
m = merge(list_freqs_ordered1, list_freqs_ordered2, by="Var1", all=TRUE)
m$Freq.x= sapply(m$Freq.x,FUN=FillMissingsWithZeros)
m$Freq.y= sapply(m$Freq.y,FUN=FillMissingsWithZeros)
m$total = m$Freq.x + m$Freq.y
m_ordered = m[order(m$total, decreasing=TRUE),]
m_ordered

#rfcv: gives CV-validated prediction performance of models with sequentially reduced number of predictors 
#(ranked by variable importance)
result <- rfcv(trainx=, trainy=, cv.fold=5, scale="log")
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#rfImpute: imputs missing values.  If rf is run using this filled in data, error rate will likely be optimistic.
#tuneRF: finds optimal mtry parameter value.
#varImpPlot

###############################################################################

#Use the lasso to obtain best subset (alpha=1 says to use lasso penalty; alpha=0 corresponds to ridge)
#penalty.factor: can apply separate penalty factors to each coefficient.  
#This is a number that multiplies lambda to allow differential shrinkage
#setting this to zero for the clinical variables should gaurenttee their inclusion, could cause overfitting.
#pmax: limit the number of variables to be non-zero
#glmnet produces the regularization path over a grid of values for the tuning parameter lambda
fit1=glmnet(x=as.matrix(processedD[cindices]),y=as.matrix(processedD[,"DEN_dum"]),
            family="binomial", alpha=1, standardize=TRUE)
plot(fit1)  #visual display of coefficient values (neat looking, though not super important)
#the s parameter sets the lambda value used for prediction
predict(fit1,newx=as.matrix(processedD[cindices]), s=0.001, type="response")
plot(fit1$lambda, fit1$df, main="Number of non-zero coefficients as function of tuning parameter") 

#cv.glmnet chooses the optimal tuning parameter based on CV.
#type.measure could also be set to "auc", which would find the lambda that minimizes the area under the ROC
#I do not see documentation which allows penality.factor to be used with cv.glmnet, 
#so I might need to find the optimal lamba myself using glmnet with CV
set.seed(456)
cvLasso = cv.glmnet(x=as.matrix(processedD[,c(cvars_noMiss, allMFs)]),
                    y=as.matrix(as.numeric(processedD[,"DEN_dum"])),
                    family='binomial', alpha=1, nfolds=10, type.measure="deviance")
cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
MFpreds = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs (MFs with non-zero coef)

#lasso did not seem to work with factor variables
cvLasso = cv.glmnet(x=as.matrix(data_numeric[,c(cvars_noMiss, allMFs)]),
                    y=as.matrix(data_numeric[,"DEN_dum"]),
                    family='binomial', alpha=1, nfolds=10, type.measure="deviance")
plot(cvLasso) #shows deviance as a function of the tuning parameter lambda
#lambda.1se if the largest value of lambda such that error is within 1 std err of the minimum
#lambda.min is the value of lambda that gives the minimum  mean CV error
cvLcoefs = as.matrix(coef(cvLasso, s="lambda.min")) #coefficients when using lambda that minimizes deviance
LassoVIP = names(cvLcoefs[which(cvLcoefs>0),][-1]) #names of selected MFs (MFs with non-zero coef)
#put in descending order the non-zero coefficients, excluding intercept
nonzeros = cvLcoefs[which(cvLcoefs>0),][-1]
LasSigX = nonzeros[order(nonzeros)]
LasSigX_MFs = names(LasSigX[grep("MZ_",names(LasSigX))])
LasSigX_MFs

#Include clinical features and take the non-clinical terms that have non-zero coefficients. 

print(paste("Selected", topX, "compounds based on ttest:"))
print(row.names(sigX))
LassoVIP
cvLasso$cvm #mean CV error
predict(cvLasso,newx=processedD[cindices], s="lambda.1se")
min(cvLasso$cvm)




############################################################### 
#### dengue (DF/DHF/DSS) vs. non-dengue (ND) prediction ####### 
############################################################### 

#create dummy for non-dengue (0) versus dengue (1)
nrows_DF = dim(wzeros_comboD)[1] #number of rows
dummyDF = vector(mode="numeric", length=nrows_DF)
dummyDF[which(wzeros_comboD["DxFinal4cat"]=="ND")] = 0
dummyDF[which(wzeros_comboD["DxFinal4cat"]!="ND")] = 1
wzeros_comboD$dummyDF = dummyDF
class(wzeros_comboD$dummyDF) 

#obtain list of imputed indicators in data that match to covarlist_genOnly
myvars = colnames(clin_full_wImputedRF1)
imputed_vars = myvars[grepl("_imputed", myvars)] #all imputed indicators in data
varlist = paste(covarlist_genOnly, "_imputed", sep="") 
imputed_vars_matched = imputed_vars[imputed_vars %in% varlist==TRUE] 
#drop redundant indicators
XD$is.cohort = (XD$Study=="Cohort")
is.diff_fun = function(var, XD) (sum(XD$is.cohort!=XD[,var]))
diffs = sapply(imputed_vars_matched, FUN=is.diff_fun, XD)
unique_imp_matched = names(diffs[diffs!=0])
cvars_noMiss = c(cvars_noMiss, unique_imp_matched)

#explanatory_vars

#get the V-fold cross-validated risk estimate for super learner.
set.seed(999)
#SL.gam, SL.glm produced errors (singularities and non-convergence)
mySL.library = c("SL.mean", "SL.knn", "SL.randomForest")
SLresultsDF = CV.SuperLearner(Y=dummyDF, X=wzeros_comboD[explanatory_vars], 
                              V=10, family=binomial(), SL.library=mySL.library,
                              method = "method.NNLS", verbose=FALSE)
SLresultsDF
summary(SLresultsDF)
coef(SLresultsDF)



#Plot the CV risk to sample size to see if it’s starting to plateau.




############################################################
###### severe dengue (DSS) vs. non-severe dengue (DHF) #####
############################################################

#drop the ND and DF folks
DSS_and_DHF = comboD[which(comboD["diagnosis"]=="DSS" | comboD["diagnosis"]=="DHF"),]
nrows_DSS = dim(DSS_and_DHF)[1] #number of rows
dummyDSS = vector(mode="numeric", length=nrows_DSS)
dummyDSS[which(DSS_and_DHF["diagnosis"]=="DHS")] = 0
dummyDSS[which(DSS_and_DHF["diagnosis"]=="DSS")] = 1
DSS_and_DHF$dummyDSS = dummyDSS
DSS_and_DHF[, c("diagnosis","dummyDSS")]
class(DSS_and_DHF$dummyDSS) 

#drop variables with one or more missing values (SL cannot deal with missings)
#better to input missing values? (when/why do missings occur with LC-MS?)
dim(DSS_and_DHF)[2]  #number of columns before dropping those with missing values
missingsDSSDHF = sapply(DSS_and_DHF, function(x) sum(is.na(x)))
DSS_and_DHF = DSS_and_DHF[,which(missingsDSSDHF==0)]
dim(DSS_and_DHF)[2]  #number of columns after dropping those with missing values
#covariates drop from 744 to 448

#obtain list of explanatory variables for super learner
myvars = colnames(DSS_and_DHF)
explanatory_vars = myvars[grepl("MZ_", myvars)]
#explanatory_vars

#get the V-fold cross-validated risk estimate for super learner.
set.seed(999)
mySL.library = c("SL.glm", "SL.randomForest")
SLresultsDSS = CV.SuperLearner(Y=dummyDSS, X=DSS_and_DHF[explanatory_vars], 
                               V=10, family=binomial(), SL.library=mySL.library,
                               method = "method.NNLS", verbose=FALSE)
SLresultsDSS
summary(SLresultsDSS)
coef(SLresultsDSS)