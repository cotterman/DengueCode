
###################### Apply simple learners to Dengue data ###################


###############################################################################
######################## Establish directories and data #######################
###############################################################################

rm(list = ls()) #start with blank slate (clear everything from workspace)

library(rpart) #Recursive partitioning for classification, regression and survival trees
library(survival)
library(LogicReg) #Logic regression (to find predictors that are boolean combos of original predictors)

#### Establish directories ####
inputsDir = "/home/carolyn/dengue_dx_not_backed_up/Dengue_data/Processed/R_inputs/"
outputsDir = "/home/carolyn/dengue_dx_not_backed_up/Dengue_data/Processed/R_outputs/"
resultsDir = "/home/carolyn/dengue_dx_not_backed_up/R_results/"


###############################################################################
######################## Load and clean data to use ###########################
###############################################################################

process_data = function(mydata){
  cindices = grep("MZ_",colnames(mydata))
  processedD = mydata
  #note: must have a return statement b/c assignment operator will not return anything (and default is 0)
  FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
  processedD[cindices] = lapply(mydata[cindices],FUN=FillMissingsWithZeros)
  return(processedD)
}

### Option 1: Use the full clinical/lab dengue data (most will have no LCMS data)
clinic_varsD = read.delim(paste(inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
source("/home/carolyn/temp_Dengue_code/clean_clinical_data.R") #produces clinical_comboD
#drop observations with unknown final dengue dx
dataf = clinical_comboD[which(!is.na(clinical_comboD$DxFinal4cat)),]
save(dataf, file=paste(outputsDir,"dataf.RData",sep=""))
load(paste(outputsDir,"dataf.RData",sep="")) #loads dataf dataframe

### Option 2: Use only observations that have LCMS data
load(paste(outputsDir,"clinical_comboD_clean.RData",sep="")) 
dataf = clinical_comboD_clean

### Option 3: Use LCMS data
load(paste(outputsDir,"comboD1_filter50n.RData",sep="")) #loads comboD1_filter50n dataframe
dataf = comboD1_filter50n

#create a version of the data which replaces values for missing compounds with zeros
dataf = process_data(dataf)


###############################################################################
######################## Obtain list of predictors  ###########################
###############################################################################

#obtain list of clinical variables that have no missing values in data
outc = "DEN_dum" #outcome is just ND vs. DEN
clinic_varsD = read.delim(paste(inputsDir,"List of clinical variables for analysis.txt", sep=""), header=TRUE, nrows=500)
clinic_vars = c(as.character(clinic_varsD[which(clinic_varsD$Use.in.ND.vs.DEN.prediction==1),
                                          "Variable.Name.in.Hospital.data"]),"age","DaysSick")
clinic_vars_missCount = sapply(dataf[clinic_vars], function(x) sum(is.na(x)))
cvars_noMiss = clinic_vars[clinic_vars_missCount==0] #use these 32 predictors

#obtain list of MFs 
cindices = grep("MZ_",colnames(dataf)) #column indices in dataset
allMFs = colnames(dataf[cindices]) #names
MFcount = length(allMFs)

#choose which predictors to use ( cvars_noMiss or c(cvars_noMiss, allMFs) )
#c_preds = cvars_noMiss
c_preds = c(cvars_noMiss, allMFs)
predn = c_preds #full names of predictor variables (same as var names in c_preds for now)


###############################################################################
############################ Ordinary CART ####################################
###############################################################################

Y=dataf[,outc]
X=dataf[,c_preds]
fit.rpart=rpart(Y ~ ., data=as.data.frame(X),
                control = rpart.control(cp = 0.001, minsplit = 20, 
                                        xval = 20, minbucket = 10), 
                method = "class") #, cost = rep(0.000001,xp)
fit.rpart
summary(fit.rpart)


###############################################################################
#################### Oridinary Logic regression ###############################
###############################################################################



###############################################################################
####################### Catepillar Learner ####################################
###############################################################################
# Based on R code from Alan (emailed 9-28-2014)

nalpha=10
alpha=seq(1.1,2,length=nalpha)

oper=NULL
xs=NULL
setsx=NULL
varn=NULL
outnames=NULL

Y=dataf[,outc]
X=dataf[,c_preds]
nx=length(c_preds)
dat.tmp=data.frame(Y,X)

#loop thru each predictor variable
for(i in 1:nx){
    #i = 7 #the 7th predictor (rash) is one of the chosen covaraites used for splitting
    Xv=dat.tmp[,(i+1)]
    Xv=as.numeric(Xv)
    Yv = dat.tmp[,1]
    dat.tmp2=as.data.frame(cbind(Yv,Xv))
    
    #finds best split (while satisfying various specified tuning parameters)
    fit.rpart=rpart(Yv ~ ., data = dat.tmp2, control = rpart.control(cp = 0.001, minsplit = 20, 
        xval = 20, maxdepth = 1, minbucket = 10), method = "class", cost = 0.000001)
    xsplt=fit.rpart$splits

    if(is.null(xsplt)) {
        xs=c(xs,NA)
        oper=c(oper,"")
    }

    if(is.null(xsplt)==F) {
        xsplt=xsplt[1,"index"]	
        x1=dat.tmp2[,2]
        ix1=as.numeric(x1>xsplt)
        tt=table(dat.tmp2[,1],ix1)
        y1=tt[2,1]/sum(tt[,1])
        y2=tt[2,2]/sum(tt[,2])
        optmp = ifelse(y1>y2,"<",">")
        xs=c(xs,xsplt)
        oper=c(oper,optmp)
        pp=predict(fit.rpart,type="class")
        dat.tmp=dat.tmp[pp==1,] 
        #next time we loop, exclude observations for which we have already predicted dengue
    }
}
#so end result is dianostic test that goes as follows: 
  #do you have a?  if yes, then positive.  if no, then do you have b?
  #if yes then positive, if no then do you have c?  if yes then positive, if no then no dengue


dd=data.frame(outnames=rep(outc,nx), setsx=predn, cut=xs, oper=oper)
dd #view results




###############################################################################
########## Logic regression on discretized covariate space ####################  
###############################################################################
# Based on R code from Alan (emailed 9-28-2014)
# Can turn into a "superlearner" by discretizing based on set of tuning parameters 

sum.na=function(x){sum(is.na(x))}

varn=c_preds
Y=dataf[,outc]
X = dataf[,c_preds]
#convert predictors into numeric variables
for (name in c_preds) {
  X[,name] = as.numeric(dataf[,name])
} 
xp=dim(X)[2]
filepre=paste("OUTis_", outc, "_PREDSETis_", predn, sep="")
filepre
newX=NULL
coln=NULL

#obtain quantiles
qt=apply(X, 2, quantile, probs=seq(0.1,0.9,0.2))

#for each predictor, create a dummy for each quartile (except for first)
for(k in 1:xp) {
    #k = 39 #use age as a test case
    #xn will be predictor values, recoded into quartiles
    xn=cut( X[,k],breaks=unique(c(min(X[,k])-0.1, qt[,k], max(X[,k]+0.1))) )
    #inds is a matrix with a dummy variable for each quartile (except the first one)
    inds <- model.matrix(~ factor(xn) - 1)[,-1]
    #create sensible names for new data columns
    nmes=colnames(inds)
    nc=nchar(nmes)
    nmes=paste(varn[k],substr(nmes,11,nc),sep="")
    coln=c(coln,nmes)
    newX=cbind(newX,inds)
}
#piece together all of our dummies
colnames(newX)=coln
ncol=dim(newX)[2]
#select=1 means to fit a single model (tree)
logicreg=logreg(Y, newX, select=1)
logicreg
#produce visual of results
tit=paste("Predictors are ",predn," Outcome is ",outc,sep="")
pdf(paste(filepre,"pdf",sep=".")) #how to get this written to specified folder??
plot.logregmodel(logicreg$model,nms=coln,nny=0.25)
mtext(tit, side = 3, line = 1,cex=1.25)
dev.off()

