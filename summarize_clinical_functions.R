
#function returns cross tabulation with col and row totals
tabletot = function(mydata, var1, var2, useNA="no"){
  t = table(mydata[,var1], mydata[,var2], useNA=useNA)
  t = cbind(t, rowSums(t))
  t = rbind(t, colSums(t))
  print(t)
}

#function takes data column and returns variance of that column.  will convert factors to numeric values first. 
var_after_conversion = function(x){
  #R's built-in variance function only works for numeric variables - convert factors and characters to numeric
  if(is.factor(x)){
    #print("make")
    variance = var(as.numeric(x), na.rm=TRUE)
  }
  else{
    #print("yes")
    variance = var(x, na.rm=TRUE)
  }
  return(variance)
}

#function considers "clinic_vars" list and returns a list of variables that take the same value for all obs in data XD
#returns empty list of none of the "clinic_vars" are constant throughout XD
get_constant_var_list = function(XD, clinic_vars){
  variances = sapply(XD[,clinic_vars], var_after_conversion)
  varnames = names(variances[which(variances==0)])
  return(varnames)
}

summarize_clinical = function(mydata){
  print(summary(mydata))
  
  #days between symptom onset and sample collection
  print(table(mydata$DaysSick))
  
  #diagnosis variables
  print("Res_final: final DENV test result (determined by lab)")
  print(table(mydata$Res_Final)) #no missings
  print("ClasificacionFinal (without regard for lab result)")
  print(table(mydata$ClasificacionFinal)) #no missings
  print("DxFinal4cat: considers all info")
  print(table(mydata$DxFinal4cat)) #no missings
  
  print("All samples: ClasificacionPrimerDia by DxFinal4cat")
  print(table(mydata$ClasificacionPrimerDia, mydata$DxFinal4cat))
  print("All samples: ClasificacionFinal by DxFinal4cat")
  print(table(mydata$ClasificacionFinal, mydata$DxFinal4cat)) 
  
  print("Restrict to just LC-MS serum, run 1")
  print(table(mydata[which(mydata$serum==1),"ClasificacionPrimerDia"], 
              mydata[which(mydata$serum==1),"DxFinal4cat"]))
  print(table(mydata[which(mydata$serum==1),"DaysSick"])) 
  
  print("restrict to just LC-MS serum, run 2")
  print(table(mydata[which(mydata$serum==2),"ClasificacionPrimerDia"], 
              mydata[which(mydata$serum==2),"DxFinal4cat"]))
  print(table(mydata[which(mydata$serum==2),"DaysSick"])) 
  
  print("restrict to just LC-MS saliva")
  print(table(mydata[which(mydata$saliva==3),"ClasificacionPrimerDia"], 
              mydata[which(mydata$saliva==3),"DxFinal4cat"])) 
  print(table(mydata[which(mydata$saliva==3),"DaysSick"])) 
  
  print("restrict to just LC-MS urine")
  print(table(mydata[which(mydata$urine==5),"ClasificacionPrimerDia"], 
              mydata[which(mydata$urine==5),"DxFinal4cat"]))
  print(table(mydata[which(mydata$urine==5),"DaysSick"])) 
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
    XD$is.DEN = XD$DHF_DSS_dum
  }
  
  #obtain list of MFs 
  cindices = grep("MZ_",colnames(XD)) #column indices in dataset
  allMFs = colnames(XD[cindices]) #names
  
  if(dim_reduce_covar==F){
    rf = randomForest(x=XD[,c(allMFs)], y=as.factor(XD[,"is.DEN"]), ntree=501, localImp=FALSE)
  }else if(dim_reduce_covar==T){
    #obtain list of clinical variables that have no missing values in data
    clinic_vars_missCount = sapply(XD[clinic_vars], function(x) sum(is.na(x)))
    cvars_noMiss = clinic_vars[clinic_vars_missCount==0]  
    rf = randomForest(x=XD[,c(cvars_noMiss, allMFs)], y=as.factor(XD[,"is.DEN"]), ntree=501, localImp=FALSE)
  }
  rf_VIPS = importance(rf, type=2) #type=1 uses the permutation method.  not sure how to interpret negative numbers
  rf_VIPS_ordered = rf_VIPS[order(rf_VIPS, decreasing=TRUE),]
  rf_VIP_MFs = rf_VIPS_ordered[grep("MZ_",names(rf_VIPS_ordered))] #exclude non-MF variables from list
  rf_VIP_MFs_formated = as.data.frame(t(t(rf_VIP_MFs)))
  colnames(rf_VIP_MFs_formated)=mycolname
  
  return(rf_VIP_MFs_formated)
}
