
### function to give us the duplicated values in input vector
get_dupIDs = function(vec, vec_name){
  #verify there are no duplicate code values in either dataset
  t1 = table(vec)
  v1dups = table(t1) 
  if(dim(v1dups)>1){
    print(paste("Counts of duplicates in", vec_name))
    print(v1dups)
    print(paste("Duplicated IDs in", vec_name))
    print(t1[which(t1>1)])
  }else{
    print("No duplicated values")
  }
}

### function to give us the ID of observations which are not found in both datasets
    # also checks ID duplicates in each dataset
get_mismatches = function(v1, v2, v1_name, v2_name){
  #v1 = respD[,"code"]
  #v1_name = "LCMS"
  #v2 = study_info[,"code"]
  #v2_name = "study info"
  
  get_dupIDs(v1, v1_name)
  get_dupIDs(v2, v2_name)
  
  #check out values that are not in both datasets
  v1_not_v2 = v1[v1 %in% v2==FALSE] #v1 values not found in v2
  v2_not_v1 = v2[v2 %in% v1==FALSE] #v2 values not found in v1
  if(length(v1_not_v2)>0){
    print(paste(v1_name,"values not found in",v2_name))
    print(v1_not_v2)
  }
  if(length(v2_not_v1)>0){
    print(paste(v2_name,"values not found in",v1_name))
    print(v2_not_v1)
  }
}

### function to read in the file indicating which study each sample belongs to ###
get_study = function(respD_only, lcms_run){
  
  # process file that contains indicator of study (cohort versus hospital)
  if(lcms_run==2){
    study_info = read.delim(paste(inputsDir,"Study classifications_batches 3 and 4.txt", sep=""), header=TRUE, nrows=200)
  } else if(lcms_run==3){
    study_info = read.delim(paste(inputsDir,"Study classifications_saliva.txt", sep=""), header=TRUE, nrows=200)
  } else if(lcms_run==5){
    study_info = read.delim(paste(inputsDir,"Study classifications_urine.txt", sep=""), header=TRUE, nrows=200)
  }
  newnames = as.character(c("Study","code","Res_Final2","OMS"))
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  fcode =    function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
  code_list = apply(X=study_info, MARGIN = 1, FUN=fcode)
  study_info["code"] = code_list
  
  #merge with lc-ms data
  mismatches = get_mismatches(v1=study_info[,"code"], v2=respD_only[,"code"], v1_name="study info", v2_name="LCMS")
  respD_new = merge(study_info[,c("Study","code")], respD_only, by="code") #this returns the intersection of each dataset
  get_dupIDs(respD_new[,"code"], "LCMS") #may still have duplicates
  
  return(respD_new)
}


### function to read in the LC-MS text file and organize it into dataframe ###
clean_LCMS = function(infile, lcms_run, roundme=FALSE, decimals=2, printme=FALSE) {
  
  #infile="LCMS_saliva_Nica_50percent.txt"
  #infile = "LCMS_urine_Nica_50percent.txt"
  #lcms_run=5
  #printme=TRUE
  #decimals=2
  #roundme=FALSE
  
  abundancesD1 = read.delim(paste(inputsDir, infile, sep=""), header=TRUE, nrows=200)
  nColumns = dim(abundancesD1)[2] #number of columns in abundance data
  #MZ = mass to charge of the compound, it usually is approximately +1 of the mass.  
  #This value distinguishes the type of metabolite.  
  #RT = retention times.
  #Resp = response (aka abundance).
  
  #obtained list of compounds (in same order in which they appear in spreadsheet columns)
  # take first row of values for variables MZ, MZ.1, MZ.2.,,,,MZ.2232
  MZ_Cols = seq(from=2, to=nColumns-2, by=3)
  #MZ_Cols #these are the column numbers that we will want to keep
  MZ_Nums = abundancesD1[1,MZ_Cols] #take the first row of data containing MZ values
  #MZ_Nums #these are the values in the specified columns
  
  #view MZ values in histogram
  MZ_Nums_reshaped1 = t(MZ_Nums) #min value is 107, max is 1427
  if(roundme==FALSE) {
    MZ_Names = paste("MZ_", sapply(MZ_Nums, as.character), sep="") #no rounding
  }
  takeDdecimals = function(x) round(x, digits=decimals)
  MZ_Nums_short = sapply(X=MZ_Nums, takeDdecimals)
  if(roundme==TRUE) {
    MZ_Names = paste("MZ_", sapply(MZ_Nums_short, as.character), sep="") #with rounding
  }
  
  #drop the RT and MZ columns and relabel the Resp columns with MZ names
  RespCols = seq(from=4, to=nColumns, by=3)
  respD = abundancesD1[,c(1,RespCols)]
  newnames = as.character(c("code", MZ_Names))
  #Note: when MZ values are repeated, R will give the repeated MZs a .1, .2, etc. suffix
     #ex: MZ_336.32513 and MZ_336.32513.1 and MZ_336.32513.2 will be in the data if MZ_336.32513 occurs 3 times
  colnames(respD) = newnames #rename columns
  
  #reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
  formatID = function(x, start, stop) sprintf("ID%04d",as.numeric(substr(as.character(x["code"]), start=start, stop=stop)))
  newIDs_1 = apply(X=respD, MARGIN = 1, FUN=formatID, start=13, stop=17)
  respD["code"] = newIDs_1
  
  #add an indicator of LCMS_run (to keep track once merged with other round of LC-MS)
  respD["LCMS_run"] = lcms_run
  
  #add an indicator for type of study (will need this later when we merge)
  if(lcms_run==1){
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run)
  } 
  
  #count  the number of duplicated MZ values
  if(printme==TRUE){
    print("Number of duplicated MZ values with no rounding of MZ values:")
    print(table(table(MZ_Nums_reshaped1))) #without rounding
    print(paste("Number of duplicated MZ values with rounding to", decimals, "decimals"))
    print(table(table(MZ_Nums_short))) #with rounding
  }
  return(respD)
  
}


### function takes abundance data as input and returns matrix of numeric MZ values (including duplicates if present)
get_MZ_vals = function(mydata){
  #obtain values of MZ numbers from processed abundance data
  temp = colnames(mydata) #all column names
  MZ_Names = temp[grep("MZ_",temp)] #take just the columns that start with "MZ_"
  #eliminate the extensions that indicate duplicate values --- these are always after the second decimal
  #(TODO: should verify there are no whole-number MZ values -- or reformat to eliminate whole numbers
  #will only be a problem if there are repeats of whole-number MZ values, as I wont detect them)
  #knock out the "MZ_" prefix and convert to number
  get_mz_nums = function(x){
    LocateDecimal2 = gregexpr(pattern ='\\.', x)[[1]][2]
    mystop = ifelse(is.na(LocateDecimal2), 50, LocateDecimal2-1)
    as.numeric(substr(x, start=4, stop=mystop))
  }
  MZ_Nums = lapply(X=MZ_Names, FUN=get_mz_nums) 
  MZ_Nums_reshaped = as.matrix(do.call(rbind.data.frame, MZ_Nums)) #convert list to matrix
  return(MZ_Nums_reshaped)
}



### function takes abundance data as input and produces histograms of MZ values
graph_MZ_frequencies = function(data1, data2, dname1, dname2, fileprefix) {
  
  MZ_Nums_reshaped1 = get_MZ_vals(data1)
  MZ_Nums_reshaped2 = get_MZ_vals(data2)
  
  
  #mybreaks = seq(min(MZ_Nums_reshaped1),max(MZ_Nums_reshaped1),10)
  #histogram(MZ_Nums_reshaped1, nint=100) #overview
  #histogram(MZ_Nums_reshaped1, xlim=c(300,400), nint=10000) #zoomed
  #histogram(MZ_Nums_reshaped1, xlim=c(330,380), nint=10000) #zoomed
  
  ### kernel density, overview
  png(paste("/home/carolyn/dengue_dx/R_results/",fileprefix,"MZ_density_full.png", sep=""))
  plot(density(MZ_Nums_reshaped1, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
       col="blue", xlab="", main="Kernel Density of MZ Values (full range)")
  par(new=TRUE)
  plot(density(MZ_Nums_reshaped2, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
       col="red", xlab="", main="") 
  legend("topright",legend=c(dname1, dname2), col=c("blue","red"),lty=1)
  dev.off()
  
  ### kernel density, zoomed in
  
  png(paste("/home/carolyn/dengue_dx/R_results/",fileprefix,"MZ_density_100.png", sep=""))
  plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
       col="blue", main="Kernel Density of MZ Values", xlab="") #kernal density
  lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
        col="red", main="", xlab="") #kernal density
  legend("topleft",legend=c(dname1, dname2), col=c("blue","red"),lty=1)
  dev.off()
  
  plot_kdensity = function(mz_start, mz_end, ylimit, fileprefix) {
    png(paste("/home/carolyn/dengue_dx/R_results/",fileprefix,"MZ_density_",mz_start,".png", sep=""))
    plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
         col="blue", main="", xlab="") #kernal density
    lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
          col="red", main="", xlab="") #kernal density
    dev.off()
  }
  plot_kdensity(mz_start=200, mz_end=300, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=300, mz_end=400, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=400, mz_end=500, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=500, mz_end=600, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=600, mz_end=700, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=700, mz_end=800, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=800, mz_end=900, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=900, mz_end=1000, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=1000, mz_end=1100, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1100, mz_end=1200, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1200, mz_end=1300, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1300, mz_end=1400, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1400, mz_end=1500, ylimit=1000, fileprefix)
  
}

### Basic summary of clinical and lab variables ###
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

### knn (run thru SuperLearner) does not like factor variables with 2 levels, so convert to numeric.
#convert all yes/no questions into factors (regression algorithms may run suboptimal if coded as numeric)
convert_factor_to_numeric = function(data, vars_to_convert){
  #print(vars_to_convert)
  for (name in vars_to_convert) {
    if(name=="Sexo"){
      data[,name] = as.numeric(data$Sexo)-1
    }
    data[,name] = as.numeric(as.character(data[,name]))
  } 
  return(data)
}