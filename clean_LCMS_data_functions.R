### function to create a version of the data which replaces values for missing compounds with zeros
  #other options: fill with min value observed for that compound,
  #fill with 1/2 of the min value, input value
fill_LCMS_blanks = function(mydata){
  cindices = grep("MZ_",colnames(mydata))
  processedD = mydata
  #note: must have a return statement b/c assignment operator will not return anything (and default is 0)
  FillMissingsWithZeros <- function(x){ x[which(is.na(x))] <- 0; return(x)}
  processedD[cindices] = lapply(mydata[cindices],FUN=FillMissingsWithZeros)
  return(processedD)
}

### function to read in the LC-MS text file and organize it into dataframe ###
clean_LCMS = function(infile, lcms_run, roundme=FALSE, decimals=2, printme=FALSE) {
  
  #infile="LCMS_saliva_Nica_50percent.txt"
  #infile = "LCMS_urine_Nica_50percent.txt"
  #lcms_run=5
  #printme=TRUE
  #decimals=2
  #roundme=FALSE
  
  abundancesD1 = read.delim(paste(lcmsCO_inputsDir, infile, sep=""), header=TRUE, nrows=200)
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
  
  #add an indicator for type of study and sample code (will need this later when we merge)
  if(lcms_run==1){
    respD["Study"] = "Hospital"
  } else {
    respD = get_study(respD, lcms_run, clinical_inputsDir)
  } 
  
  #count  the number of duplicated MZ values
  if(printme==TRUE){
    print("Number of duplicated MZ values with no rounding of MZ values:")
    print(table(table(MZ_Nums_reshaped1))) #without rounding
    print(paste("Number of duplicated MZ values with rounding to", decimals, "decimals"))
    print(table(table(MZ_Nums_short))) #with rounding
  }
  
  #create a version of the data which replaces values for missing compounds with zeros
  respD = fill_LCMS_blanks(respD)
  
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

### function gets list of patient IDs that are present in specified data
get_IDs_in_common = function(resp_D1D2, respD3_filter50n, respD5_filter50n){
  IDs_in_resp_D1_D2 = resp_D1D2[,c("code","Study","LCMS_run")]
  temp = merge(IDs_in_resp_D1_D2[,c("code","Study","LCMS_run")], 
               respD3_filter50n[,c("code","Study","LCMS_run")], by=c("code","Study"), all=T)
  temp$serum = temp$LCMS_run.x
  temp$LCMS_run.x = NULL
  temp$saliva = temp$LCMS_run.y
  temp$LCMS_run.y = NULL
  IDs_in_resp_all = merge(temp, respD5_filter50n[,c("code","Study","LCMS_run")], by=c("code","Study"), all=T)
  IDs_in_resp_all$urine = IDs_in_resp_all$LCMS_run
  IDs_in_resp_all$LCMS_run = NULL
  #are all saliva and urine samples are from same people? yes
  table(IDs_in_resp_all$urine, IDs_in_resp_all$saliva, useNA="ifany")
  #are all saliva and serum samples from the same people? just 1 (might be mistake -- ID1408 in cohort, which was already problematic)
  table(IDs_in_resp_all$serum, IDs_in_resp_all$saliva, useNA="ifany")
  return(IDs_in_resp_all)
}