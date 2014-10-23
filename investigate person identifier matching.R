

#select main directory in which to find subfolders
homeDir = "/srv/scratch/carolyn/" #on Amold from Amold's perspective
#homeDir = "/home/carolyn/dengue_dx/" #on my PC

clinical_inputsDir = paste(homeDir, "lab_and_clinical_data/Cleaned/", sep="")
lcmsCO_inputsDir = paste(homeDir, "lcms_data_processed_in_CO/Cleaned/", sep="")
outputsDir = paste(homeDir, "intermediate_data/", sep="")
resultsDir = paste(homeDir, "Results/", sep="")

source(paste(codeDir, "clean_data_functions.R", sep=""))


### see who has an LCMS .d file but is not in Natalia's LCMS text file (and vice versa) ###

#Natalia's LCMS data -- 88 observations
respD1_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_first batch.txt", lcms_run=1, printme=TRUE)
#.d files that I was able to process with msconvert and python -- 89 observations
#load(paste(outputsDir,"df1_from_python_withRdata.RData", sep="")) #will load comboD1_filter50n into R
load(paste(outputsDir,"df_from_python.RData", sep="")) #will load respD into R -- codes 0086 and 1208 in here but not in CO's LCMS
clinical_hospitD[which(clinical_hospitD$code==86 | clinical_hospitD$code==1208 ),c("Res_Final","ClasificacionFinal")] #these codes are in clinical/lab data, too
get_mismatches(respD1_filter50n$code, respD$X1, "CO's D1", "my D1")

#Natalia's LCMS data -- 83 observations
respD2_filter50n = clean_LCMS(infile="LCMS_serum_Nica_50percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) 
#LCMS .d files that I was able to process with msconvert and python -- 68 observations (why so few?)
load(paste(outputsDir,"df2_from_python_withRdata.RData", sep="")) 
  #loads comboD1_filter50n - missing "ID1503" "ID4015" "ID4592" "ID5666" "ID6680" "ID7194" "ID7744", explained by indeterminant Dx
load(paste(outputsDir,"df2_from_python.RData", sep="")) 
  #loads respD (75 obs) - missing "ID0763" "ID1518" "ID1519" "ID3635" "ID4499" "ID5349" "ID5370" "ID6453" (try converting .d files again)
get_mismatches(comboD1_filter50n$code, respD$X1, "with R", "before merge")
get_mismatches(respD2_filter50n$code, respD$X1, "CO's D2", "my D2")
respD2_filter50n[which(respD2_filter50nD$code==5349 | respD2_filter50n$code==5370 ),c("Res_Final","ClasificacionFinal")]

### function to read in the LC-MS text file and organize it into dataframe ###
fake_clean_LCMS = function(infile, lcms_run, roundme=FALSE, decimals=2, printme=FALSE) {
  
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
    #respD = get_study(respD)
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

respD2_filter50n = fake_clean_LCMS(infile="MS abundances_50percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples
respD2_filter10n = fake_clean_LCMS(infile="MS abundances_10percent_batches 3 and 4.txt", lcms_run=2, printme=TRUE) #83 samples
table(table(respD2_filter10n$code))
table(table(respD2_filter50n$code))

# process file that contains indicator of study (cohort versus hospital)
study_info = read.delim(paste(inputsDir,"Study classifications_batches 3 and 4.txt", sep=""), header=TRUE, nrows=200)
newnames = as.character(c("Study","code","Res_Final2","OMS"))
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
fcode =    function(x) sprintf("ID%04d",as.numeric(as.character(x["code"])))
code_list = apply(X=study_info, MARGIN = 1, FUN=fcode)
study_info["code"] = code_list
table(table(study_info$code))

#merge with lc-ms data
respD_new10 = merge(study_info, respD2_filter10n, by="code", all=TRUE)

respD_new50 = merge(study_info[,c("Study","code")], respD2_filter50n, by="code", all=TRUE)
