#Run with R 3.1.1
#To run remotely using Rstudio, go to URL: amold.lbl.gov:8787/

################################################################################
####### Prepare data from MZmine so it works with R prediction functions #######
################################################################################

#### Establish directories ####

#directory containing code
codeDir = "/srv/scratch/carolyn/Dengue_code/" #on Amold from Amold's perspective
#codeDir = "/home/carolyn/temp_Dengue_code/" #on myPC

#select main directory in which to find subfolders
homeDir = "/srv/scratch/carolyn/" #on Amold from Amold's perspective
#homeDir = "/home/carolyn/dengue_dx/" #on my PC

clinical_inputsDir = paste(homeDir, "lab_and_clinical_data/Cleaned/", sep="")
lcmsMZmine_inputsDir = paste(homeDir, "output_files_MZmine/", sep="")
outputsDir = paste(homeDir, "intermediate_data/", sep="")
resultsDir = paste(homeDir, "Results/", sep="")


#read in csv output created by running MZmine batch file
mzmineD1 = read.csv(paste(lcmsMZmine_inputsDir, "aligned_peaks_t1.csv", sep=""), header=TRUE, nrows=200)

#Keep only info on peak area (in future will want to do with peak height as well) and rowID
  #Will name the columns according to the rowID (can lookup corresponding RT and MZ values as desired)
    #Using row ID rather than MZ values to avoid one-to-many mapping when MZ value is not unique
temp = colnames(mzmineD1) #all column names
ID_cols = temp[grep(".mzML.filtered.filtered.peak.area",temp)] 
mzmineSmall = mzmineD1[ , c("row.ID", ID_cols)]

#Rename columns so that they are IDXXXX where XXXX is 4-digit patient code
replace_blah_with_ID = function(x) substr( sub("Nicaserhilic", "ID", x), 1, 6)
ID_Names = sapply(X=ID_cols, replace_blah_with_ID )
ID_Names
colnames(mzmineSmall) = c("rowID",ID_Names)

#reshape data so each row represents a different patient while columns indicate abundance info
mzminet1 = as.data.frame(t(mzmineSmall))[-1,]
replace_V_with_MZ = function(x) sub("V", "MZ_", x)
MZ_Names = sapply(X=colnames(mzminet1), replace_V_with_MZ)
colnames(mzminet1) = MZ_Names
mzminet1$code = rownames(mzminet1)
mzmine_clean = mzminet1[c("code",MZ_Names)] #reorder columns

  
#### now can use same code as for prediction with LCMS from python ####


