###############################################################################
####  Prepare LCMS mass hunter data for use in prediction_in_python.py ########
###############################################################################

import os, sys
import re
import cPickle as pickle
import numpy as np
import scipy as sp
import pandas as pd


def get_RP_batch1_crosswalk(clinDir):
    """Read in and clean data that will link reverse phase LCMS serum batch 1 to clinical data
    """
    #Note: Kristof's batch 1 data contains some cohort samples 
          #use info I've gathered (added to .csv file) on Study to make merge
    sampleMap = pd.read_csv(clinDir + "RP_batch1_sample_merge_info.csv", sep=',')
    #get patient ID and sample code to be compatible with other data
    sampleMap = sampleMap.assign(code = sampleMap['Code'].astype(str))
    sampleMap.code = "ID" + sampleMap.code.str.rjust(4,'0') #91 patients
    sampleMap = sampleMap.assign(Cod_Nin = sampleMap['Sample'].str.lower())
    #replace Cod_Nin with missings for hospital obs since clinical hospital data lacks them
        #will have merge problems later if not set to missing in both datasets
    sampleMap.loc[sampleMap.Study=="Hospital","Cod_Nin"] = np.nan
    return sampleMap

def main():
    
    inDir =  "/srv/scratch/ccotter/intermediate_data/"
    outDir = "/srv/scratch/ccotter/py_out/"
    clinDir = "/srv/scratch/ccotter/lab_and_clinical_data/Cleaned/" 

    #D1 is serum, D3 is saliva, D5 is urine (all normal phase)
    #RPD1 is reverse phase serum from first batch
    sample = 'RPD1_fill'
    if sample == 'D1':
        inFname = "comboD1_filter50n_wImpRF1.txt"
        outFname = "MassHuntNP" #name of pickled data frame to create
    elif sample == 'D3':
        inFname = "comboD3_filter50n_wImpRF1.txt"
        outFname = "SalivaMH" #name of pickled data frame to create
    elif sample == 'D5':
        inFname = "comboD5_filter50n_wImpRF1.txt"
        outFname = "UrineMH" #name of pickled data frame to create
    elif sample == 'RPD1_noFill':
        inFname = "091715_deng_QC_temp_nofillpeaks_copy.csv"
        outFname = "MassHuntRP_noFill"
    elif sample == 'RPD1_fill':
        inFname = "091715_deng_QC_temp_YESfillpeaks.csv" 
        outFname = "MassHuntRP_fill"
    elif sample == 'RPD1_isotop':
        inFname = "result_CAMERA_091715_deng_QC_temp_YESfillpeaks.csv"
        outFname = "MassHuntRP_isotop"
        
    # prepare data so it will properly merge with clinical data on code, Study, Cod_Nin

    if sample in ['D1','D3','D5']:
        df = pd.read_csv(inDir + inFname, sep='\t')
        print "Column names: " , list(df)
        print "Number of obs and of vars: " , df.shape
        # keep only LCMS variables and vars for merging
        keepervars = ["code","Study","Cod_Nin"]
        for var in df.columns.values:
            y = re.findall("MZ_.*", var)
            assert len(y) in (0,1)
            if len(y)==1: keepervars.append(y[0])
        print "Number of MZ values: " , len(keepervars) - 2
        df_final = df[keepervars]

    elif sample in ['RPD1_noFill', 'RPD1_fill', 'RPD1_isotope']:

        df = pd.read_csv(inDir + inFname, sep=',')

        ## Explore data ##
        #print "Column names: " , list(df) #284 columns -- one for each person-sample-run
        print "Number of obs and of vars: " , df.shape #rows are for mz, rt, npeaks values
        #are mz values currently unique? yes
        df['mz'].isnull().values.any() #no null values
        len(df['mz'].unique()) == len(df['mz']) #true

        ## Fix data ##
        #transpose so each row is person-sample-run and each column is unique mz value
        df2 = df.transpose()
        #change names of columns so they are called "MZ_" followed by mz number
        df2.columns = ["MZ_"+ str(x) for x in range(1,len(df2.columns)+1)] 
        #eliminate rows that were product of the transpose
        df2 = df2[11:]
        #add column indicating personID
        #Naming convention of dataframe row indices: 
            #RunOrder_PatientID_Unimportant_RepNumber. 
        sample = df2.index
        skip_start = [x[x.find('_')+1: ] for x in sample]
        df2['code'] = ['ID'+ (x[ :x.find('_')]).rjust(4,'0') for x in skip_start]
        uniqIDs = pd.unique(df2.code.ravel())
        print "uniqIDs: " , uniqIDs
        len(uniqIDs) #92 unique values
        #replace missing intensity values with zeros
        df2 = df2.where( (pd.notnull(df2)), 0 )
        #take mean value for each person and collapse to 1 row per person
            #Note: There are three technical reps for each patient.
        df2 = df2.groupby(df2.code).mean()
        
        #get crosswalk from patient codes to all necessary info for merging with clinical data
        sampleMap = get_RP_batch1_crosswalk(clinDir)

        #merge study info to LCMS data
        df_final = pd.merge(df2, sampleMap, left_index = True, right_on="code") #91 rows
        print "df_final: ", df_final.shape
        myinter = set(sampleMap.code).intersection(set(uniqIDs))
        test = set(uniqIDs) - myinter #forget about ID341. -- already ahve ID0341 in data
        print "test: " , myinter, test
      
    # save for later use (employ pd.read_pickle() to load)
    print "Creating pickled data frame"
    df_final.to_pickle(os.path.join(outDir, outFname)) 


if __name__ == '__main__':
    main()
