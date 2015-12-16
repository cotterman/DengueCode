###############################################################################
############ Run prediction_in_python (loop thru parameter settings) ##########
###############################################################################

import subprocess

#Set parameters
    #default for all boolean flags is False (ex: '' means False)
    #specifying flag turns indicators to true (ex: '--run_MainAnalysis')
run_MainAnalysis = '--run_MainAnalysis'
run_SL = '--run_SL'
plot_MainAnalysis = '' #'--plot_MainAnalysis' 
run_BestSubset = '' #'--run_BestSubset'
run_FDR = '' #'--run_FDR'
NoInitialDHF = '' #should generally be true - applies only to is.DHF_DSS analysis
onlyLCMSpatients = '--onlyLCMSpatients'
testlib = '' #leave as '' to run with full SL library 

#Loop thru parameters options
ps = []
for outcome in ['is.DEN']:
    if outcome=='is.DHF_DSS':
        NoOFI_list = [''] #options for is.DHF_DSS: '--NoOFI',''
    else:
        NoOFI_list = [''] #the only sensible value for is.DEN analysis
    for NoOFI in NoOFI_list:              
        
  #"NPbins50x50","MassHuntNP","RPbins50x50","MassHuntRP_fill","MassHuntRP_noFill","MassHuntRP_isotope"
        for inLCMSData in ["MassHuntNP","MassHuntRP_noFill"]:
        #for inLCMSData in ["SalivaMH", "UrineMH"]:

            predictor_desc  = "covarlist_all"
            #for predictor_desc in ["covarlist_all","covarlist_noUltraX",
            #        "covarlist_CohortRestrict","covarlist_genOnly"]:

            for include_clinvars in ['','--include_clinvars']:
                include_imp_dums = ''
                #implement 3 versions: no imp_dums, with imp_dums, and imp_dums_only
                #for include_imp_dums in ['', '--include_imp_dums']:
                if include_imp_dums=='--include_imp_dums':
                    imp_dums_only_list = ['','--imp_dums_only']
                else:
                    imp_dums_only_list = ['']
                imp_dums_only = ''

                for include_LCMSvars in ['','--include_LCMSvars']:
                #for imp_dums_only in imp_dums_only_list:

                    for correctData in ['0','1','2']:

                        #string args together
                        args = " ".join([run_MainAnalysis, run_SL, plot_MainAnalysis, 
                                run_BestSubset, run_FDR,
                                include_clinvars, include_LCMSvars, onlyLCMSpatients,
                                NoOFI, NoInitialDHF, include_imp_dums, imp_dums_only,
                                '--outcome', outcome, '--inLCMSData', inLCMSData, 
                                '--correctData', correctData, '--testlib', testlib,
                                '--predictor_desc', predictor_desc])

                        #begin each python instance
                            #Note: processes will be run in parallel
                            #beware of overloading machine with too many processes
                                #ex: n_jobs (in grid search) also spawns jobs 
                        #does not make sense to include no predictors
                        if include_clinvars == '' and include_LCMSvars == '':
                            continue
                        else:
                            p = subprocess.Popen("python prediction_in_python.py " + args, shell=True)
                            ps.append(p)

#wait till each process is finished before exiting
for p in ps: p.wait()
    


