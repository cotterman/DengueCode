###############################################################################
############ Run prediction_in_python (loop thru parameter settings) ##########
###############################################################################

import subprocess

#Set parameters
    #default for all boolean flags is False (ex: '' means False)
    #specifying flag turns indicators to true (ex: '--run_MainAnalysis')
run_MainAnalysis = '--run_MainAnalysis'
plot_MainAnalysis = '' #'--plot_MainAnalysis' 
NoInitialDHF = '--NoInitialDHF' #should generally be true (for prediction)
onlyLCMSpatients = '--onlyLCMSpatients'
#Loop thru parameters options
ps = []
for outcome in ['is.DHF_DSS']:
    if outcome=='is.DHF_DSS':
        NoOFI_list = ['--NoOFI',''] #2 options for is.DHF_DSS
    else:
        NoOFI_list = [''] #the only sensible value for is.DEN analysis
    for NoOFI in NoOFI_list:

        for inLCMSData in ["RPbins50x50", "NPbins50x50"]:

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

                    #string args together
                    args = " ".join([run_MainAnalysis, plot_MainAnalysis, 
                            include_clinvars, include_LCMSvars, onlyLCMSpatients,
                            NoOFI, NoInitialDHF, include_imp_dums, imp_dums_only,
                            '--outcome', outcome, '--inLCMSData', inLCMSData,
                            '--predictor_desc', predictor_desc])

                    #begin each python instance
                        #Note: processes will be run in parallel
                        #beware of overloading machine with too many processes
                            #ex: n_jobs (in grid search) also spawns jobs 
                    p = subprocess.Popen("python prediction_in_python.py " + args, shell=True)
                    ps.append(p)

#wait till each process is finished before exiting
for p in ps: p.wait()
    


