
#eliminate variables that are missing for more than half of observations (n=1624) -- todo
#temp_vars_missCount = sapply(clinical_full_clean[,as.character(temp$variable.name.in.final.data)], function(x) sum(is.na(x)))
#temp_vars_missCount
#cvars_fewMiss = temp[temp_vars_missCount<=812] #40 
#clinicCount = length(cvars_fewMiss)
#cat("\n Clinical variables with no missings \n")
#print(cvars_fewMiss) #view 
#cat("\n Clinical variables with >= 812 missings \n")
#print(temp[temp_vars_missCount>812]) #view 
#temp = cvars_fewMiss