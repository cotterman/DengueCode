
### Basic summary of clinical and lab variables ###

#### messing around with colors ###

# These are the "Tableau 20" colors as RGB.  
tableau20 = cbind(rgb(31/255, 119/255, 180/255), rgb(174/255, 199/255, 232/255), rgb(255/255, 127/255, 14/255), rgb(255/255, 187/255, 120/255),  
             rgb(44/255, 160/255, 44/255), rgb(152/255, 223/255, 138/255), rgb(214/255, 39/255, 40/255), rgb(255/255, 152/255, 150/255),  
             rgb(148/255, 103/255, 189/255), rgb(197/255, 176/255, 213/255), rgb(140/255, 86/255, 75/255), rgb(196/255, 156/255, 148/255),  
             rgb(227/255, 119/255, 194/255), rgb(247/255, 182/255, 210/255), rgb(127/255, 127/255, 127/255), rgb(199/255, 199/255, 199/255),  
             rgb(188/255, 189/255, 34/255), rgb(219/255, 219/255, 141/255), rgb(23/255, 190/255, 207/255), rgb(158/255, 218/255, 229/255)) 
n = 20
par(mfrow=c(1,1))
pie(rep(1,n), col=tableau20, main="tableau20") #view color wheel


###############################################################################
########### Graphs and tables that use hospital and cohort data ###############
###############################################################################

d2 = clin12_full_clean
dh = clin12_full_clean[which(clin12_full_clean$Study=="Hospital"),]
dc = clin12_full_clean[which(clin12_full_clean$Study=="Cohort"),]

#combine 7 and more days into 1 group
dh[,"plot_DaysSick"] = dh$DaysSick
dh[which(dh$plot_DaysSick>7),"plot_DaysSick"] = 7
dc[,"plot_DaysSick"] = dc$DaysSick
dc[which(dc$plot_DaysSick>7),"plot_DaysSick"] = 7

#plot DaysSick for cohort and hospit on same graph
png(paste(resultsDir,"DaysSick_CandH.png", sep=""), width=800, height=500)
par(cex=1.6)
mybreaks = c(seq(0,7,1))
#round the value labels to hundreths place
dh_labs = as.character(sprintf("%.2f",hist(dh$plot_DaysSick, breaks=mybreaks)$density))
dc_labs = as.character(sprintf("%.2f",hist(dc$plot_DaysSick, breaks=mybreaks)$density))
hist(dh$plot_DaysSick, col=rgb(0,0,1,.5), breaks=mybreaks, freq=F, axes=F, labels=dh_labs,
    ylim=c(0,.5), ylab="Proportion", xlab="Days since fever onset", main="")
hist(dc$plot_DaysSick, col=rgb(1,0,0,.25), breaks=mybreaks, freq=F, labels=dc_labs, add=T)
axis(side=2)
axis(side=1, at=c(seq(.5,6.5,1)), labels=c(seq(1,6,1),"7+"))
legend(x="topright",legend=c("Hospital","Cohort"), col=c(rgb(0,0,1,.5),rgb(1,0,0,.25)), 
       pch=19, bty="n")
dev.off()

#table of initial and final diagnoses
table(dh$WHO_initial_given, dh$WHOFinal4cat, useNA="always")  #initial Dx calculated by Douglas

#plot of proportion showing severe symptoms upon arrival, by day of illness
dh_severe = dh[which(dh$WHOFinal3cat=="DHF_DSS"),] #240 with eventual severe dengue
dh_severe$initialSevere = as.integer(dh_severe$WHO_initial_given=="DHF" |
                             dh_severe$WHO_initial_given=="DSS")
collapsed = summaryBy(initialSevere ~ plot_DaysSick, FUN=mean, data=dh_severe)
png(paste(resultsDir,"proportion_initialSevere.png", sep=""), width=800, height=500)
par(mar=c(5,6,4,2), cex=1.6)
bp = barplot(collapsed$initialSevere.mean,
        ylab="Patients displaying severe symptoms \n upon arrival (proportion)", xlab="Days since fever onset",
        col=rgb(0,0,1,.5), names.arg=c(seq(1,6,1),"7+"))
#add value labels
text(x=bp, y=collapsed$initialSevere.mean+.05, sprintf("%.2f",collapsed$initialSevere.mean))
dev.off()


###############################################################################
########### Numbers and graphs to make for each dataset sepeartely ############
###############################################################################

## chose which data to use here ##
#all
df = clin12_full_clean
#hospital only 
#df = clin12_full_clean[which(clin12_full_clean$Study=="Hospital"),]
#cohort only 
#df = clin12_full_clean[which(clin12_full_clean$Study=="Cohort"),]

## Summarize lab and clinical data ##
#summarize_clinical(df)


###############################################################################
###### Data description section of Dissertation Ch2 -- sample sizes ##########
###############################################################################

#number of children
length(unique(df$code)) #3578

#time period
min(df$FTM) #"2004-09-01"
max(df$FTM) #"2015-03-13"

#ages
summary(df$age) #6 months = 15.89
sum(df$age<.5) #none

table(df$is.DEN)

print("WHOFinal4cat")
table(df$WHOFinal4cat)
table(df[-which(df$serum == 1),"WHOFinal4cat"])

print("DaysSick")
table(df$DaysSick)
table(df[-which(df$serum == 1),"DaysSick"])

print("ClasificacionPrimerDia by WHOFinal4cat")
print(table(clinical_D1_clean$ClasificacionPrimerDia, clinical_D1_clean$WHOFinal4cat))
table(df[-which(df$serum == 1),]$ClasificacionPrimerDia, df[-which(df$serum == 1),]$WHOFinal4cat)

print("Initial Dx of severe dengue patients, by DaysSick")
table(df[which(df$WHOFinal4cat == "DHF" | df$WHOFinal4cat == "DSS"),]$DaysSick, 
      df[which(df$WHOFinal4cat == "DHF" | df$WHOFinal4cat == "DSS"),]$ClasificacionPrimerDia)

print("Number of patients with no missing clinical information (restricting to vars collected at both hospital and clinic settings)")
num_miss_by_patient = apply(df[,c(covarlist_CohortRestrict)], MARGIN=1, FUN=function(x) sum(is.na(x)))
table(num_miss_by_patient)


###############################################################################
############### Summary of clinical variables for dissertation ################
###############################################################################

#create binary indicators for the general variables that were not binary
#XD_binary = create_binary_variables(df) #function created in "clean_data_functions.R"

## get list of clinical variables to include (returns list of variable.name.in.final.data)
varlist = get_clinic_var_list(clinic_varsD, outcome="either", 
                                        eliminate_vars_with_missings=F, eliminate_constant_vars=T, 
                                        eliminate_vars_with_minXnomiss=50,
                                        XD=df, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T)

## overview of these variables
contents(df[,varlist])

#compare the serotypes for hospital v cohort patients
th = table(dh[which(dh$is.DEN==1), "PCR"])
tc = table(dc[which(dc$is.DEN==1), "PCR"])
compare = as.table(rbind(th[2:4], tc[2:4]))
chisq.test(compare) #significantly different -- more type 3s in hospital

## subset the data with variable info to include only rows that we care about
varinfo = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% varlist==T),]

## Output table for including in dissertation -- description of variables
outTable = varinfo[order(varinfo$CC_broadcat_sort, varinfo$CC_cat_sort, varinfo$CC_name), ] #sort 
outTable = outTable[,c("CC_name","CC_category","CC_type","CC_description","in.cohort.data")] #restrict
#open this file in excel and copy/paste into lyx table using middle-click mouse button
write.csv(outTable, paste(resultsDir,"variable_descriptions.csv", sep=""), row.names = FALSE)


## view variable absence in cohort data
table(varinfo[which(varinfo$in.cohort.data==0),]$Variable.Category, varinfo[which(varinfo$in.cohort.data==0),]$CC_type)
sum(varinfo$in.cohort.data==0, na.rm=T) #38 (so 47 present)
## view variable categories and types (binary, continuous, etc.)
    #in displaying output, put "Laboratory-General", "Laboratory-Clinical lab" and "Laboratory-Virological" with "Laboratory-Clinical lab"
    #note that I consider "daysSick" to be "General - symptom" while it had been in "Demographics/Gen Info" group
table(varinfo$Variable.Category, varinfo$CC_type)
table(varinfo$CC_broadcat, varinfo$CC_type) #less detailed
table(varinfo$CC_category, varinfo$CC_type) #slightly more detailed than what I present in qual, but same basic groupings
#of the variables in varlist, list those which cannot be used to predict DHF/DSS  
varinfo[which(varinfo$Use.in.DF.vs.DHF.DSS.prediction!=1),c("variable.name.in.final.data","Variable.Category")]
#of the variables in varlist, list those which cannot be used to predict ND vs DENV  
varinfo[which(varinfo$Use.in.ND.vs.DEN.prediction!=1),c("variable.name.in.final.data","Variable.Category")]




###############################################################################
############## Graphs for dissertation appendix ###############################
###############################################################################
df = dh #toggle this

## Demographics/ general - age and gender
png(paste(resultsDir,"Age_and_gender_hospit.png", sep=""), width=800, height=500) #pngs are more portable but don't scale (else use svg)
par(mfrow=c(1,1), cex=1.85)
max_age = max(trunc(df[,"age"]))
mybreaks = c(seq(-1,max_age,1))
p1 = hist(trunc(df[which(df$Sexo=="female"),"age"]), breaks = mybreaks, freq=T)
p2 = hist(trunc(df[which(df$Sexo=="male"),"age"]), breaks = mybreaks, freq=T)
plot(p1, xlab = "Age", main="", xlim=c(-.5,max_age), ylim=c(0,100), col=rgb(.5,.5,.5,.65), axes=F)
plot(p2, col=rgb(0,0,1,.5), add=T, axes=F)
axis(side=2)
axis(side=1, at=c(seq(-.5,max_age-.5,1)), labels=c(seq(0,max_age,1)))
legend(x="topright",legend=c("Female","Male"), col=c(rgb(.5,.5,.5,.65),rgb(0,0,1,.5)), 
       pch=19, bty="n")
dev.off()


## General symptoms, binary and categorical 
  #todo: add the remaining categorical variable -- is.pulse_rapid (for Pulse)

#remove Epigastralgia since it was not included in latest hospital data installment :()
var_names = c(as.character(varinfo[which(varinfo$Variable.Category=="Gen Sign/Symptom" &
                                              varinfo$CC_type=="binary" & 
                                              varinfo$variable.name.in.final.data != "Epigastralgia"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as per)
get_mean = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals = sapply(df[c(var_names)], FUN=get_mean)
#add is.pulse_rapid (based on categorical variable pulso)
mean_val = sapply(df[c("is.pulse_rapid")], FUN=mean, na.rm=T)
mean_vals = c(mean_val, mean_vals)
#create better names to display (CC_name)
myd = data.frame(mean_vals, names(mean_vals))
namefinder = merge(myd, varinfo[,c("CC_name","variable.name.in.final.data")], 
                   by.x="names.mean_vals.", by.y="variable.name.in.final.data")
addon = myd[which(myd$names.mean_vals.=="is.pulse_rapid"),]
addon$CC_name = "Rapid pulse"
addon = rbind(namefinder, addon)
toplot = addon[order(addon$mean_vals),]
#graph
png(paste(resultsDir,"General_symtoms_binary_hospit.png", sep=""), width=800, height=800)
par(mfrow=c(1,1),las=2,mar=c(5,10,4,2), cex=1.9) # make label text perpendicular to axis # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(toplot$mean_vals, horiz=T, names.arg=toplot$CC_name, xlim=c(0,1), xlab="Proportion affirmative", 
        beside=T, xaxt="n", col=rgb(0, 0, 1,.5))
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=toplot$mean_vals+0.08, y=bplt, labels = sprintf("%.2f",toplot$mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

## General symptoms, continuous
var_names = c(as.character(varinfo[which(varinfo$Variable.Category=="Gen Sign/Symptom" &
                                        varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"General_symtoms_continuous_hospit.png", sep=""), width=800, height=500)
par(mfrow=c(2,3), mar=c(5, 5, 4, 2) + 0.1, mai=c(1,.3,.1,.2), cex.lab=2, cex.axis=1.8) 
for(var in var_names){
  myname = varinfo[which(varinfo$variable.name.in.final.data==var),"CC_name"]
  hist(df[,var], xlab=myname, col=rgb(0, 0, 1,.5), main="", las=0, ylab="", freq=F, )
}
dev.off()


## General indicators of hemorrhaging, binary and categorical
var_names1 = c(as.character(varinfo[which(varinfo$Variable.Category=="Hem Sign/Symptom" &
                                           varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as percent)
get_mean1 = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals1 = sapply(df[var_names1], FUN=get_mean1) #previously used with XD_binary
#the only gen-hem variable that is not already binary is the torniquete one so I created binary form of it
var_names2 = c("is.torniquete10plus","is.torniquete20plus")
get_mean2 = function(x) mean(as.numeric(x), na.rm=T)
mean_vals2 = sapply(df[var_names2], FUN=get_mean2) #previously used with XD_binary
mean_vals = c(mean_vals1, mean_vals2)
var_names = c(var_names1, var_names2)
#rename and sort
myd = data.frame(mean_vals, names(mean_vals))
namefinder = merge(myd, varinfo[,c("CC_name","variable.name.in.final.data")], 
                   by.x="names.mean_vals.", by.y="variable.name.in.final.data")
#add rows to varinfo so we keep the binary torniquete vars
toadd = myd[which(myd$names.mean_vals.=="is.torniquete10plus" |
                           myd$names.mean_vals.=="is.torniquete20plus"),]
toadd$CC_name = toadd$names.mean_vals.
toplot_prelim = rbind(namefinder, toadd)
toplot = toplot_prelim[order(toplot_prelim$mean_vals),]
#graph
plot.new()
png(paste(resultsDir,"General_hemorrhage_binary_hospit.png", sep=""), width=800, height=800)
par(mfrow=c(1,1),las=2, cex=1.9) # make label text perpendicular to axis
par(mar=c(5,10,4,2)) # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(toplot$mean_vals, horiz=T, names.arg=toplot$CC_name, xlim=c(0,1), 
               xlab="Proportion affirmative", 
        cex.names=1, beside=T, xaxt="n", col=rgb(0, 0, 1,.5))
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=toplot$mean_vals+0.08, y=bplt, labels = sprintf("%.3f",toplot$mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

## Blood lab count variables, continuous
plot.new()
var_names = c(as.character(varinfo[which( (varinfo$CC_category=="Blood - count" | varinfo$CC_category=="Urine - count") &
                                            varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"Blood_counts_continuous_hospit.png", sep=""), width=880, height=550)
par(mfrow=c(2,4))
for(var in var_names){
  myname = varinfo[which(varinfo$variable.name.in.final.data==var),"CC_name"]
  hist(df[,var], xlab=myname, cex.lab=2, cex.axis=1.7, col=rgb(0, 0, 1,.5), main="", las=0, ylab="", freq=F, )
}
dev.off()
## Blood lab chemistry variables, continuous
plot.new()
var_names = c(as.character(varinfo[which(varinfo$CC_category=="Blood - chemistry" &
                                           varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"Blood_chemistry_continuous_part1.png", sep=""), width=800, height=500)
par(mfrow=c(2,4))
for(var in var_names[1:8]){
  myname = varinfo[which(varinfo$variable.name.in.final.data==var),"CC_name"]
  hist(df[,var], xlab=myname, cex.lab=2, cex.axis=1.7, col=rgb(0, 0, 1,.5), main="", las=0, ylab="", freq=F, )
}
dev.off()
png(paste(resultsDir,"Blood_chemistry_continuous_part2.png", sep=""), width=880, height=550)
par(mfrow=c(3,3))
for(var in var_names[9:17]){
  myname = varinfo[which(varinfo$variable.name.in.final.data==var),"CC_name"]
  hist(df[,var], xlab=myname, cex.lab=2, cex.axis=1.7, col=rgb(0, 0, 1,.5), main="", las=0, ylab="", freq=F, )
}
dev.off()

## Ultrasound, continuous
plot.new()
var_names = c(as.character(varinfo[which((varinfo$CC_category=="Ultrasound"|varinfo$CC_category=="X-ray") &
                                           varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"ultrasound_continuous.png", sep=""), width=800, height=300)
par(mfrow=c(1,3))
for(var in var_names){
  myname = varinfo[which(varinfo$variable.name.in.final.data==var),"CC_name"]
  hist(df[,var], xlab=myname, cex.lab=2, cex.axis=1.7, col=rgb(0, 0, 1,.5), main="", las=0, ylab="", freq=F, )
}
dev.off()
### Ultrasound and xray, binary
plot.new()
var_names = c(as.character(varinfo[which((varinfo$CC_category=="Ultrasound"|varinfo$CC_category=="X-ray") &
                                           varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as percent)
get_mean = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals = sapply(df[var_names], FUN=get_mean)
#create better names to display (CC_name)
myd = data.frame(mean_vals, names(mean_vals))
namefinder = merge(myd, varinfo[,c("CC_name","variable.name.in.final.data")], 
                   by.x="names.mean_vals.", by.y="variable.name.in.final.data")
toplot = namefinder[order(namefinder$mean_vals),]
png(paste(resultsDir,"ultrasoundXray_binary.png", sep=""), width=800, height=450)
par(mfrow=c(1,1), las=2) # make label text perpendicular to axis
par(mar=c(5,12,4,2), cex=1.5) # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(toplot$mean_vals, horiz=T, names.arg=toplot$CC_name, xlim=c(0,1), xlab="Proportion affirmative", 
        cex.names=1, beside=T, xaxt="n", col=rgb(0, 0, 1,.5))
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=toplot$mean_vals+0.08, y=bplt, labels = sprintf("%.3f",toplot$mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

###############################################################################
############## Missingness explorations for the D #############################
###############################################################################

## missingness
# "38 variables never missing for any of the hospital patients while 10 variables are never missing across the combined data"
df = dc #toggle this
var_missCount = sapply(df[varlist], function(x) sum(is.na(x)))
table(var_missCount)
var_nomissCount = sapply(df[varlist], function(x) sum(!is.na(x)))
table(var_nomissCount)
#variables with no missings
var_noMiss = c(as.character(varlist[var_missCount==0])) 
length(var_noMiss) 
varinfo[which(varinfo$variable.name.in.final.data %in% var_noMiss==T), c("CC_name","variable.name.in.final.data","Variable.Category")]
#variables with some missings 
var_someMiss = c(as.character(varlist[var_missCount>0]))
length(var_someMiss)

var_someMiss_info = varinfo[which(varinfo$variable.name.in.final.data %in% var_someMiss==T), c("variable.name.in.final.data","Variable.Category")]
var_someMiss_gen = c(as.character(var_someMiss_info[which(var_someMiss_info$Variable.Category=="Demographics/Gen Info" | 
                                                            var_someMiss_info$Variable.Category=="Gen Sign/Symptom" | 
                                                            var_someMiss_info$Variable.Category=="Laboratory-General" |
                                                            var_someMiss_info$Variable.Category=="Hem Sign/Symptom"),c("variable.name.in.final.data")]))
var_someMiss_lab = c(as.character(var_someMiss_info[which(var_someMiss_info$Variable.Category=="Laboratory-Clinical lab" | 
                                                            var_someMiss_info$Variable.Category=="Laboratory-Virological" |
                                                            var_someMiss_info$Variable.Category=="Laboratory-Urine analysis"),c("variable.name.in.final.data")]))
var_someMiss_ultraX = c(as.character(var_someMiss_info[which(var_someMiss_info$Variable.Category=="Laboratory-Ultrasound" | 
                                                               var_someMiss_info$Variable.Category=="Laboratory-Ultrasound and X-ray" | 
                                                               var_someMiss_info$Variable.Category=="Laboratory-X ray"),c("variable.name.in.final.data")]))
table(df[which(df$WHOFinal4cat!="ND"),"IR"], useNA=c("always")) #missing 37/964
table(df[which(df$WHOFinal4cat!="ND"),"PCR"], useNA=c("always")) #missing 131/964

## graphs of number of missing values per variable, organized by variable category ##

#general cat
named_vals = as.data.frame(var_missCount[var_someMiss_gen])
named_vals$variable.name.in.final.data = rownames(named_vals)
named_vals$val = var_missCount[var_someMiss_gen]
namefinder = merge(named_vals, varinfo[,c("CC_name","variable.name.in.final.data")], by="variable.name.in.final.data")
toplot = namefinder[order(namefinder$val),]
#plot
png(paste(resultsDir,"Clinical_vars_missingness_gen_hospit.png", sep=""), width=800, height=500)
par(mfrow=c(1,1), las=1, cex=2, mar=c(5, 10, 4, 2) + 0.1) # make label text perpendicular to axis
barplot(toplot$val, names.arg=toplot$CC_name, horiz=T, beside=T, col=rgb(44/255, 160/255, 44/255,.5), main="", xlim=c(0,nrow(df)))
dev.off()

#blood and urine lab
named_vals = as.data.frame(var_missCount[var_someMiss_lab])
named_vals$variable.name.in.final.data = rownames(named_vals)
named_vals$val = var_missCount[var_someMiss_lab]
namefinder = merge(named_vals, varinfo[,c("CC_name","variable.name.in.final.data")], by="variable.name.in.final.data")
toplot = namefinder[order(namefinder$val),]
#plot
png(paste(resultsDir,"Clinical_vars_missingness_bloodlab_hospit.png", sep=""), width=800, height=1000)
par(mfrow=c(1,1), las=1, cex=2.2, mar=c(5, 10, 4, 2) + 0.1) # make label text perpendicular to axis
barplot(toplot$val, names.arg=toplot$CC_name, horiz=T, beside=T, col=rgb(44/255, 160/255, 44/255,.5), main="", xlim=c(0,nrow(df)))
dev.off()

#ultra sound and x-ray
named_vals = as.data.frame(var_missCount[var_someMiss_ultraX])
named_vals$variable.name.in.final.data = rownames(named_vals)
named_vals$val = var_missCount[var_someMiss_ultraX]
namefinder = merge(named_vals, varinfo[,c("CC_name","variable.name.in.final.data")], by="variable.name.in.final.data")
toplot = namefinder[order(namefinder$val),]
#plot
png(paste(resultsDir,"Clinical_vars_missingness_ultraX_hospit.png", sep=""), width=800, height=500)
par(mfrow=c(1,1), las=1, cex=2.2, mar=c(5, 10, 4, 2) + 0.1) # make label text perpendicular to axis
barplot(toplot$val, names.arg=toplot$CC_name, horiz=T, beside=T, col=rgb(44/255, 160/255, 44/255,.5), main="", xlim=c(0,nrow(df)))
dev.off()


#Variables available per patient
#hospit
patient_nomissCounth = apply(dh[varlist], MARGIN=1, function(x) sum(!is.na(x)))
playh = table(patient_nomissCounth)
var_counts_decreasingh = sort(names(playh), decreasing=TRUE)
cumsum(playh[var_counts_decreasingh])
#cohort
patient_nomissCountc = apply(dc[varlist], MARGIN=1, function(x) sum(!is.na(x)))
playc = table(patient_nomissCountc)
var_counts_decreasingc = sort(names(playc), decreasing=TRUE)
cumsum(playc[var_counts_decreasingc])
#graph
png(paste(resultsDir,"Clinical_vars_perPatient_CandH.png", sep=""), width=800, height=500)
par(mar=c(5,6,4,2)+.1, cex=1.6) #default is  c(5, 4, 4, 2) + 0.1
plot(y=cumsum(play[var_counts_decreasingh]), x=var_counts_decreasingh, type="s", 
     col=rgb(0,0,1,1), xlim=c(0,85), ylim=c(0,nrow(dc)),
     xlab="Number of variables (non-missing)", ylab="Number of observations")
lines(y=cumsum(play[var_counts_decreasingc]), x=var_counts_decreasingc, type="s", 
      col=rgb(1,0,0,1), xlab="", ylab="")
legend(x="topright",legend=c("Hospital","Cohort"), col=c(rgb(0,0,1,1),rgb(1,0,0,1)), 
       lty=1, bty="n")
dev.off()


###############################################################################
########################## Scratch code #######################################
###############################################################################


## compare the liver enlargment variables (expected to be correlated) ##
svg(paste(resultsDir,"higado_explore.svg", sep="")) #pngs are more protable but don't scale (else use svg)
par(mfrow=c(2,1)) 
hist(df[which(df$Higado<3),"Hepatomegalia_mm"], 
     xlim=c(0,200), main="For Higado<3", xlab="Hepatomegalia_mm")
hist(df[which(df$Higado>=3),"Hepatomegalia_mm"], 
     xlim=c(0,200), main="For Higado>=3", xlab="Hepatomegalia_mm")
par(1,1)
dev.off()


## FYI these variables consist of most general symptoms and hemmorrhaging variables, plus 7 blood counts
varinfo_88noMiss = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% covarDEN_88noMiss==T),]


if(T==F){
  vars_continuous = c(as.character(varinfo[which(varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
  vars_binary = c(as.character(varinfo[which(varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
  
  vars_continuous = c("Higado","Hepatomegalia_mm","freq_card","Higado")
  par(mfrow=c(2,2)) 
  for(var in vars_continuous){
    hist(df[,var], xlab=var, col=3)
  }
  par(1,1)
  
  #This will work (long way)
  h1 = hist(df[,"age"])
  h2 = hist(df[,"DaysSick"])
  h3 = hist(df[,"Presion_Arterial_Dias"])
  h4 = his(df[,"freq_card"])
  print(h1, position=c(0, .6, 1, 1), more=TRUE)
  print(h2, position=c(0, 0, 1, .4))
  
  par(mfrow=c())
  
  #correlation heatmap
  #http://strengejacke.wordpress.com/2013/04/18/examples-for-sjplotting-functions-including-correlations-and-proportional-tables-with-ggplot-rstats/
    #alternative: http://www.cookbook-r.com/Graphs/Correlation_matrix/ --- this one is better for black and white publications
  vars_continuous = c("Higado","Hepatomegalia_mm","freq_card","Higado","freq_resp")
  sjp.corr(df[,vars_continuous], hideLegend=F,, decimals = 2) #this works
  vars_binary = c("Nariz","Hematoma","Hematuria","Hematemesis","Hemoptisis","Melena" )
  sjp.corr(df[,vars_binary], hideLegend=F,showCorrelationPValues = F, decimals = 2) #need to convert binary vars to numeric
  
  ## performance barchart
  #http://www.thecoatlessprofessor.com/programming/creating-stacked-barplot-and-grouped-barplot-in-r-using-base-graphics-no-ggplot2
  par(mar = c(5.1, 4.1, 4.1, 7.1), xpd = TRUE)
  barplot(df, col = heat.colors(length(rownames(df))), width = 2, beside = TRUE)
  legend("topright", inset = c(-0.25, 0), fill = heat.colors(length(rownames(df))), 
         legend = rownames(df))
  
  ## transparent overlapping graphs
  #http://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  multiplot(h1, h2, h3, h4, cols=2) #did not work
}

multhist(list(dh$DaysSick, dc$DaysSick))
sjp.grpfrq(as.factor(d2$DaysSick), d2$Study)
hist(t(cbind(dh$DaysSick, dc$DaysSick)), col=c("red","lightblue"), plot=T)

plotH(x=collapsed$plot_DaysSick, y=collapsed$initialSevere.mean,
      ylab="Proportion displaying severe symptoms upon arrival", xlab="Days since fever onset",
      col=rgb(0,0,1,.5), axes=NULL)


