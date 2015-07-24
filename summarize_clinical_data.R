
### Basic summary of clinical and lab variables for qual ###

####### Find newer, v2 of this code !!!!!!!!!!!!!!!!!! ############









## Summarize lab and clinical data ##
summarize_clinical(clinical_full_clean)

###############################################################################
###### Data description section of Qual presentation -- sample sizes ##########
###############################################################################

print("DxFinal4cat")
table(clinical_D1_clean$DxFinal4cat)
table(clinical_full_clean[-which(clinical_full_clean$serum == 1),"DxFinal4cat"])

print("DaysSick")
table(clinical_D1_clean$DaysSick)
table(clinical_full_clean[-which(clinical_full_clean$serum == 1),"DaysSick"])

print("ClasificacionPrimerDia by DxFinal4cat")
print(table(clinical_D1_clean$ClasificacionPrimerDia, clinical_D1_clean$DxFinal4cat))
table(clinical_full_clean[-which(clinical_full_clean$serum == 1),]$ClasificacionPrimerDia, clinical_full_clean[-which(clinical_full_clean$serum == 1),]$DxFinal4cat)

print("Initial Dx of severe dengue patients, by DaysSick")
table(clinical_full_clean[which(clinical_full_clean$DxFinal4cat == "DHF" | clinical_full_clean$DxFinal4cat == "DSS"),]$DaysSick, 
      clinical_full_clean[which(clinical_full_clean$DxFinal4cat == "DHF" | clinical_full_clean$DxFinal4cat == "DSS"),]$ClasificacionPrimerDia)

print("Number of patients with no missing clinical information (restricting to vars collected at both hospital and clinic settings)")
num_miss_by_patient = apply(clinical_full_clean[,c(covarlist_CohortRestrict)], MARGIN=1, FUN=function(x) sum(is.na(x)))
table(num_miss_by_patient)


###############################################################################
############### Summary of clinical variables for qual ########################
###############################################################################

#create binary indicators for the general variables that were not binary
XD_binary = create_binary_variables(clinical_full_clean) #function created in "clean_data_functions.R"

## get list of clinical variables to include (returns list of variable.name.in.final.data)
varlist = get_clinic_var_list(clinic_varsD, outcome="either", 
                                        eliminate_vars_with_missings=F, eliminate_constant_vars=T, 
                                        eliminate_vars_with_minXnomiss=50,
                                        XD=clinical_full_clean, restrict_to_cohort_vars=F, restrict_to_hospit_vars=T)

## overview of these variables
contents(clinical_full_clean[,varlist])
## subset the data with variable info to include only rows that we care about
varinfo = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% varlist==T),]
## view variable absence in cohort data
table(varinfo[which(varinfo$in.cohort.data==0),]$Variable.Category, varinfo[which(varinfo$in.cohort.data==0),]$CC_type)
sum(varinfo$in.cohort.data==0, na.rm=T) #38 (so 47 present)
## view variable categories and types (binary, continuous, etc.)
    #in displaying output, put "Laboratory-General", "Laboratory-Clinical lab" and "Laboratory-Virological" with "Laboratory-Clinical lab"
    #note that I consider "daysSick" to be "General - symptom" while it had been in "Demographics/Gen Info" group
table(varinfo$Variable.Category, varinfo$CC_type)
table(varinfo$CC_category, varinfo$CC_type) #slightly more detailed than what I present in qual, but same basic groupings
#of the variables in varlist, list those which cannot be used to predict DHF/DSS  
varinfo[which(varinfo$Use.in.DF.vs.DHF.DSS.prediction!=1),c("variable.name.in.final.data","Variable.Category")]
#of the variables in varlist, list those which cannot be used to predict ND vs DENV  
varinfo[which(varinfo$Use.in.ND.vs.DEN.prediction!=1),c("variable.name.in.final.data","Variable.Category")]

## missingness
var_missCount = sapply(clinical_full_clean[varlist], function(x) sum(is.na(x)))
table(var_missCount)
var_nomissCount = sapply(clinical_full_clean[varlist], function(x) sum(!is.na(x)))
table(var_nomissCount)
#variables with no missings (33)
var_noMiss = c(as.character(varlist[var_missCount==0])) #
varinfo[which(varinfo$variable.name.in.final.data %in% var_noMiss==T), c("variable.name.in.final.data","Variable.Category")]
#variables with some missings (52)
var_someMiss = c(as.character(varlist[var_missCount>0]))
var_someMiss_info = varinfo[which(varinfo$variable.name.in.final.data %in% var_someMiss==T), c("variable.name.in.final.data","Variable.Category")]
#For now I include jandice as a "general" variable, though it could be recategorized as blood lab (currently unsure how diagnosis is made)
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
table(clinical_full_clean[which(clinical_full_clean$DxFinal4cat!="ND"),"IR"], useNA=c("always")) #missing 37/964
table(clinical_full_clean[which(clinical_full_clean$DxFinal4cat!="ND"),"PCR"], useNA=c("always")) #missing 131/964
#graph of number of missing values per variable, organized by variable category
png(paste(resultsDir,"Clinical_vars_missingness.png", sep=""), width=960, height=600)
par(mfrow=c(1,3), las=2, cex.lab=2) # make label text perpendicular to axis
par(mar=c(5,10,4,2)) # increase y-axis margin. c(bottom, left, top, right) 
barplot(var_missCount[var_someMiss_gen], horiz=T, cex.names=1.2, beside=T, col=4, main="General")
barplot(var_missCount[var_someMiss_lab], horiz=T, cex.names=1.2, beside=T, col=4, main="Blood lab")
barplot(var_missCount[var_someMiss_ultraX], horiz=T, cex.names=1.2, beside=T, col=4, main="Ultrasound & X-ray")
dev.off()

#Variables available per patient
patient_nomissCount = apply(clinical_full_clean[varlist], MARGIN=1, function(x) sum(!is.na(x)))
play = table(patient_nomissCount)
var_counts_decreasing = sort(names(play), decreasing=TRUE)
cumsum(play[var_counts_decreasing])
#graph
png(paste(resultsDir,"Clinical_vars_perPatient.png", sep=""), width=800, height=500)
par(mar=c(5,6,4,2)+.1) #default is  c(5, 4, 4, 2) + 0.1
plot(y=cumsum(play[var_counts_decreasing]), x=var_counts_decreasing, type="s", cex.lab=2, cex.axis=2, col=4,
     xlab="Number of non-missing variables", ylab="Number of observations")
dev.off()

## compare the liver enlargment variables (expected to be correlated) ##
svg(paste(resultsDir,"higado_explore.svg", sep="")) #pngs are more protable but don't scale (else use svg)
par(mfrow=c(2,1)) 
hist(clinical_full_clean[which(clinical_full_clean$Higado<3),"Hepatomegalia_mm"], 
     xlim=c(0,200), main="For Higado<3", xlab="Hepatomegalia_mm")
hist(clinical_full_clean[which(clinical_full_clean$Higado>=3),"Hepatomegalia_mm"], 
     xlim=c(0,200), main="For Higado>=3", xlab="Hepatomegalia_mm")
par(1,1)
dev.off()

## Demographics/ general (todo: increase font size of axis labeling)
png(paste(resultsDir,"Demographics.png", sep=""), width=800, height=500) #pngs are more portable but don't scale (else use svg)
par(mfrow=c(1,2)) 
plot(clinical_full_clean[,"Sexo"], main="Gender", col=4, ylab="Frequency", cex.lab=1.5, cex.axis=1.5, cex.main=2)
hist(clinical_full_clean[,"age"], main="Age", xlab="", col=4, cex.lab=1.5, cex.axis=1.5, cex.main=2,)
par(1,1)
dev.off()

## General symptoms, binary and categorical (todo: add categorical)
var_names = c(as.character(varinfo[which(varinfo$Variable.Category=="Gen Sign/Symptom" &
                                              varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as per)
get_mean = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals = sapply(clinical_full_clean[var_names], FUN=get_mean)
#todo: sort in order of mean value 
#todo: create better names to display (cc_name)
#graph
png(paste(resultsDir,"General_symtoms_binary.png", sep=""), width=800, height=500)
par(mfrow=c(1,1),las=2) # make label text perpendicular to axis
par(mar=c(5,10,4,2)) # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(mean_vals, horiz=T, names.arg=c(var_names), xlim=c(0,1), xlab="Proportion affirmative", 
        cex.names=0.8, beside=T, xaxt="n", col=4)
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=mean_vals+0.05, y=bplt, labels = sprintf("%.3f",mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

## General symptoms, continuous
var_names = c(as.character(varinfo[which(varinfo$Variable.Category=="Gen Sign/Symptom" &
                                        varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"General_symtoms_continuous.png", sep=""), width=800, height=500)
par(mfrow=c(2,3), mai=c(1,.3,.1,.2), cex.lab=2, cex.axis=1.5) 
for(var in var_names){
  hist(clinical_full_clean[,var], xlab=var, col=4, main="", las=0, ylab="", freq=F, )
}
dev.off()


## General indicators of hemorrhaging, binary and categorical
var_names1 = c(as.character(varinfo[which(varinfo$Variable.Category=="Hem Sign/Symptom" &
                                           varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as percent)
get_mean1 = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals1 = sapply(XD_binary[var_names1], FUN=get_mean1)
#the only gen-hem variable that is not already binary is the torniquete one so I created binary form of it
var_names2 = c("is.torniquete10plus","is.torniquete20plus")
get_mean2 = function(x) mean(as.numeric(x), na.rm=T)
mean_vals2 = sapply(XD_binary[var_names2], FUN=get_mean2)
mean_vals = c(mean_vals1, mean_vals2)
var_names = c(var_names1, var_names2)
#todo: sort in order of mean value 
#todo: create better names to display (cc_name)
#graph
plot.new()
png(paste(resultsDir,"General_hemorrhage_binary.png", sep=""), width=800, height=500)
par(mfrow=c(1,1),las=2) # make label text perpendicular to axis
par(mar=c(5,10,4,2)) # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(mean_vals, horiz=T, names.arg=c(var_names), xlim=c(0,1), xlab="Proportion affirmative", 
        cex.names=1, beside=T, xaxt="n", col=4)
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=mean_vals+0.05, y=bplt, labels = sprintf("%.3f",mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

## Blood lab count variables, continuous
plot.new()
var_names = c(as.character(varinfo[which(varinfo$CC_category=="Blood - count" &
                                            varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"Blood_counts_continuous.png", sep=""), width=880, height=550)
par(mfrow=c(2,4))
for(var in var_names){
  hist(clinical_full_clean[,var], xlab=var, cex.lab=2, cex.axis=1.5, col=4, main="", las=0, ylab="", freq=F, )
}
dev.off()
## Blood lab chemistry variables, continuous
plot.new()
var_names = c(as.character(varinfo[which(varinfo$CC_category=="Blood - chemistry" &
                                           varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"Blood_chemistry_continuous_part1.png", sep=""), width=800, height=500)
par(mfrow=c(2,4))
for(var in var_names[1:8]){
  hist(clinical_full_clean[,var], xlab=var, cex.lab=2, cex.axis=1.5, col=4, main="", las=0, ylab="", freq=F, )
}
dev.off()
png(paste(resultsDir,"Blood_chemistry_continuous_part2.png", sep=""), width=880, height=550)
par(mfrow=c(3,3))
for(var in var_names[9:17]){
  hist(clinical_full_clean[,var], xlab=var, cex.lab=2, cex.axis=1.5, col=4, main="", las=0, ylab="", freq=F, )
}
dev.off()

## Ultrasound, continuous
plot.new()
var_names = c(as.character(varinfo[which((varinfo$CC_category=="Ultrasound"|varinfo$CC_category=="X-ray") &
                                           varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
png(paste(resultsDir,"ultrasound_continuous.png", sep=""), width=800, height=300)
par(mfrow=c(1,3))
for(var in var_names){
  hist(clinical_full_clean[,var], xlab=var, cex.lab=2, cex.axis=1.5, col=4, main="", las=0, ylab="", freq=F, )
}
dev.off()
### Ultrasound and xray, binary
plot.new()
var_names = c(as.character(varinfo[which((varinfo$CC_category=="Ultrasound"|varinfo$CC_category=="X-ray") &
                                           varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
#get mean values (coded 0/1 so this is same as percent)
get_mean = function(x) mean(as.numeric(x)-1, na.rm=T)
mean_vals = sapply(clinical_full_clean[var_names], FUN=get_mean)
png(paste(resultsDir,"ultrasoundXray_binary.png", sep=""), width=800, height=300)
par(mfrow=c(1,1), las=2) # make label text perpendicular to axis
par(mar=c(5,12,4,2)) # increase y-axis margin. c(bottom, left, top, right) 
bplt = barplot(mean_vals, horiz=T, names.arg=c(var_names), xlim=c(0,1), xlab="Proportion affirmative", 
        cex.names=1, beside=T, xaxt="n", col=4)
axis(side=1, at=seq(0,1,.2), las=0) #this is how I got x-axis labels to be parallel
text(x=mean_vals+0.06, y=bplt, labels = sprintf("%.3f",mean_vals), xpd=TRUE) #add value labels to bars
dev.off()

###############################################################################
########################## Scratch code #######################################
###############################################################################

## FYI these variables consist of most general symptoms and hemmorrhaging variables, plus 7 blood counts
varinfo_88noMiss = clinic_varsD[which(clinic_varsD$variable.name.in.final.data %in% covarDEN_88noMiss==T),]


if(T==F){
  vars_continuous = c(as.character(varinfo[which(varinfo$CC_type=="continuous"),"variable.name.in.final.data"]))
  vars_binary = c(as.character(varinfo[which(varinfo$CC_type=="binary"),"variable.name.in.final.data"]))
  
  vars_continuous = c("Higado","Hepatomegalia_mm","freq_card","Higado")
  par(mfrow=c(2,2)) 
  for(var in vars_continuous){
    hist(clinical_full_clean[,var], xlab=var, col=3)
  }
  par(1,1)
  
  #This will work (long way)
  h1 = histogram(clinical_full_clean[,"age"])
  h2 = histogram(clinical_full_clean[,"DaysSick"])
  h3 = histogram(clinical_full_clean[,"Presion_Arterial_Dias"])
  h4 = histogram(clinical_full_clean[,"freq_card"])
  print(h1, position=c(0, .6, 1, 1), more=TRUE)
  print(h2, position=c(0, 0, 1, .4))
  
  par(mfrow=c())
  
  #correlation heatmap
  #http://strengejacke.wordpress.com/2013/04/18/examples-for-sjplotting-functions-including-correlations-and-proportional-tables-with-ggplot-rstats/
    #alternative: http://www.cookbook-r.com/Graphs/Correlation_matrix/ --- this one is better for black and white publications
  source("sjPlotCorr.R")
  vars_continuous = c("Higado","Hepatomegalia_mm","freq_card","Higado","freq_resp")
  sjp.corr(clinical_full_clean[,vars_continuous], hideLegend=F,showCorrelationPValues = F, decimals = 2) #this works
  vars_binary = c("Nariz","Hematoma","Hematuria","Hematemesis","Hemoptisis","Melena" )
  sjp.corr(clinical_full_clean[,vars_binary], hideLegend=F,showCorrelationPValues = F, decimals = 2) #need to convert binary vars to numeric
  
  ## performance barchart
  #http://www.thecoatlessprofessor.com/programming/creating-stacked-barplot-and-grouped-barplot-in-r-using-base-graphics-no-ggplot2
  par(mar = c(5.1, 4.1, 4.1, 7.1), xpd = TRUE)
  barplot(data, col = heat.colors(length(rownames(data))), width = 2, beside = TRUE)
  legend("topright", inset = c(-0.25, 0), fill = heat.colors(length(rownames(data))), 
         legend = rownames(data))
  
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


