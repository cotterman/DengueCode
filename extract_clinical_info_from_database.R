
#############################################################################
################### Add hospital data from scratch  #########################
#############################################################################

## Admission and General Information ##

#This data contains just one obs per patient
Datos_Admision = read.delim(paste(clinical_inputsDir,"Datos_Admision.txt", sep=""), header=TRUE, nrows=35000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Datos_Admision["code"] = apply(X=Datos_Admision, MARGIN = 1, FUN=fcode, var="Codigo")
#clean variables
Datos_Admision$Fecha_Estudio = as.Date(Datos_Admision$Fecha_Estudio, format = "%m/%d/%Y")
Datos_Admision$Hora_Estudio = substr(as.character(Datos_Admision$Hora_Estudio),12, 20)
Datos_Admision$Time_Estudio = chron(times = Datos_Admision$Hora_Estudio)
Datos_Admision$Fecha_Nac = as.Date(Datos_Admision$Fecha_Nac, format = "%m/%d/%Y")
Datos_Admision$age = as.numeric(Datos_Admision$Fecha_Estudio - Datos_Admision$Fecha_Nac)/365.25
#keep variables contained in the clinic_varsH list, plus identifiers
get_mismatches(clinic_varsH, colnames(Datos_Admision), "prelim", "added")
matches = intersect(clinic_varsH, colnames(Datos_Admision)) #only 6 vars in common ("Dolor_Abdominal", "Tos", "Cefalea", "Rash", "Dolor_ocular", "Vomito")
#there are many more variables that are part of clinic_varsH but with a prefix (e.g., SC.Hipotension.","EFC.Hemoconcentracion)
#where did these variables come from?  (I don't see corresponding questions on the "Admission and General Information" form)
#merge with other cleaned hospital data and compare values -- more obs in Datos_Admision
t = merge(clinical_hospitD_clean, Datos_Admision, by="code", all=T)
#View(subset(t, is.na(FTM))) #115 extra obs in data, all from admission year 2014
table(t$Dolor_Abdominal.x, t$Dolor_Abdominal.y) #not a super strong correspondence (strange that some admission "yes" responses are "no"s in 24hr data...)
#compare date of admission and birthdate values
t = merge(clinical_hospitD, Datos_Admision, by="code", all=T)
test = t$Fecha_Nac.y - t$Fecha_Nac.x
table(test) #perfect match
test = t$FTM - t$Fecha_Estudio
table(test) #different for 4 out of 1650 obs; identical for others


## Clinical History ##  

# This data contains just one observation per person and matches perfectly with Datos_Admision data
Historia_Clinica = read.delim(paste(clinical_inputsDir,"Historia_Clinica.txt", sep=""), header=TRUE, nrows=35000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Historia_Clinica["code"] = apply(X=Historia_Clinica, MARGIN = 1, FUN=fcode, var="Codigo")
get_mismatches(clinic_varsH, colnames(Historia_Clinica), "prelim", "added")
intersect(clinic_varsH, colnames(Historia_Clinica)) #only 9 vars in common ("Dolor_Abdominal","Tos", "Perdida_Apetito","Rash", "Encias", "Nariz", "Hematemesis", "Melena", "Vaginal")
#merge with admissions data -- a perfect match (n=1765)
t = merge(Historia_Clinica, Datos_Admision, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, Historia_Clinica, by="code", all=T)
table(t$Dolor_Abdominal.x, t$Dolor_Abdominal.y) #strong correlation but more "yes" responses in "24h" data

## Laboratory and Office ##

Lab = read.delim(paste(clinical_inputsDir,"Lab_Clinico_Pacientes_Hospitalizados.txt", sep=""), header=TRUE, nrows=35000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Lab["code"] = apply(X=Lab, MARGIN = 1, FUN=fcode, var="Codigo")
#sort by patient and date/time
Lab$FechaH = as.Date(Lab$FechaH, format = "%m/%d/%Y")
Lab[which(Lab$Turno=="6AM a 6PM"),"shift"] = 1
Lab[which(Lab$Turno=="6PM a 6AM"),"shift"] = 2
t = Lab[order(Lab$code, Lab$FechaH, Lab$shift),]
#verify only one observation per unique combo of codigo, fechaH, Turno
sum(duplicated(t[,c("code", "FechaH", "shift")])) #sum of zero means all unique.  good.
#replace 9999 values with missings
t[,c("Albumina", "Plaquetas","Proteina_tot","Hematocrito")][t[,c("Albumina", "Plaquetas","Proteina_tot","Hematocrito")]==9999] = NA
#keep summary of variables that feed into final Dx algorithm for each patient
t_min = aggregate(t[,c("Albumina", "Plaquetas","Proteina_tot")], by=list(t$code), FUN=function(x) min(x, na.rm=TRUE))
#Aggretate gives inf when observation has all missing values.  Want to replace these with NA
t_min[,c("Albumina", "Plaquetas","Proteina_tot")][t_min[,c("Albumina", "Plaquetas","Proteina_tot")]==Inf] = NA
#take max level of Hematocrito (this var takes values btwn 15 and 58)
t_max = aggregate(t[,c("Hematocrito")], by=list(t$code), FUN=function(x) max(x, na.rm=TRUE))
names(t_max)[names(t_max)=="x"] <- "Hematocrito" #rename
t_max[,c("Hematocrito")][t_max[,c("Hematocrito")]==-Inf] = NA
#combine the max and min data
Lab_SummObs = merge(t_min, t_max, by="Group.1", all=T)
#keep just the first observation for each patient.  This method isnt great for large data, but works here
t = by(t, t$code, head, n=1) #creates list
Lab_1stOb = do.call("rbind", as.list(t)) #turns list into data frame
get_mismatches(clinic_varsH, colnames(Lab_1stOb), "prelim", "added")
intersect(clinic_varsH, colnames(Lab_1stOb)) # 10 in common ("Albumina","RelA.G","Colesterol","Globulina","HDL","Hemoglobina","LDL.C","Plaquetas","Proteina_tot","Leucocitos")
#merge with admissions data -- Lab contains a subset of patients 
t = merge(Historia_Clinica[,c("code","Perdida_Apetito")], Lab_1stOb, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, Lab_1stOb, by="code", all=T)
cor(t$Albumina.y, t$Albumina.x, use="complete.obs") #correlation of .58)

## Clinical Information ##

PH = read.delim(paste(clinical_inputsDir,"Pacientes_Hospitalizados.txt", sep=""), header=TRUE, nrows=35000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
PH["code"] = apply(X=PH, MARGIN = 1, FUN=fcode, var="Codigo")
#sort by patient and date/time
PH$FechaH = as.Date(PH$FechaH, format = "%m/%d/%Y")
PH[which(PH$Turno=="6AM a 6PM"),"shift"] = 1
PH[which(PH$Turno=="6PM a 6AM"),"shift"] = 2
t = PH[order(PH$code, PH$FechaH, PH$shift),]
#verify only one observation per unique combo of codigo, fechaH, Turno
sum(duplicated(t[,c("code", "FechaH", "shift")])) #sum of zero means all unique.  good.
get_mismatches(clinic_varsH, colnames(t), "prelim", "added")
intersect(clinic_varsH, colnames(t))  #28 in common, including repeats from above, like tos, vomito, rash
DHF_vars = intersect(DxDHF_varsH, colnames(t)) #10 binary variables currently coded with values "N" and "S"
#Purpura Hematoma Venopuncion Hipermenorrea Vaginal Subconjuntival Hemoptisis Melena Hematemesis Ascitis
additional_binary_vars =c("Epistaxis", "Encia", "Edema", "periorbitaria","Facial","Miembros_inferiores","Hidrocele","Generalizado","Fovea",
                          "Derrame_Pleural","Derrame_.pleural_izq_Rx","Derrame_.pleural_der_Rx","Derrame_Clinico","Derecho","Izquierdo",
                          "Esquimosis","Petequia")
intersect(DxDSS_varsH, colnames(t)) #none
#aggregate data by patient
apply(t[,c(DHF_vars, additional_binary_vars)], 2, FUN=function(x) table(x))
for (name in c(DHF_vars, additional_binary_vars)) {
  #print(name)
  t$temp = NULL
  t[which(t[name]=="S"),"temp"]=1
  t[which(t[name]!="S"),"temp"]=0
  class(t$temp)
  #drop old version and rename new version
  t[name] = NULL
  t[name] = t$temp
  t$temp = NULL
}
apply(t[,c(DHF_vars, additional_binary_vars)], 2, FUN=function(x) table(x))
t$is.torniquete20plus = ifelse(t$Prueba_torniquete==20,1,0)
PH_SummObs = aggregate(t[,c(DHF_vars, additional_binary_vars, "is.torniquete20plus")], by=list(t$code), FUN=function(x) max(x, na.rm=TRUE))
#keep just the first observation for each patient.  This method isnt great for large data, but works here
t = by(t, t$code, head, n=1) #creates list
PH_1stOb = do.call("rbind", as.list(t)) #turns list into data frame
#merge with lab data -- contains exact same 1389 patients as lab 
t = merge(Lab_1stOb[,c("code","CPK")], PH_1stOb, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, PH_1stOb, by="code", all=T)
table(t$Tos.x, t$Tos.y) #there are yes values in new clinical info that are no in 24hr data 
table(t$Dia_Defervesencia) #this data appears to be from outside of the 24 hr window (be careful with taking 1st obs)


## Vital Signs Report ##

#This data has multiple observations per patient, day, and shift (and even times)
PH2 = read.delim(paste(clinical_inputsDir,"Paciente_Hospitalizado_2.txt", sep=""), header=TRUE, nrows=35000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
PH2["code"] = apply(X=PH2, MARGIN = 1, FUN=fcode, var="Codigo")
#sort by patient and date/time
PH2$FechaH = as.Date(PH2$FechaH, format = "%m/%d/%Y")
PH2[which(PH2$Turno=="6AM a 6PM"),"shift"] = 1
PH2[which(PH2$Turno=="6PM a 6AM"),"shift"] = 2
PH2$Hora_PH2 = substr(as.character(PH2$Hora),12, 20)
PH2$Hora_PH2 = chron(times = PH2$Hora_PH2)
t = PH2[order(PH2$code, PH2$FechaH, PH2$shift, PH2$Hora_PH2),]
ranks = with(t, ave(code, code, FUN=function(x) rank(x, ties.method="first")))
intersect(DxDHF_varsH, colnames(t)) #none
DSS_vars = intersect(DxDSS_varsH, colnames(t)) #"Presion_Arterial_Sist""Presion_Arterial_Dias""Palidez""Sudoracion""Escalofrio""Llenado_Capilar"
apply(t[,c(DSS_vars)], 2, function(x) table(x, useNA="always")) #take a look at variable values
#clean relevant variables 
#replace 9999 values with missings
t[,c("Presion_Arterial_Sist","Presion_Arterial_Dias","Temperatura")][t[,c("Presion_Arterial_Sist","Presion_Arterial_Dias","Temperatura")]==9999] = NA
for (name in c("Palidez","Sudoracion","Escalofrio")) {
  #print(name)
  t$temp = NULL
  t[which(t[name]=="S"),"temp"]=1
  t[which(t[name]!="S"),"temp"]=0
  class(t$temp)
  #drop old version and rename new version
  t[name] = NULL
  t[name] = t$temp
  t$temp = NULL
}
t$ExtremidadesFrias = ifelse(t$Extremidades=="Frialdad distal" |  
                               t$Extremidades=="Frias y pegajosas", 1, 0) #educated guess
t$is.PoorCapillaryRefill = NULL
t[which(t["Llenado_Capilar"]=="> 2 seg" | t["Llenado_Capilar"]==">3\"" | t["Llenado_Capilar"]=="3\""),"is.PoorCapillaryRefill"]=1
t[which(t["Llenado_Capilar"]=="< 2 seg" | t["Llenado_Capilar"]=="2" | t["Llenado_Capilar"]=="2\""),"is.PoorCapillaryRefill"]=0
t$is.NarrowPulse = ifelse(t$Presion_Arterial_Sist - t$Presion_Arterial_Dias <= 20, 1, 0)
t$is.pulse_danger = ifelse(t$Pulso=="R" | t$Pulso=="N", 1, 0) #danger if rapid or not palpable
#add age since it is used in definition of hypotension
t = merge(Datos_Admision[,c("code","age")], t, by="code", all.y=T) #keep only obs in PH2 data
t$is.hypotension = ifelse( (t$Presion_Arterial_Sist<80 & t$age<5) | (t$Presion_Arterial_Sist<90 & t$age>=5), 1, 0)
PH2_SummObs = aggregate(t[,c("Palidez","Sudoracion","Escalofrio","is.PoorCapillaryRefill","ExtremidadesFrias",
                             "is.hypotension","is.NarrowPulse","is.pulse_danger","Temperatura")], 
                        by=list(t$code), FUN=function(x) max(x, na.rm=TRUE))
PH2_1stOb = t[ranks==1,] #INCORRECT WAY TO CHOSE EARLIEST OBS -- REDO IF USING THIS DATA!!
#note that none of these variables (aside from temperature) are used in definition of DHF, though more are used for DSS
intersect(clinic_varsH, colnames(PH2_1stOb)) # "Escalofrio","Presion_Arterial_Dias","freq_card","Palidez","Llenado_Capilar","Pulso","freq_resp","Sudoracion","Presion_Arterial_Sist","Temperatura"  
#merge with lab data -- contains exact same 1389 patients as lab 
t = merge(Lab_1stOb[,c("code","CPK")], PH2_1stOb, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, PH2_1stOb, by="code", all=T)
cor(t$Temperatura.x, t$Temperatura.y, use="complete.obs") 



## External Consultation Form ##

#This data has one observation per patient/date
Consulta_Externa = read.delim(paste(clinical_inputsDir,"Consulta_Externa.txt", sep=""), header=TRUE, nrows=5000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Consulta_Externa["code"] = apply(X=Consulta_Externa, MARGIN = 1, FUN=fcode, var="Codigo")
#sort by patient and date
Consulta_Externa$fechacons = as.Date(Consulta_Externa$fechacons, format = "%m/%d/%Y")
t = Consulta_Externa[order(Consulta_Externa$code, Consulta_Externa$fechacons),]
#verify only one observation per unique combo of codigo, fechacons
#note: when fechacons (from Consulta_Externa) = Fecha_Estudio (from Datos_Admision), then dia=1 (from Constulta_Externa)
sum(duplicated(t[,c("code", "fechacons")])) #sum of zero.  good.
ranks = with(t, ave(code, code, FUN=seq)) #works so long as already sorted
## prepare to aggregate using variables relevant for DHF and DSS diagnoses
matches = intersect(clinic_varsH, colnames(t)) #35 variables in common, including Artralgia, Ascitis, Tos, Vaginal, Melena, etc.
DHF_vars = intersect(DxDHF_varsH, colnames(t)) #10 binary variables currently coded with values "N","S", and "D"
DSS_vars = intersect(DxDSS_varsH, colnames(t)) #"Presion_Arterial_Sist" "Presion_Arterial_Dias" "Llenado_Capilar"
apply(t[,c(DHF_vars)], 2, function(x) table(x, useNA="always")) #take a look at variable values
apply(t[,c(DSS_vars)], 2, function(x) table(x, useNA="always")) #take a look at variable values
additional_binary_vars =c("Epistaxis", "Encia", "Edema", "periorbitaria","Facial","Miembros_inferiores","Hidrocele","Generalizado","Fovea",
                          "Derrame_.pleural_izq_Rx","Derrame_.pleural_der_Rx","Derrame_Clinico","Esquimosis","Petequia")
apply(t[,c(DHF_vars, additional_binary_vars)], 2, FUN=function(x) table(x))
#clean relevant variables 
#replace 9999 values with missings
t[,c("Presion_Arterial_Sist","Presion_Arterial_Dias")][t[,c("Presion_Arterial_Sist","Presion_Arterial_Dias")]==9999] = NA
for (name in c(DHF_vars, additional_binary_vars)) {
  #print(name)
  t$temp = NULL
  t[which(t[name]=="S"),"temp"]=1
  t[which(t[name]=="N"),"temp"]=0
  class(t$temp)
  #drop old version and rename new version
  t[name] = NULL
  t[name] = t$temp
  t$temp = NULL
}
t$is.torniquete20plus = ifelse(t$Prueba_torniquete==20,1,0)
table(t$Llenado_Capilar)
t$is.PoorCapillaryRefill = NULL
t[which(t["Llenado_Capilar"]=="> 2 seg" | t["Llenado_Capilar"]==">3\"" | t["Llenado_Capilar"]=="3\""),"is.PoorCapillaryRefill"]=1
t[which(t["Llenado_Capilar"]=="< 2 seg" | t["Llenado_Capilar"]=="2"),"is.PoorCapillaryRefill"]=0
table(t$Llenado_Capilar, t$is.PoorCapillaryRefill, useNA="always")
t$is.NarrowPulse = ifelse(t$Presion_Arterial_Sist - t$Presion_Arterial_Dias <= 20, 1, 0)
#add age since it is used in definition of hypotension
t = merge(Datos_Admision[,c("code","age")], t, by="code", all=T)
t$is.hypotension = ifelse( (t$Presion_Arterial_Sist<80 & t$age<5) | (t$Presion_Arterial_Sist<90 & t$age>=5), 1, 0)
#create data with aggregation of vars across time
Consulta_Externa_SummObs = aggregate(t[,c(DHF_vars,additional_binary_vars,"is.torniquete20plus","is.PoorCapillaryRefill","is.NarrowPulse","is.hypotension")], 
                                     by=list(t$code), FUN=function(x) max(x, na.rm=TRUE))
#create data containing just first obs per patient
Consulta_Externa_1stOb = t[ranks==1,] 
get_mismatches(clinic_varsH, colnames(Consulta_Externa_1stOb), "prelim", "added")
#merge with admissions data -- perfectly matches all 1765 patients 
t = merge(Datos_Admision[,c("code","Sexo")], Consulta_Externa_1stOb, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, Consulta_Externa_1stOb, by="code", all=T)
table(t$Ascitis.x, t$Ascitis.y, useNA="no") #good correlation with simply more cases reported in the 24hr data


## Additional lab info and possibly redundant info with other sheets ##

#This contains a ton of info that other CRF forms are already accounting for...what is this all about?
Lab_Consulta_Externa = read.delim(paste(clinical_inputsDir,"Lab_Clinico_Consulta_Externa.txt", sep=""), header=TRUE, nrows=5000) 
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Lab_Consulta_Externa["code"] = apply(X=Lab_Consulta_Externa, MARGIN = 1, FUN=fcode, var="Codigo")
#sort by patient and date
t = Lab_Consulta_Externa[order(Lab_Consulta_Externa$code, Lab_Consulta_Externa$Dia),]
#verify only one observation per unique combo of codigo, dia
sum(duplicated(t[,c("code", "Dia")])) #sum of zero.  good.
ranks = with(t, ave(code, code, FUN=seq)) #works so long as already sorted
Lab_Consulta_Externa_1stOb = t[ranks==1,] 
intersect(clinic_varsH, colnames(Lab_Consulta_Externa_1stOb)) #"Albumina","RelA.G","Colesterol","Globulina","HDL","Hemoglobina","LDL.C","Plaquetas","Proteina_tot" "Leucocitos" 
DHF_vars = intersect(DxDHF_varsH, colnames(t)) #"Proteina_tot" "Albumina" "Plaquetas" (same as in Lab_Clinico_Pacientes_Hospitalizados.txt)
intersect(DxDSS_varsH, colnames(t)) #none
#replace 9999 values with missings
t[,c("Albumina", "Plaquetas","Proteina_tot","Hematocrito")][t[,c("Albumina", "Plaquetas","Proteina_tot","Hematocrito")]==9999] = NA
t_min = aggregate(t[,c("Albumina", "Plaquetas","Proteina_tot")], by=list(t$code), FUN=function(x) min(x, na.rm=TRUE))
#Aggretate gives inf when observation has all missing values.  Want to replace these with NA
t_min[,c("Albumina", "Plaquetas","Proteina_tot")][t_min[,c("Albumina", "Plaquetas","Proteina_tot")]==Inf] = NA
#take max level of Hematocrito
t_max = aggregate(t[,c("Hematocrito")], by=list(t$code), FUN=function(x) max(x, na.rm=TRUE))
names(t_max)[names(t_max)=="x"] <- "Hematocrito" #rename
t_max[,c("Hematocrito")][t_max[,c("Hematocrito")]==-Inf] = NA
#combine the max and min data
Lab_Consulta_Externa_SummObs = merge(t_min, t_max, by="Group.1", all=T)
#merge with admissions data -- perfectly matches all 1765 patients 
t = merge(Datos_Admision[,c("code","Sexo")], Lab_Consulta_Externa_1stOb, by="code", all=T)
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, Lab_Consulta_Externa_1stOb, by="code", all=T)
cor(t$Proteina_tot.x, t$Proteina_tot.y, use="complete.obs") #negative correlation!!

#Composite info -- largely redundent with "Additional outcomes_hospital study.txt" 
Resumen_Alta = read.delim(paste(clinical_inputsDir,"Resumen_Alta.txt", sep=""), header=TRUE, nrows=35000)
#reformat the sample ID variable so that it is "ID" followed by 4-digit character variable
Resumen_Alta["code"] = apply(X=Resumen_Alta, MARGIN = 1, FUN=fcode, var="Codigo")
#verify only one observation per code
sum(duplicated(Resumen_Alta[,c("code")])) #sum of zero.  good.
#contains all patients in admision data
t = merge(Datos_Admision[,c("code","Sexo")], Resumen_Alta, by="code", all=T)
#we see here that only the Manejo="H"(hospitalizado) patients got Lab data, not the Manejo="E"(externo) patients
t = merge(Lab_1stOb[,c("code","CPK")], Resumen_Alta[which(Resumen_Alta$Manejo=="E"),], by="code", all=F)
intersect(clinic_varsH, colnames(Resumen_Alta)) #no overlap
#compare with indicators in "24 hr" data
t = merge(clinical_hospitD_clean, Resumen_Alta, by="code", all=T)
table(t$WHORevisedFinal, t$DxNuevaClasificacion, useNA="no") #poor correspondance(!)
table(t$ClasificacionFinal, t$Dx_Clinico_Final, useNA="no") #poor correspondence, with ClasificacionFinal more often severe dengue

## Other data that we do not currently use ##
#Info on oral rehydration therapy (may want to incorporate into analysis in future) 
#Paciente_Hosp_Balance_Liquidos = read.delim(paste(clinical_inputsDir,"Paciente_Hosp_Balance_Liquidos.txt", sep=""), header=TRUE, nrows=35000) 
#Not sure which form this corresponds to.  Only contains a couple of variables.
#PH_3 = read.delim(paste(clinical_inputsDir,"Paciente_Hospitalizado_3.txt", sep=""), header=TRUE, nrows=35000) 


################################################################
############ Combine the patient-summarized data ###############
################################################################

# subset of patients (1389): Lab_SummObs, PH_SummObs
summObs_subset = merge(PH_SummObs, Lab_SummObs, by="Group.1", all=T)
summObs_subset$InAllNewD = 1
colnames(summObs_subset) #vars are unique

# all patients (1765): Lab_Consulta_Externa_SummObs, Consulta_Externa_SummObs, PH2_SummObs
summObs = merge(Lab_Consulta_Externa_SummObs, Consulta_Externa_SummObs, by="Group.1", all=T)
colnames(summObs) #vars are unique
summObs = merge(summObs, PH2_SummObs, by="Group.1", all=T)
colnames(summObs) #3 overlapping vars between Consulta_Externa_SummObs and PH2_SummObs (is.NarrowPulse, is.PoorCapillaryRefill, is.hypotension)
table(summObs$is.NarrowPulse.x, summObs$is.NarrowPulse.y, useNA="always") 
table(summObs$is.PoorCapillaryRefill.x, summObs$is.PoorCapillaryRefill.y, useNA="always")
table(summObs$is.hypotension.x, summObs$is.hypotension.y, useNA="always")
summObs = return_min_or_max(summObs, "is.NarrowPulse", max)
summObs = return_min_or_max(summObs, "is.PoorCapillaryRefill", max)
summObs = return_min_or_max(summObs, "is.hypotension", max)
colnames(summObs)

#combine long and short data
summObs_all = merge(summObs, summObs_subset, by="Group.1", all=T)
colnames(summObs_all)[order(colnames(summObs_all))] #nearly all variables come from 2 different sources and with diff values

#add age info and also hemoconcentration from Datos Admision (no other apparent sourse of hemoconcentration)
names(summObs_all)[names(summObs_all)=="Group.1"] <- "code"  #rename Group.1 to code
summObs_all = merge(Datos_Admision[,c("code","age","Sexo","EFC.Hemoconcentracion.")], summObs_all, by="code", all=T)

#take min values for these vars
summObs_all = return_min_or_max(summObs_all, "Plaquetas", min)
summObs_all = return_min_or_max(summObs_all, "Proteina_tot", min)
summObs_all = return_min_or_max(summObs_all, "Albumina", min)

#meet definition of thrombocytopenia if <100,000 platelets per mm3.  (Remember Plaquetas variable is /1000)
summObs_all$is.thrombocytopenia = summObs_all$Plaquetas<=100
summObs_all[which(is.na(summObs_all$Plaquetas)), "is.thrombocytopenia"] = NA
#these are defined according to appearance in definition of DHF
summObs_all$is.Hypoproteinemia = summObs_all$Proteina_tot<4.2  # <4.2 g/dl
summObs_all[which(summObs_all$age>=2),"is.Hypoproteinemia"] = summObs_all[which(summObs_all$age>=2),"Proteina_tot"]<6 #looser criteria if age>=2 
summObs_all$is.Hypoalbuminemia = summObs_all$Albumina<2  # <2 g/dl
summObs_all[which(summObs_all$age>=1),"is.Hypoalbuminemia"] = summObs_all[which(summObs_all$age>=1),"Albumina"]<3.5  #looser criteria if age>=1

#obtain list of variables for which we have a .x and .y version
vars_to_try = colnames(summObs_all)[grep("\\.x",colnames(summObs_all))]
vars_no_suffix = lapply(vars_to_try, function(var) gsub("\\.x","",var))
#for each of the qualifying variables, take the max value across the .x and .y versions
for (var in vars_no_suffix){
  print(var)
  summObs_all = test_then_return_min_or_max(summObs_all, var, max)  
  print(table(summObs_all[var]))
}
colnames(summObs_all)


###############################################################################
###################### Re-create final diagnosis variable #####################
###############################################################################

#apply definition of HematocritoElev
summObs_all$HematocritoElev = with(summObs_all, ifelse((age<=2 & Hematocrito>=40) | (age>2 & age<15 & Hematocrito>=42) | 
                                                         (age>=15 & Sexo=="F" & Hematocrito>=46) | (age<=2 & Sexo=="M" & Hematocrito>=50), 1, 0))
summObs_all$Edema_ = with(summObs_all, ifelse(Edema+periorbitaria+Facial+Miembros_inferiores+Hidrocele+Generalizado+Fovea > 1, 1, 0))
summObs_all$Derrame_ = with(summObs_all, ifelse(Derrame_Pleural+Derrame_.pleural_izq_Rx+Derrame_.pleural_der_Rx+Derrame_Clinico+Derecho+Izquierdo > 1, 1, 0))
#rename for compatibility with other data
summObs_all = rename(summObs_all, c("Epistaxis"="Nariz",
                                    "Encia"="Encias",
                                    "Petequia"="PetequiasEspontaneas",
                                    "Esquimosis"="Equimosis",
                                    "EFC.Hemoconcentracion."="Hemoconcentracion",
                                    "is.NarrowPulse" = "is.Estrechamiento",
                                    "is.PoorCapillaryRefill" = "Llenado_Capilar")) 

temp = summObs_all
#DHF
temp[is.na(temp)]<-0 #TEMPORARILY code numeric NAs as zeros.  Note: wanted to restrict to numeric vars but had trouble with apply(temp, 2, FUN=is.numeric) code
temp[is.na(temp)]<-F #TEMPORARILY code factor variable NAs as False.
temp$Hem1 = with(temp, ifelse(PetequiasEspontaneas+Equimosis+Purpura+Hematoma>0,1,0)) #check this?? Petechia, equimosis, purpura, or hematoma
temp$Hem2 = with(temp, ifelse(Venopuncion+Hipermenorrea+Vaginal+Subconjuntival+Hemoptisis+Nariz+Encias>0,1,0)) # Menorrhagia, vaginal bleeding, subconjunctival bleeding, hemoptysis, epistaxis, or gingivorrhagia
temp$Hem3 = with(temp, ifelse(Melena+Hematemesis>0,1,0)) #Melena or hematemesis
#Pleural Effusion is Derrame_ #check this?? PleuralEffusionLeftRx,  PleuralEffusionRightRx, PulmonaryEdemaRx, PleuralEffusionClinical, PleuralEffusionRightUS or PleuralEffusionLeftUS
#make sure the "Edema" variable includes "Edema, periorbital, facial, limbs, hydrocele, generalized and fovea"
#If (Fever = Yes) AND  ([Tourniquet test] = Yes or [HEM1] = Yes or [HEM2] = Yes or [HEM3] = Yes) AND (Thrombocyto-penia = Yes) 
#  AND ([Hemoconcentration] = Yes or [ElevHematocrit] = Yes  or [PlerualEffusion] = Yes or [Ascites] = Yes or [Hypoprotein-emia] = Yes or [Hypoalbumin-emia] =Yes or [Edema] = Yes), then [DHF] = Yes
temp$Leakage = with(temp, ifelse(Derrame_ + Ascitis + Edema_ + is.Hypoproteinemia + is.Hypoalbuminemia 
                                 + Hemoconcentracion + HematocritoElev > 0, 1, 0))
temp$DHFfinal = with(temp, ((is.torniquete20plus | Hem1 | Hem2 | Hem3) &
                              is.thrombocytopenia &
                              Leakage )) #may want to add is.fever (why does DHF/DSS description sheet not include it?)
#DSS (is.Estrechamiento, is.hypotension, is.pulse_danger, is.coldclammy created in "create_binary_variables")
#Llenado_Capilar is supposedly for if refill rate is > 1 sec but definition wants refill > 2 sec
temp$is.coldclammy = with(temp, ifelse(ExtremidadesFrias + Palidez + Sudoracion + Escalofrio >0 , 1, 0))

temp$DSSfinal = with(temp, ifelse( DHFfinal==1 & (is.Estrechamiento | is.hypotension)
                                   & (is.pulse_danger | is.coldclammy | Llenado_Capilar), 1, 0 ))
temp$myClasificacionFinal = with(temp, ifelse(DSSfinal==1, "DSS",
                                              ifelse(DHFfinal==1,"DHF","DF")))


###############################################################################
################### Add Dx info provided and compare ##########################
###############################################################################

DxCompare = merge(temp, clinical_hospitD_clean[,c("code","ClasificacionFinal")], by="code", all=F)
tabletot(DxCompare, "myClasificacionFinal", "ClasificacionFinal") #very good though not perfect (19 of 1650 are divergent --- 1%)

get_mismatches(colnames(summObs_all),DxDHF_varsH,"in_data","needed") 
#PH has (note: also using info for most of these from Consulta_Externa): 
  # Epistaxis (="Nariz")
  # Encia (="Encias"?)
  # Prueba_torniquete (="Torniquete")
  # Petequia = "PetequiasEspontaneas"
  # Esquimosis (="Equimosis")
  # Edema_ = (Edema+periorbitaria+Facial+Miembros_inferiores+Hidrocele+Generalizado+Fovea > 1) 
  # Darrame_ = ("Derrame_Pleural","Derrame_.pleural_izq_Rx","Derrame_.pleural_der_Rx","Derrame_Clinico","Derecho","Izquierdo") #not a perfect match with definition
#Lab and Lab_Consulta_Externa have "Hematocrito" (for "HematocritoElev")  
#Datos Admision has EFC.Hemoconcentracion. (for "Hemoconcentracion"?) 
get_mismatches(colnames(summObs_all),DxDSS_varsH,"in_data","needed") #
#"Presion_Arterial_Sist","Presion_Arterial_Dias","ExtremidadesFrias",
# is.PoorCapillaryRefill = "Llenado_Capilar"  Done.

