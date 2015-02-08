


### function takes abundance data as input and produces histograms of MZ values
graph_MZ_frequencies = function(data1, data2, dname1, dname2, fileprefix) {
  
  MZ_Nums_reshaped1 = get_MZ_vals(data1)
  MZ_Nums_reshaped2 = get_MZ_vals(data2)
  
  
  #mybreaks = seq(min(MZ_Nums_reshaped1),max(MZ_Nums_reshaped1),10)
  #histogram(MZ_Nums_reshaped1, nint=100) #overview
  #histogram(MZ_Nums_reshaped1, xlim=c(300,400), nint=10000) #zoomed
  #histogram(MZ_Nums_reshaped1, xlim=c(330,380), nint=10000) #zoomed
  
  ### kernel density, overview
  png(paste(resultsDir,fileprefix,"MZ_density_full.png", sep=""))
  plot(density(MZ_Nums_reshaped1, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
       col="blue", xlab="", main="Kernel Density of MZ Values (full range)")
  par(new=TRUE)
  plot(density(MZ_Nums_reshaped2, bw=.00001), xlim=c(100,1700), ylim=c(0,300),
       col="red", xlab="", main="") 
  legend("topright",legend=c(dname1, dname2), col=c("blue","red"),lty=1)
  dev.off()
  
  ### kernel density, zoomed in
  
  png(paste(resultsDir,fileprefix,"MZ_density_100.png", sep=""))
  plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
       col="blue", main="Kernel Density of MZ Values", xlab="") #kernal density
  lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(100,200), ylim=c(0,3000),
        col="red", main="", xlab="") #kernal density
  legend("topleft",legend=c(dname1, dname2), col=c("blue","red"),lty=1)
  dev.off()
  
  plot_kdensity = function(mz_start, mz_end, ylimit, fileprefix) {
    png(paste(resultsDir,fileprefix,"MZ_density_",mz_start,".png", sep=""))
    plot(density(MZ_Nums_reshaped1, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
         col="blue", main="", xlab="") #kernal density
    lines(density(MZ_Nums_reshaped2, bw=.000001),xlim=c(mz_start,mz_end), ylim=c(0,ylimit),
          col="red", main="", xlab="") #kernal density
    dev.off()
  }
  plot_kdensity(mz_start=200, mz_end=300, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=300, mz_end=400, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=400, mz_end=500, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=500, mz_end=600, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=600, mz_end=700, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=700, mz_end=800, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=800, mz_end=900, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=900, mz_end=1000, ylimit=3000, fileprefix)
  plot_kdensity(mz_start=1000, mz_end=1100, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1100, mz_end=1200, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1200, mz_end=1300, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1300, mz_end=1400, ylimit=1000, fileprefix)
  plot_kdensity(mz_start=1400, mz_end=1500, ylimit=1000, fileprefix)
  
}


####### Compare frequences of MZ values across datasets (histograms) ########

graph_MZ_frequencies(respD1_filter50n, respD2_filter50n, "Round1_filter50", "Round2_filter50", "f50_")
graph_MZ_frequencies(respD1_filter10n, respD2_filter10n, "Round1_filter10", "Round2_filter10", "f10_")


######### Other checks and descriptives of abundance data #############

#missingness (exclude "LCMS_run","Study", and "code" columns)
missingsD = sapply(respD1_filter50n[grep("MZ_",colnames(respD1_filter50n))], function(x) sum(is.na(x))) 
#missingsD[which(missingsD==0)] #examine the compounds that are never missing
fileprefix = "R1_filter50"
png(paste(resultsDir,fileprefix,"_missingMZ.png", sep=""))
histogram(missingsD, xlim=c(0,30), xlab="number of observations with missing values within each MZ value",
          main="LC-MS Run 1 with 50% filter")
dev.off()
summary(missingsD)
table(missingsD)
sum(missingsD==0)/dim(respD1_filter50n)[2] #fraction of compounds that have no missings
sum(missingsD<2)/dim(respD1_filter50n)[2] #fraction of compounds that have 0 or 1 missings
quantile(missingsD, c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99)/100) #more percentiles
#do we have any zeros for abundance?
zerosD = sapply(resp_D1D2, function(x) sum(x==0, na.rm=TRUE))
sum(zerosD!=0, na.rm=TRUE) #no zeros
