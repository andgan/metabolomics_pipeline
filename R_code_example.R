########################################################################################
# This program describes feature detection, quality control and spectra generation for #
# PIVUS data. Real data are not shared, but the skelthon of the program can be applied #
# to every UPLC-MS/MS data. Three data sources are used: 1. Raw data (.CDF format)     #
# 2. ID to identify samples from same individual (XXXnamesXXX) 3. Phenotypic data to   #
# identify variables of unwanted variabiloty that need to be adjusted.                 #
########################################################################################

### SUMMARY ####
# MODULE 1:
#  1.Peaks detection
#  2.Peaks alignment
#  3. Peaks grouping
#  4. Filling missing features
#
# MODULE 2:
#  0. Separate and save idMS and idMS/MS data
#  1. Log2 transform
#  2. Calculate total intesity for each sample and check for outlier
#  3. ANOVA-type normalization
#  4. Average between duplicates
# 
# CREATE idMS and idMS/MS spectra
#
# Specfun.R function
######################



### LOAD LIBRARIES ###
library(xcms)
library(pcaMethods)


#########################
####### MODULE 1 ########
### FEATURE DETECTION ###
#########################

# Get raw data
raw_pivus <-list.files("path_to_raw_data", pattern=".CDF", recursive = FALSE, full.names=TRUE)

### 1.Peaks detection ###
print("Peak detection")
xset<-xcmsSet(raw_pivus, method = "centWave", ppm=25 , peakwidth=c(2:15), snthresh=8, mzCenterFun="wMean", integrate=2, mzdiff=0.05 , prefilter=c(1,5))
## ! Use 'nSlaves' option to make this parallel. To use nSlaves you need to load library(Rmpi) ##

#########
## If needed is possible to use just one sample to check the peak detection quality ##
## check quality of peak detection in 'peak_detection.pdf'
one_sample <- xcmsRaw("path_to_raw_data/one_sample.CDF")
pdf("peak_detection.pdf")
f1 <- findPeaks(one_sample,"path_to_raw_data/one_sample.CDF", method = "centWave", ppm=25 , peakwidth=c(2:15), snthresh=8, mzCenterFun="wMean", integrate=2, mzdiff=0.05 , prefilter=c(1,5), sleep=0.0001, fitgauss=T)
dev.off()
#########

### 2.Peaks alignment ###
xset2 <- retcor(xset,  method="obiwarp",plottype="deviation") 

# Plot a function of the retention time deviation for each sample
rmsd<-sapply(1:length(xset2@rt$corrected),function(x) {sqrt(sum((xset2@rt$raw[[x]]-xset2@rt$corrected[[x]])**2))})
plot(rmsd)

### 3. Peaks grouping ###
# Group retention-time-corrected peaks and print plots for visual inspection of grouping quality in 'grouping.pdf'
pdf("grouping.pdf")
xset3 <- group(xset2, bw=2, minfrac=0.03, max= 100, mzwid=0.01, sleep=0.0001)
dev.off()

### 4. Filling missing features ####
xset4 <- fillPeaks.chrom(xset3)

# Get features x samples matrix
xset5 <- groupval(xset4, value="into")
rownames(xset5) <- groupnames(xset4,mzdec=3, rtdec=3)


#######################
####### MODULE 2 ######
### QUALITY CONTROL ###
#######################

### 0. Here idMS are separated from idMS/MS, only idMS are kept for future analysis ###
# Save only idMS samples (they end with '01')
d_cdf1 <- xset5[,substr(colnames(xset5),nchar(colnames(xset5))-1,nchar(colnames(xset5)))=="01"]

# Save to create idMS spectra later
d_cdf1a <- t(d_cdf1)
d_cdf1ab <- cbind(row.names(d_cdf1a),d_cdf1a)
colnames(d_cdf1ab)[1] <- "Sample.ID"
colnames(d_cdf1ab) <- gsub("M","",colnames(d_cdf1ab))
colnames(d_cdf1ab) <- gsub("T","_",colnames(d_cdf1ab))

write.csv(d_cdf1ab,file=paste(path,"MSdata.csv",sep=""), quote=F, row.names=F)
 

# Save only idMS/MS samples (they end with '02')
d_cdf2 <- xset5[,substr(colnames(xset5),nchar(colnames(xset5))-1,nchar(colnames(xset5)))=="02"]

# Save to create idMS/MS spectra later
d_cdf2a <- t(d_cdf2)
d_cdf2ab <- cbind(row.names(d_cdf2a),d_cdf2a)
colnames(d_cdf2ab)[1] <- "Sample.ID"
colnames(d_cdf2ab) <- gsub("M","",colnames(d_cdf2ab))
colnames(d_cdf2ab) <- gsub("T","_",colnames(d_cdf2ab))

write.csv(d_cdf2ab,file=paste(path,"MSMSdata.csv",sep=""), quote=F, row.names=F)


### 1. Log2 transform ### 

d_cdf1_n <-  apply(d_cdf1,2,as.numeric)
# substitute 0 with min(d_cdf1), if intensity is 0, then it gets the lower possible value
d_cdf1_n[d_cdf1_n == 0] <- min(d_cdf1_n[d_cdf1_n != 0])
Ld_cdf1_n <- log2(d_cdf1_n)

### 2. Calculate total intesity for each sample and check for outlier ###

TIC<-colSums(Ld_cdf1_n, na.rm=T)

# Plot to identify outliers
pdf("outliers.pdf")
par(mfrow=c(1,2))
boxplot(TIC, ylab="Extracted peaks area")
hist(TIC, breaks=100, main="", xlab="Extracted peaks area")
dev.off()

# Where XXXXX is the TIC value for outliers exclusion
Ld_cdf1_n2 <- Ld_cdf1_n[,!(colnames(Ld_cdf1_n) %in% c(colnames(Ld_cdf1_n)[TIC<XXXXX]))]

### 3. ANOVA-type normalization ###
# Where XXXnamesXXX is a vector with the id if the samples (so that samples from the same individual have the same id)
metabo_num <-  apply(t(Ld_cdf1_n2),2,as.numeric)
pc <- pca(metabo_num,method="ppca", nPcs=4)

### ! Identify variables highly associated with PC1 and PC2, you might want to remove their effect from the data ##

# For example, we run a linear regression for each feature, adjusting by season and storage time 
lmres <- function(x){resid(lm(as.numeric(x)~as.numeric(season)+as.numeric(storage), na.action=na.exclude))}
Ld_cdf1_n3 <- t(sapply(data.frame(t(Ld_cdf1_n2)),lmres))
colnames(Ld_cdf1_n3) <- XXXnamesXXX
rownames(Ld_cdf1_n3) <- rownames(d)


### 4. Average between duplicates ###

# Loop, select two injections from the same individuals, average them and calculate correlation
D_S_C <- NULL
CR <- NULL
for (j in unique(substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2)))
	{
		# Select two samples form same individual
		d_s <- Ld_cdf1_n3[,substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2)==j]
		
		if (!is.null(dim(d_s)[2]))
		{
			# Averge
			d_s_c <- rowMeans(d_s)
			# Check correlations
			cr <- cor(d_s, method="spearman", use="complete.obs")[upper.tri(cor(d_s, method="spearman", use="complete.obs"))]
			# Transform correlation
			zcr <- 0.5*log((1+cr)/(1-cr))
		}else
		{
			d_s_c <- d_s
			zcr <- NA
		}
		D_S_C <- rbind(D_S_C,d_s_c)
		CR <- c(CR,zcr)
	}
	
rownames(D_S_C) <- unique(substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2))

# Make it a data.frame
Ld_cdf1_n4 <- data.frame(cbind(rownames(D_S_C),CR, D_S_C), stringsAsFactors=F)
colnames(Ld_cdf1_n4)[1:2] <- c("id","zcorr")

write.table(Ld_cdf1_n4,file="final_dataset.txt", row.names=F, col.names=T, quote=F, sep="\t")


#######################################
### CREATE idMS and idMS/MS spectra ###
#######################################


t <- strsplit(gsub("M","", colnames(Ld_cdf1_n4_res)[c(-1,-2)]),"T")
M <- as.numeric(sapply(t,"[",1))
T <- as.numeric(sapply(t,"[",2))

# Write the list of all the features with a spectra
write.table(paste(formatC(M,format="f",digits=3),"_",formatC(T,format="f",digits=3), sep=""), file="features_for_spectra.csv",row.names=F,quote=F, col.names=F, sep=",")


# idMSMS and DDA Parameter selection 
maxlabel<-12 				# maximum number of features labeled on spectrum plots (10 default)
rt_dev<-2					# retention time window (seconds) on either side of feature rt for idMSMS spectral reconstruction 
library<-"pivus_small"			# name your spectral library - ensure quotes on either side
cor_filter<-"rval"  			# can be either "rval" or "pval" - ensure quotes on either side
rcut<-0.5					# r value cutoff for idMSMS spectral reconstruction - MS or MSMS
pcut<-0.01					# optional pval cutoff - select either rcut or pcut for "cor_filter"; caution: pval can reflect either positive of negative correlational relationships		

import_feat_list<-"TRUE"		# logical, should R import csv file with one column of feature names, no header
FileName<-"features_for_spectra.csv"		###file name of feature list for spectra plotting

# read isMS and isMS/MS data saved at the beginning of module 2
if(import_feat_list<-"TRUE") {
sig_feat_list<-as.vector(read.csv(FileName, header=FALSE, check.names=FALSE)[,1])
MSdata<-read.csv(file="MSdata.csv", header=TRUE, check.names=FALSE)
MSMSdata<-read.csv(file="MSMSdata.csv", header=TRUE, check.names=FALSE)
feature<-names(MSdata)[2:dim(MSdata)[2]]
}

## ! remember to check that all the samples are both in idMS and idMS/MS, otherwise keep only those in common ##

sig_list<-unique(sig_feat_list)
mass_Rt<-data.frame(sapply(sig_list, strsplit, sig_list, split="_"))
sig_feat_mz<-as.numeric(t(mass_Rt[1,]))
sig_feat_rt<-as.numeric(t(mass_Rt[2,]))
sig_feat<-data.frame(sig_feat_mz, sig_feat_rt)
names(sig_feat)<-c("mz", "rt")
row.names(sig_feat)<-sig_list
libName<-paste(library, ".msp", sep="")
feature<-names(MSdata[2:dim(MSdata)[2]])
mass_Rt_full<-data.frame(sapply(feature, strsplit, feature, split="_"), check.names=FALSE)
mass_Rt_full[1:2,1:4]
mz<-as.numeric(t(mass_Rt_full[1,]))
rt<-as.numeric(t(mass_Rt_full[2,]))

# Load the specfun function, see end of the program
source("specfun.R")

# Loop each feature
lapply(1:length(sig_list), specfun)



################################
##########SPECFUN.R ############
### SAVE IT AS .R FILE #########
### IT IS USED TO GENERATE #####
### idMS and isMS/MS spectra ###
################################

## begin specfun
specfun<- function (i) {
targetmz<-sig_feat[i,"mz"]
targetrt<-sig_feat[i,"rt"]

##subset MSdata by retention time similarity to target feature for idMS recronstruction
## get all the features in a certain rt window (define by rt_dev)
rt_index<-intersect(which(rt>(targetrt-rt_dev)), which(rt<(targetrt+rt_dev)))
if(length(rt_index)>1) MS<-MSdata2[,rt_index+1] else MS<-as.data.frame(MSdata2[,rt_index+1])
if(length(rt_index)>1) MSMS<-MSMSdata2[,rt_index+1] else MSMS<-as.data.frame(MSMSdata2[,rt_index+1])
names(MS)<-names(MSdata2)[rt_index+1]
names(MSMS)<-names(MSMSdata2)[rt_index+1]




# Calculate the correlation between the feature of interest and the features in the same rt window
corfun<-function (x) {cor.test(MS[,sig_list[i]], MS[,x], exact=TRUE, na.rm=TRUE)$p.value}
corfun2<-function (x) {cor.test(MS[,sig_list[i]], MS[,x], exact=TRUE, na.rm=TRUE)$estimate}
if(cor_filter=="pval") pval<-as.numeric(sapply(1:dim(MS)[2], corfun))
if(cor_filter=="rval")  rval<-as.numeric(sapply(1:dim(MS)[2], corfun2))
if (cor_filter=="rval") idMS_index<-which(rval>rcut) else idMS_index<-which(pval<pcut)
if(length(idMS_index) > 1) {
	idMS<-MS[,idMS_index] 
	names(idMS)<-names(MS)[idMS_index]
	idMS_mass_Rt<-data.frame(sapply(names(idMS), strsplit, names(idMS), split="_"))
	idMS_sig_feat_mz<-as.numeric(t(idMS_mass_Rt[1,]))
	idMS_sig_feat_rt<-as.numeric(t(idMS_mass_Rt[2,]))
	idMS_abundance<-colMeans(idMS)
	idMS_spec<-data.frame(idMS_sig_feat_mz, idMS_sig_feat_rt, idMS_abundance, rval[idMS_index])
	names(idMS_spec)<-c("mz", "rt", "intensity", "rval")
	idMS_spec<-idMS_spec[order(idMS_spec[,"intensity"], decreasing=TRUE),]
	row.names(idMS_spec)<-names(idMS)
	} else
   {idMS_spec<-data.frame(targetmz, targetrt, 100, 1)
	row.names(idMS_spec)<-row.names(sig_feat)[i]
	names(idMS_spec)<-c("mz", "rt", "intensity", "rval")}

##write idMS spectra to msp library file
if(length(idMS_index) > 1) {
idMS_npeaks<-dim(idMS_spec)[1]
write(paste("Name: ", row.names(sig_feat)[i], " MS spectrum", sep=""), libName, append=TRUE)
write(paste("SYNON: $:00in-source", sep=""), libName, append=TRUE)
write(paste("SYNON: $:04 0", sep=""), libName, append=TRUE)
write(paste("Comment: r-value cutoff = ", rcut, sep=""), libName, append=TRUE)
write(paste("Comment: rt window = ", rt_dev, " sec" , sep=""), libName, append=TRUE)
write(paste("Comment:", round(idMS_spec[1,"mz"], digits=4),  round(idMS_spec[2,"mz"], digits=4),  round(idMS_spec[3,"mz"], digits=4),  round(idMS_spec[4,"mz"], digits=4),round(idMS_spec[5,"mz"], digits=4), round(idMS_spec[6,"mz"], digits=4), round(idMS_spec[7,"mz"], digits=4)), file=libName, append= TRUE)
write(paste("Num Peaks:", idMS_npeaks), file=libName, append= TRUE)
idMS_writefun<- function (j) {
	ion<- paste(round(idMS_spec[j, "mz"], digits=4), " ", round(idMS_spec[j, "intensity"]), ";", sep="")
		write(ion, file=libName, append= TRUE) } ##end writefun
	lapply(1:dim(idMS_spec)[1], idMS_writefun)
write("", file=libName, append= TRUE)
} 



#idMSMS below:
corfun<-function (x) {cor.test(MS[,sig_list[i]], MSMS[,x], exact=TRUE, na.rm=TRUE)$p.value}
corfun2<-function (x) {cor.test(MS[,sig_list[i]], MSMS[,x], exact=TRUE, na.rm=TRUE)$estimate}
pval<-as.numeric(sapply(1:dim(MSMS)[2], corfun))
rval<-as.numeric(sapply(1:dim(MSMS)[2], corfun2))
if (cor_filter=="rval") idMSMS_index<-which(rval>rcut) else idMSMS_index<-which(pval<pcut)

if(length(idMSMS_index) == 1) {  idMSMS<-data.frame(MSMS[,idMSMS_index], check.names=FALSE)
		idMSMS_mass_Rt<-data.frame(sapply(names(MSMS)[idMSMS_index], strsplit, names(idMSMS), split="_"), check.names=FALSE)	
		idMSMSrt<-as.numeric(as.vector(idMSMS_mass_Rt[2,1]))
		idMSMSmz<-as.numeric(as.vector(idMSMS_mass_Rt[1,1]))
		idMSMS_spec<-data.frame(idMSMSrt, idMSMSmz,  colMeans(idMSMS), rval[idMSMS_index])
		names(idMSMS_spec)<-c("mz", "rt", "intensity", "rval")		
		row.names(idMSMS_spec)<-names(MSMS)[idMSMS_index]} else
			if  (length(idMSMS_index) > 1)  { idMSMS<-data.frame(MSMS[,idMSMS_index], check.names=FALSE)
		idMSMS_mass_Rt<-data.frame(sapply(names(MSMS)[idMSMS_index], strsplit, names(idMSMS), split="_"), check.names=FALSE)	
#		idMSMS_mass_Rt<-data.frame(sapply(names(idMSMS), strsplit, names(idMSMS), split="_"))
		idMSMS_sig_feat_mz<-as.numeric(t(idMSMS_mass_Rt[1,]))
		idMSMS_sig_feat_rt<-as.numeric(t(idMSMS_mass_Rt[2,]))
		idMSMS_abundance<-colMeans(idMSMS)
		idMSMS_spec<-data.frame(idMSMS_sig_feat_mz, idMSMS_sig_feat_rt, idMSMS_abundance, rval[idMSMS_index], check.names=FALSE)
		names(idMSMS_spec)<-c("mz", "rt", "intensity", "rval")		
		idMSMS_spec<-idMSMS_spec[order(idMSMS_spec[,"intensity"], decreasing=TRUE),]	} 			



##write idMSMS spectra to msp library file
if(length(idMSMS_index) > 0) {
idMSMS_npeaks<-dim(idMSMS_spec)[1]
write(paste("Name: ", row.names(sig_feat)[i], " idMSMS spectrum", sep=""), libName, append=TRUE)
write(paste("SYNON: $:00in-source", sep=""), libName, append=TRUE)
write(paste("SYNON: $:04 0", sep=""), libName, append=TRUE)
write(paste("Comment: r-value cutoff = ", rcut, sep=""), libName, append=TRUE)
write(paste("Comment: rt window = ", rt_dev, " sec" , sep=""), libName, append=TRUE)
write(paste("Comment:", round(idMSMS_spec[1,"mz"], digits=4),  round(idMSMS_spec[2,"mz"], digits=4),  round(idMSMS_spec[3,"mz"], digits=4),  round(idMSMS_spec[4,"mz"], digits=4),  
	round(idMSMS_spec[5,"mz"], digits=4), round(idMSMS_spec[6,"mz"], digits=4), round(idMSMS_spec[7,"mz"], digits=4)), file=libName, append= TRUE)
write(paste("Num Peaks:", idMSMS_npeaks), file=libName, append= TRUE)
idMSMS_writefun<- function (j) {
	ion<- paste(round(idMSMS_spec[j, "mz"], digits=4), " ", round(idMSMS_spec[j, "intensity"]), ";", sep="")
		write(ion, file=libName, append= TRUE) } ##end writefun
	lapply(1:dim(idMSMS_spec)[1], idMSMS_writefun)
write("", file=libName, append= TRUE)
} 




##plot reconstructed MS spectrum and selected DDA spectrum
#maxx<-if(length(idMS_index)==0) {1.2*max(idMSMS_spec[,"mz"])} else {1.2*max(idMSMS_spec[,"mz"], idMS_spec[,"mz"])}
if(any(ls()=="idMSMS_spec")==TRUE) maxx<-round(1.2*max(idMSMS_spec[,"mz"], idMS_spec[,"mz"])) else
	maxx<-round(1.2*max(idMS_spec[,"mz"])) 

parent<-sig_feat_list[i]
pdf(file=paste(row.names(sig_feat)[i], ".pdf", sep=""))
par(mfrow=c(2,1), pty="m", mar=c(0,4,3,0), omi=c(1,0.2,0,0.2))

if(length(idMS_index) > 1) {
plot(idMS_spec[,"mz"], idMS_spec[,"intensity"], type="h", xaxt="n", ann=FALSE, yaxs="i", xlim=c(40, maxx), ylim=c(0, max(idMS_spec[,"intensity"])*1.3))
title(main = paste("MS spectrum:", targetmz, "@", targetrt, "sec"))
text(idMS_spec[1:maxlabel,"mz"], idMS_spec[1:maxlabel,"intensity"], label=round(idMS_spec[1:maxlabel,"mz"], digits=3), pos=3, srt=90, offset=1, cex=0.85*(idMS_spec[,"rval"])^1.5)
}  else {
plot(50, 100, type="n", xaxt="n", ann=FALSE, yaxs="i", xlim=c(40, maxx), ylim=c(0, 100))
text(mean(c(40, maxx)), 50, label="no correlating fragments", pos=1, offset=0.2)
title(main = paste("MS spectrum:", targetmz, "@", targetrt, "sec"))
}


if(length(idMSMS_index) > 0) {
maxx<-round(1.2*max(idMSMS_spec[,"mz"], idMS_spec[,"mz"]))
plot(idMSMS_spec[,"mz"], idMSMS_spec[,"intensity"], type="h", ann=FALSE, yaxs="i", xlim=c(40, maxx), ylim=c(0, max(idMSMS_spec[,"intensity"])*1.3))
title(main = paste("idMSMS spectrum:", targetmz, "@", targetrt, "sec"))
text(idMSMS_spec[1:maxlabel,"mz"], idMSMS_spec[1:maxlabel,"intensity"], label=round(idMSMS_spec[1:maxlabel,"mz"], digits=3), pos=3, srt=90, offset=1, cex=0.85*(idMSMS_spec[,"rval"])^1.5)
}  else {
plot(50, 100, type="n", xaxt="n", ann=FALSE, yaxs="i", xlim=c(40, maxx), ylim=c(0, 100))

text(mean(c(40, maxx)), 50, label="no correlating idMSMS fragments", pos=1, offset=0.2)
title(main = paste("idMSMS spectrum:", targetmz, "@", targetrt, "sec"))
}
mtext("m/z", side=1, line=2.5, outer=TRUE, font=2)
mtext("intensity", side=2, line=-0.5, outer=TRUE, font=2)
dev.off()
gc(verbose=FALSE)
}  ##end specfun

