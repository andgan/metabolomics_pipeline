##############################################################################################
# This program describes feature detection, quality control and spectra generation for       #
# PIVUS data. Real data are not shared, but the skelthon of the program can be applied       #
# to every UPLC-MS/MS data. Three data sources are used:                                     # 
# 1. Raw data (.CDF format)                                                                  #
# 2. ID to identify samples from same individual (XXXnames_injXXX): it should be an id for   #
#    each individual with this format: id_inj where inj is the injection number (1,2).       #    
#    For example: 123_1 refers to individual 123, injection one;                             # 
#								  123_2 refers to individual 123, injection 2.                               #
# 3. Phenotypic data to identify variables of unwanted variabilty that need to be adjusted.  #
##############################################################################################

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
#  5. Feature-level correlation between duplicates
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
# Where XXXnames_injXXX is a vector with the id of the individuals and the injection (replicate) number; (so that samples from the same individual have the same id, but different injection)
metabo_num <-  apply(t(Ld_cdf1_n2),2,as.numeric)
pc <- pca(metabo_num,method="ppca", nPcs=4)

### ! Identify variables highly associated with PC1 and PC2, you might want to remove their effect from the data ##

# For example, we run a linear regression for each feature, adjusting by season and storage time 
lmres <- function(x){resid(lm(as.numeric(x)~as.numeric(season)+as.numeric(storage), na.action=na.exclude))}
Ld_cdf1_n3 <- t(sapply(data.frame(t(Ld_cdf1_n2)),lmres))
colnames(Ld_cdf1_n3) <- XXXnames_injXXX 
rownames(Ld_cdf1_n3) <- rownames(d)


### 4. Average between duplicates ###

# Loop, select two injections from the same individuals, average them
D_S_C <- NULL
for (j in unique(substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2)))
	{
		# Select two samples form same individual
		d_s <- Ld_cdf1_n3[,substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2)==j]
		
		if (!is.null(dim(d_s)[2]))
		{
			d_s_c <- rowMeans(d_s)
		}
		else
		{
			d_s_c <- d_s
		}
		D_S_C <- rbind(D_S_C,d_s_c)
	}
	
rownames(D_S_C) <- unique(substr(colnames(Ld_cdf1_n3),1,nchar(colnames(Ld_cdf1_n3))-2))

# Make it a data.frame
Ld_cdf1_n4 <- data.frame(cbind(rownames(D_S_C), D_S_C), stringsAsFactors=F)
colnames(Ld_cdf1_n4)[1] <- c("id")


### 5. Feature-level correlation between duplicates ###
## If not all the samples have duplicates, select only those samples with duplicates
CRM <- NULL
CRT <- NULL
for (k in 1:nrow(Ld_cdf1_n3))
{
	cc <- Ld_cdf1_n3[k,]
	cc1 <- cc[substr(colnames(Ld_cdf1_n3),nchar(colnames(Ld_cdf1_n3)),nchar(colnames(Ld_cdf1_n3)))=="1"]
	cc1o <- cc1[order(names(cc1))]
	cc2 <- cc[substr(colnames(Ld_cdf1_n3),nchar(colnames(Ld_cdf1_n3)),nchar(colnames(Ld_cdf1_n3)))=="2"]
	cc2o <- cc2[order(names(cc2))]
	
	CR <- cor(cc2o,cc1o)
	CR <- cor.test(as.numeric(cc2o),as.numeric(cc1o), alternative="greater")	
	CRM <- c(CRM,CR$estimate)
	CRT <- c(CRT,CR$p.value)
}


## Plot distribution between correlations
## Plot the line corresponding to 5% FDR, features below this level have too low correlation and are excluded
pdf("dist_corr_features.pdf")
p <- density(CRM)
plot(p, xlab="Pearson Correlation", main="Distribution correlation between features")
abline(v=CRM[which.min(abs(p.adjust(CRT,method="BH")-0.05))])
dev.off()

### EXCLUDE FEATURES WITH POOR CORRELATION ###

Ld_cdf1_n5 <- Ld_cdf1_n4[,c(TRUE,as.logical(CRM > CRM[which.min(abs(p.adjust(CRT,method="BH")-0.05))]))]

write.table(Ld_cdf1_n5,file="final_dataset.txt", row.names=F, col.names=T, quote=F, sep="\t")


#######################################
### CREATE idMS and idMS/MS spectra ###
#######################################


t <- strsplit(gsub("M","", colnames(Ld_cdf1_n5)[c(-1)]),"T")
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

