#######################################################################################################
# This function is used by the program R_code_example to generate reconstructed idMS, idMS/MS spectra # 
#######################################################################################################


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

