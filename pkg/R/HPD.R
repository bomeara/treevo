HPD<-function(particleDataFrame, percent=0.95, returnData=F){
	library(coda, quietly=T)
	summary<-vector("list")
	Ints<-matrix(nrow=(dim(particleDataFrame)[2]-6), ncol=4)
	colnames(Ints)<-c("mean", "sd", paste("LowerHPD_", percent, sep=""), paste("UpperHPD_", percent, sep=""))	
	rownames(Ints)<-names(particleDataFrame[7: dim(particleDataFrame)[2]])
	subpDF<-subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame$generation))[7:dim(particleDataFrame)[2]] 
	for(i in 1:dim(subpDF)[2]){
		if(sd(subpDF[,i]) != 0) {

			Ints[i,1]<-mean(subpDF[,i]) #not weighted
			Ints[i,2]<-sd(subpDF[,i])
			Ints[i,3]<-HPDinterval(as.mcmc(subpDF[i]), prob=percent)[1] #returns lower HPD		
			Ints[i,4]<-HPDinterval(as.mcmc(subpDF[i]), prob=percent)[2] #returns upper HPD	
			if(returnData){
				summary[[i]]<-subpDF[,i][-c(which(subpDF[,i]<Ints[i,3]), which(subpDF[,i]>Ints[i,4]))]		
			}		
		}
	}
	if(returnData){
		summary$summary<-Ints
		return(summary)
	}
	else{return(Ints)}
}