CredInt<-function(particleDataFrame, percent=0.95) { 
	generation<-NULL #to appease R CMD CHECK with subset
	PercentTail<-(1-percent)/2
	Ints<-matrix(nrow=(dim(particleDataFrame)[2]-6), ncol=4)
	colnames(Ints)<-c("mean", "sd", "LowerCI", "UpperCI")
	rownames(Ints)<-names(particleDataFrame[7: dim(particleDataFrame)[2]])
	subpDF<-subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame$generation))[7:dim(particleDataFrame)[2]] 
	for(i in 1:dim(subpDF)[2]){
		if(sd(subpDF[,i]) != 0) {
			Ints[i,1]<-mean(subpDF[,i]) #not weighted
			Ints[i,2]<-sd(subpDF[,i])
			Ints[i,3]<-quantile(subpDF[,i], probs=0+PercentTail)
			Ints[i,4]<-quantile(subpDF[,i], probs=1-PercentTail)
		}
	}
	return(Ints)
}