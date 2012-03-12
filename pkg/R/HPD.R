HPD<-function(particleDataFrame, percent=0.95, returnData=F){
	summary<-vector("list")
	Ints<-matrix(nrow=(dim(particleDataFrame)[2]-6), ncol=4)
	colnames(Ints)<-c("mean", "sd", "LowerHPD", "UpperHPD")	
	rownames(Ints)<-names(particleDataFrame[7: dim(particleDataFrame)[2]])
	subpDF<-subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame$generation))[7:dim(particleDataFrame)[2]] 
	for(i in 1:dim(subpDF)[2]){
		if(sd(subpDF[,i]) != 0) {
			a<-c()
			b<-c()
			lower<-c()
			upper<-c()
			interv<-c()
			sortedValues<-sort(subpDF[,i])
			numberInTails<-floor(length(sortedValues)*(1-percent))
			for(j in 1:numberInTails){
				a<-append(a, j)
				b<-append(b, length(sortedValues)-j)
				lower<-append(lower, sortedValues[j])
				upper<-append(upper, sortedValues[length(sortedValues)-(numberInTails-j)])
				interv<-append(interv, sortedValues[length(sortedValues)-(numberInTails-j)]-sortedValues[j]) 
			}
			intervalTable<-cbind(a, b, lower, upper, interv)
			#print(intervalTable)
			Ints[i,1]<-mean(subpDF[,i]) #not weighted
			Ints[i,2]<-sd(subpDF[,i])
			Ints[i,3]<-as.numeric(intervalTable[which(intervalTable[,5]==min(intervalTable[,5])),3])
			Ints[i,4]<-as.numeric(intervalTable[which(intervalTable[,5]==min(intervalTable[,5])),4])
			if (returnData){
				low<-as.numeric(intervalTable[which(intervalTable[,5]==min(intervalTable[,5])),1])+1 #add one because you want the next value over
				high<-as.numeric(intervalTable[which(intervalTable[,5]==min(intervalTable[,5])),2])-1
				summary$returnedData<-sortedValues[low:high]
				print(summary)
				if(sd(subpDF[,i]) != 0){
					summary$returnedData<-NA
				}
			}
				
		}
	}
	if(returnData){
		summary$summary<-Ints
		return(summary)
	}
	else{return(Ints)}
}
