PairwiseESS<-function(particleDataFrame) {
#Combine doRun$particleDataFrame results and test effective sample size

	library(coda)
	source("/Users/Barb/Desktop/treevo/pkg/R/ESSpairs.R")
	# particleDataFrame can be single or a list of particleDataFrames (1:n)
	x<-particleDataFrame
	
	if(class(x)=="data.frame"){
		warning("ESS on a single run should be high; consider combining several runs")
		data1<-subset(x[which(x[,6]>0),], generation==max(x[,1])) 
		ESS<-effectiveSize(data1[,7:dim(data1)[2]])
		ESSmatrix<-matrix(c(as.numeric(ESS)))
		rownames(ESSmatrix)<-names(ESS)
		colnames(ESSmatrix)<-"ESS"		
	}
	
	if(class(x)=="list"){
		pairs<-as.matrix(ESSpairs(length(x)))
		ESS<-vector("list")
		ESSmatrix<-matrix(nrow=dim(set3[[1]])[2]-6, ncol=dim(pairs)[2])
		rownames(ESSmatrix)<-names(set2[[1]])[7:length(names(set2[[1]]))]
		colnames(ESSmatrix)<-seq(1:dim(ESSmatrix)[2])
		data1<-vector("list")
		for (list in 1:length(x)) {
			data1[[list]]<-subset(x[[list]][which(x[[list]][,6]>0),], generation==max(x[[list]]$generation)) 
		}	
		for (combination in 1:dim(pairs)[2]){
			runNames<-c()
			data2<-c()
			runNames<-c()
			for (run in 1:dim(pairs)[1]){
				if (pairs[run, combination]==1){
					data2<-rbind(data2, data1[[run]])
					runNames<-append(runNames, run)
				}	
			}	
			ESS[[combination]]<-effectiveSize(data2[,7:dim(data2)[2]])
			names(ESS)[[combination]]<-paste(runNames, collapse=".")
			ESSmatrix[,combination]<-c(as.numeric(ESS[[combination]]))
			colnames(ESSmatrix)[combination]<-paste(runNames, collapse=".")
		}
	}
	return(ESSmatrix)
}