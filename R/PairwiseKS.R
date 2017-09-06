#' Pairwise Kolmogorov-Smirnov test
#' 
#' This function calculates Kolmogorov-Smirnov on results.
#' 
#' 

#' @param particleDataFrame List of particleDataFrames from separate runs

#' @return Returns a matrix with Kolmogorov-Smirnov values of all pairwise runs

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished


#' @name intrinsicModels
#' @rdname intrinsicModels
#' @export
PairwiseKS<-function(particleDataFrame) {
#Combine doRun$particleDataFrame results and performs Kolmogorov-Smirnov tests
# particleDataFrame should be a list of particleDataFrames (1:n)

	#source("/Users/Barb/Desktop/treevo/pkg/R/pairings.R")
	
		
	if(class(particleDataFrame)=="data.frame"){
		warning("KS test requires comparison of two runs")
	}
	
	if(class(particleDataFrame)=="list"){
		x<-vector("list")
		for (list in 1:length(particleDataFrame)) {
			x[[list]]<-subset(particleDataFrame[[list]][which(particleDataFrame[[list]][,6]>0),], particleDataFrame[[list]]$generation==max(particleDataFrame[[list]]$generation)) 
		}	
	
		KSMatrixList<-vector("list", dim(x[[1]])[2]-6)
		names(KSMatrixList)<-names(x[[1]][7:dim(x[[1]])[2]])
		for (m in 1:length(KSMatrixList)){
		param<-6+m
		KSMatrix<-matrix(nrow=length(x), ncol=length(x))
		rownames(KSMatrix)<-paste("run", 1:length(x), sep="")
		colnames(KSMatrix)<-paste("run", 1:length(x), sep="")
			for (i in 1:dim(KSMatrix)[1]){
				for (j in 1:dim(KSMatrix)[2]){
					KSMatrix[i, j]<-suppressWarnings(round(ks.test(x[[i]][, param],x[[j]][, param])$p.value, digits=4)) #make diag
				}
			}
		KSMatrixList[[m]]<-KSMatrix
		}
	}
	return(KSMatrixList)
}
