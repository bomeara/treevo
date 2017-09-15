#' Pairwise Effective Sample Size
#'
#' This function calculates Effective Sample Size (ESS) on results.  Performs the best when results are
#' from multiple runs.
#'
#'

#' @param particleDataFrame particleDataFrame can be a single data frame or a
#' list of data frames

#' @return Returns a matrix with ESS values of all pairwise runs

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#'
#' breakThisExample

#' @name PairwiseESS
#' @rdname PairwiseESS
#' @export
PairwiseESS<-function(particleDataFrame) {
#Combine doRun$particleDataFrame results and test effective sample size
  ESSmatrix <- NA
	#library(coda)
	# particleDataFrame can be single or a list of particleDataFrames (1:n)
	x<-particleDataFrame

	if(class(x)=="data.frame"){
		warning("ESS on a single run should be high; consider combining several runs")
		data1<-subset(x, x$generation==max(x[,1]))
		data1<-data1[which(data1$weight>0),]

		ESS<-coda::effectiveSize(data1[,7:dim(data1)[2]])
		ESSmatrix<-matrix(c(as.numeric(ESS)))
		rownames(ESSmatrix)<-names(ESS)
		colnames(ESSmatrix)<-"ESS"
	}

	if(class(x)=="list"){
		pairs<-as.matrix(pairings(length(x)))
		ESS<-vector("list")
		ESSmatrix<-matrix(nrow=dim(x[[1]])[2]-6, ncol=dim(pairs)[2])
		rownames(ESSmatrix)<-names(x[[1]])[7:length(names(x[[1]]))]
		colnames(ESSmatrix)<-seq(1:dim(ESSmatrix)[2])
		data1<-vector("list")
		for (list in 1:length(x)) {
			data1[[list]]<-subset(x[[list]][which(x[[list]][,6]>0),], x[[list]]$generation==max(x[[list]]$generation))
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
