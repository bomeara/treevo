#' Bayesian Credible Interval
#' 
#' This function calculates credible interval for each free parameter.
#' 
#' 

#' @param particleDataFrame A \code{particleDataFrame} object, as found among the output from \code{\link{doRun}} functions.

#' @param percent Probability content of the highest probability density.

#' @return Returns a matrix with weighted mean, standard deviation, upper and lower credible
#' intervals for each free parameter.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#' 
#' data(simRunExample)
#' credibleInt(results$particleDataFrame)
#' 

# was named CredInt originally

#' @name credibleInt
#' @rdname credibleInt
#' @export
credibleInt<-function(particleDataFrame, percent=0.95) {
	# ugh ugh
	#generation<-NULL #to appease R CMD CHECK
	# yes??? I think this is right, not sure
	generation<-particleDataFrame$generation
	#
	PercentTail<-(1-percent)/2
	Ints<-matrix(nrow=(dim(particleDataFrame)[2]-6), ncol=4)
	colnames(Ints)<-c("mean", "sd", "LowerCI", "UpperCI")
	rownames(Ints)<-names(particleDataFrame[7: dim(particleDataFrame)[2]])
	subpDF<-subset(particleDataFrame[which(particleDataFrame$weight>0),],
		generation==max(particleDataFrame$generation),drop=FALSE)[7:dim(particleDataFrame)[2]]
	for(i in 1:dim(subpDF)[2]){
		if(length(subpDF[,i])>1){
			if(sd(subpDF[,i]) != 0) {
				Ints[i,1]<-mean(subpDF[,i]) #not weighted
				Ints[i,2]<-sd(subpDF[,i])
				Ints[i,3]<-quantile(subpDF[,i], probs=0+PercentTail)
				Ints[i,4]<-quantile(subpDF[,i], probs=1-PercentTail)
				}
			}
		}	
	return(Ints)
	}
