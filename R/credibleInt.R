# WE MAY WANT TO REMOVE THIS FUNCTION


#* Bayesian Credible Interval
#* 
#* This function calculates credible intervals for each free
#* parameter, essentially as quantiles. These will work as long as the data
#* has a single mode - multimodal data may have overly large
#* quantiles (to encompass multiple modes), in which case
#* calculating credible intervals as \code{\link{highestDensityRegion}}
#* may be a better idea.
#* 
#* 

#* @param particleDataFrame A \code{particleDataFrame} object, as found among the output from \code{\link{doRun}} functions.

#* @param alpha Probability content of the highest probability density.

#* @return Returns a matrix with weighted mean, standard deviation, upper and lower credible
#* intervals for each free parameter.

#* @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#* @examples
#* 
#* data(simRunExample)
#* credibleInt(results[[1]]$particleDataFrame)
#* 

# was named CredInt originally

#* @name credibleInt
#* @rdname credibleInt
#* @export
credibleInt<-function(particleDataFrame, alpha=0.95) {
	# ugh ugh
	#generation<-NULL #to appease R CMD CHECK
	# yes??? I think this is right, not sure
	generation<-particleDataFrame$generation
	#
	alphaTail<-(1-alpha)/2
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
				Ints[i,3]<-quantile(subpDF[,i], probs=0+alphaTail)
				Ints[i,4]<-quantile(subpDF[,i], probs=1-alphaTail)
				}
			}
		}	
	return(Ints)
	}
