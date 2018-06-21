#' Summarize Posterior Distribution for a Free Parameter
#' 

#' This function summarizes the posterior distribution from the final generation for all free parameters, outputting the mean, standard deviation and Highest Posterior Density (at a 0.8 alpha) for each parameter.
#' 

#' @param particleDataFrame \code{particleDataFrame} output from \code{doRun}

#' @param alpha Probability content of the highest posterior density (HPD).

#' @param returnData Option to return data that falls within HPD interval.

#' @return Returns a matrix with weighted mean, standard deviation, upper and lower HPD for
#' each free parameter.

#' @author Brian O'Meara and Barb Banbury

#' @seealso \code{\link{HPDinterval}} in package \code{coda}

# @references O'Meara and Banbury, unpublished


#' @examples
#' 
#' data(simRunExample)
#' 
#' highestDensityRegion(results[[1]]$particleDataFrame, alpha=0.95, returnData=FALSE)
#' 



#' @name
#' @rdname
#' @export


summarizePosterior<-function(particleDataFrame, alpha=0.95, returnData=FALSE){
  	# ugh ugh
	#generation<-NULL #to appease R CMD CHECK
	# yes??? I think this is right, not sure
	generation<-particleDataFrame$generation
	maxGen<-max(particleDataFrame$generation)
#	library(coda, quietly=TRUE)
	summary<-vector("list")
	# find particles from the last generation
	subpDF<-particleDataFrame[particleDataFrame$generation==maxGen
		& particleDataFrame$weight>0,]
	#
	# now only get the parameter estimates
	subpDF<-as.data.frame(subpDF[,7:dim(particleDataFrame)[2]])
	#parNames<-colnames(subpDF)
	#
	for(i in 1:dim(subpDF)[2]){
		if(length(subpDF[,i])>1){
			if(sd(subpDF[,i], na.rm=TRUE) == 0) {
				subpDF<-subpDF[,-i]
				}
			}
		}
	Ints<-matrix(nrow=(dim(subpDF)[2]), ncol=4)
	colnames(Ints)<-c("mean", "sd")	
	rownames(Ints)<-names(subpDF)
	for(i in 1:dim(subpDF)[2]){
		if(length(subpDF[,i])>1){
			if(sd(subpDF[,i]) != 0) {
				Ints[i,1]<-mean(subpDF[,i]) #not weighted
				Ints[i,2]<-sd(subpDF[,i])
				
				


				
				HPD<-highestDensityRegion(dataVector=data,alpha=alpha)
				
				
				Ints[i,3]<-codaHPD[1] #returns lower HPD		
				Ints[i,4]<-codaHPD[2] #returns upper HPD	
				
				
				if(returnData){
					summary[[i]]<-subpDF[,i][-c(which(subpDF[i]<Ints[i,3]), which(subpDF[,i]>Ints[i,4]))]
					}		
				}
			}
		}
	if(returnData){
		summary$summary<-Ints
		res<-summary
	}else{
		res<-as.data.frame(Ints)
		}
	return(res)
	}