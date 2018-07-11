#' Summarize Posterior Distribution for a Free Parameter
#' 

#' This function summarizes the posterior distribution from the final generation for all free parameters,
#' outputting the mean, standard deviation and Highest Posterior Density (at a 0.8 alpha) for each parameter.
#' 

#' @param particleDataFrame A \code{particleDataFrame} object, as found among the output from \code{\link{doRun}} functions.



# @param returnData Option to return data that falls within HPD interval.

#' @return Returns a list, with each element of the list as a vector containing the weighted mean, standard deviation, and followed by the highest density intervals (e.g. the highest posterior density intervals). Because posterior estimates of parameter values may be multimodal, multiple sets of bounds may be reported for complex posterior distributions.

#* @return Returns a matrix with weighted mean, standard deviation, upper and lower credible
#* intervals for each free parameter.

#' @author Brian O'Meara and Barb Banbury

#' @seealso \code{\link{HPDinterval}} in package \code{coda}

# @references O'Meara and Banbury, unpublished


#' @examples
#' 
#' data(simRunExample)
#' 
#' highestDensityInterval(results[[1]]$particleDataFrame, alpha=0.95)
#' 



#' @name
#' @rdname
#' @export



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

summarizePosterior<-function(particleDataFrame, alpha=0.95){

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
				
				


				
				HPD<-highestDensityInterval(dataVector=data,alpha=alpha)
				
				
				Ints[i,3]<-codaHPD[1] #returns lower HPD		
				Ints[i,4]<-codaHPD[2] #returns upper HPD	
				
				
				if(returnData){
					summary[[i]]<-subpDF[,i][-c(which(subpDF[i]<Ints[i,3]), which(subpDF[,i]>Ints[i,4]))]
					}		
				}
			}
		}
	###########################
	#if(returnData){
	#	summary$summary<-Ints
	#	res<-summary
	#}else{
	#	res<-as.data.frame(Ints)
	#	}
	###########################
	res<-as.data.frame(Ints)
	return(res)
	}
	
