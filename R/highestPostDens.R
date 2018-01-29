#' Highest Posterior Density
#' 
#' This function calculates the highest posterior density (HPD) for each freely varying parameter.
#' 
#' 

#' @param particleDataFrame \code{particleDataFrame} output from \code{doRun}

#' @param percent Probability content of the highest posterior density (HPD).

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
#' highestPostDens(results$particleDataFrame, percent=0.95, returnData=FALSE)
#' 

#' @name highestPostDens
#' @rdname highestPostDens
#' @export
highestPostDens<-function(particleDataFrame, percent=0.95, returnData=FALSE){
  	# ugh ugh
	#generation<-NULL #to appease R CMD CHECK
	# yes??? I think this is right, not sure
	generation<-particleDataFrame$generation
	
#	library(coda, quietly=TRUE)
	summary<-vector("list")
	subpDF<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame$generation))[7:dim(particleDataFrame)[2]])
	for(i in 1:dim(subpDF)[2]){
		if(length(subpDF[,i])>1){
			if(sd(subpDF[,i], na.rm=TRUE) == 0) {
				subpDF<-subpDF[,-i]
				}
			}
		}
	Ints<-matrix(nrow=(dim(subpDF)[2]), ncol=4)
	colnames(Ints)<-c("mean", "sd", paste("LowerHPD_", percent, sep=""), paste("UpperHPD_", percent, sep=""))	
	rownames(Ints)<-names(subpDF)
	for(i in 1:dim(subpDF)[2]){
		if(length(subpDF[,i])>1){
			if(sd(subpDF[,i]) != 0) {
				Ints[i,1]<-mean(subpDF[,i]) #not weighted
				Ints[i,2]<-sd(subpDF[,i])
				Ints[i,3]<-coda::HPDinterval(coda::as.mcmc(subpDF[,i]), prob=percent)[1] #returns lower HPD		
				Ints[i,4]<-coda::HPDinterval(coda::as.mcmc(subpDF[,i]), prob=percent)[2] #returns upper HPD	
				if(returnData){
					summary[[i]]<-subpDF[,i][-c(which(subpDF[i]<Ints[i,3]), which(subpDF[,i]>Ints[i,4]))]
					}		
				}
			}
		}
	if(returnData){
		summary$summary<-Ints
		return(summary)
	}
	else{return(as.data.frame(Ints))}
}
