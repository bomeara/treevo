#' Highest Posterior Density
#' 
#' This function calculates the highest posterior density (HPD) for each freely
#' varying parameter, using \code{\link{HPDinterval}} in package \code{coda}.
#' 


# dealing with multimodal distributions

# quantiles don't work - based on ECDF, so they can only be a single interval
z<-sample(c(rnorm(50,1,2),rnorm(100,50,3)))

# how to get out the actual values?

# could do this using just density from stats
density(z)

alpha<-0.8

zDensOut <- density(z)
zDensity <- zDensOut$y/sum(zDensOut$y)
inHPD<-cumsum(-sort(-zDensity))<=alpha
# now reorder
inHPD<-inHPD[order(order(-zDensity))]
colDens<-rep(1,length(zDensity))
colDens[inHPD]<-2
plot(zDensOut$x,zDensity,col=colDens)

# let's make this a function

#' @param dataVector A vector of data, which can be reasonably assumed to
#' represent independent, identically-distributed random variables, such as
#' estimates of a parameter from a Bayesian MCMC.

#' @param prob The threshold used for defining the highest density frequency
#' cut-off.

#' @param ... Additional arguments passed to \code{density}. 
#' A user may want to mess with
#' this to adjust bandwidth, et cetera.

highestPostDensity<-function(dataVector, prob, ...){
	densOut <- density(dataVector, ...)
	densityScaled <- densOut$y/sum(densOut$y)
	#
	# count max number of ties
	maxTies<-max(table(densityScaled))
	# stop if more than half the dataset is tied
	if(maxTies>(length(dataVector)/2)){
		stop("Values of distribution are more than half tied with each other, may be flat")}
	#
	inHPD<-cumsum(-sort(-densityScaled))<=prob
	# now reorder
	inHPD<-inHPD[order(order(-densityScaled))]
	# get breaks
	startInt<-densOut$x[c(inHPD[1],
		sapply(2:length(inHPD),function(x) inHPD[x] & !inHPD[x-1])
		)]
	endInt<-densOut$x[c(
		sapply(1:(length(inHPD)-1),function(x) inHPD[x] & !inHPD[x+1])
		,inHPD[length(inHPD)]
		)]
	resMatrix<-cbind(startInt,endInt)
	return(resMatrix)
	}

highestPostDensity(z,prob=0.8)




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
#' highestPostDensity(results[[1]]$particleDataFrame, percent=0.95, returnData=FALSE)
#' 

#' @name highestPostDensity
#' @rdname highestPostDensity
#' @export


highestPostDensityCODA<-function(particleDataFrame, percent=0.95, returnData=FALSE){
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
		res<-summary
	}else{
		res<-as.data.frame(Ints)
		}
	return(res)
	}
