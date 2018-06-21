#' Highest Density Interval
#' 
#' This function calculates highest density intervals (HDIs) for a given univariate vector.
#' Typically, this is done to obtain highest posterior density (HPD) for each freely varying
#' parameter estimated in the posterior of a Bayesian MCMC analysis. If these intervals are
#' calculated for more than one variable, they are referred to instead as regions. 
#' 
#' By default,  is done by fitting a kernal density estimate (KDE) via R function \code{density}
#' with default bandwidth, rescaling the kernal to a density, sorting intervals along the KDE
#' by that density, and then summing these values from largest to smallest, until the desired
#' alpha is reached. This algorithm is quick, and accounts for potentially multimodal distributions,
#' or those with complex shapes, unlike unimodal intervals, such as quantiles, or the
#' \code{\link{HPDinterval}} in package \code{coda}.
#' 
#' Alternatively, a user can opt to use function \code{\link{HPDinterval}} from package \code{coda}
#' to calculate highest density intervals. These will work as long as the data has a single mode
#' - data with multiple modes may have overly large quantiles (to encompass those multiple modes),
#' resulting in overly wide estimates of coverage. 
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

#' @param alpha The threshold used for defining the highest density frequency
#' cut-off.

#* @param alpha Probability content of the highest probability density.

#' @param alpha Probability content of the highest posterior density (HPD).

#' @param alpha Probability density content of the 


#' @param coda Default is \code{FALSE}. If \code{TRUE}, unimodal highest density regions will
#' instead be calculated using \code{\link{HPDinterval}} from package \code{coda}.

#' @param ... Additional arguments passed to \code{density}. 
#' A user may want to mess with
#' this to adjust bandwidth, et cetera.

#' @examples
#' 
#' data(simRunExample)
#' 
#' highestDensityRegion(results[[1]]$particleDataFrame, alpha=0.95, returnData=FALSE)
#' 

highestDensityRegion(z,alpha=0.8)






#' @name highestDensityRegion
#' @rdname highestDensityRegion
#' @export

highestDensityRegion<-function(dataVector, alpha, coda, ...){
	# test that its a vector

	if(coda){
		codaHPD<-getHPDcoda(data,alpha)
		resMatrix<-matrix(codaHPD[1:2],1,2)
	}else{
		#
		densOut <- density(dataVector, ...)
		densityScaled <- densOut$y/sum(densOut$y)
		#
		# count max number of ties
		maxTies<-max(table(densityScaled))
		# stop if more than half the dataset is tied
		if(maxTies>(length(dataVector)/2)){
			stop("Values of distribution are more than half tied with each other, may be flat")}
		#
		inHPD<-cumsum(-sort(-densityScaled))<=alpha
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
		}
	# name the columns and rows
		# paste("LowerHPD_", alpha, sep=""), paste("UpperHPD_", alpha, sep="")
		
		
	return(resMatrix)
	}



# internal function for getting unimodal HDR using coda
getHPDcoda<-function(data,alpha){
	codaHPD<-coda::HPDinterval(coda::as.mcmc(subpDF[,i]), prob=alpha)
	return(codaHPD)
	}
