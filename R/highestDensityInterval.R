#' Highest Density Interval
#' 
#' This function calculates highest density intervals (HDIs) for a given univariate vector.
#' Typically, this is done to obtain highest posterior density (HPD) for each freely varying
#' parameter estimated in the posterior of a Bayesian MCMC analysis. If these intervals are
#' calculated for more than one variable, they are referred to instead as regions. 
#' 
#' By default, HDI calculation is preformed by fitting
#' a kernal density estimate (KDE) via R function \code{density}
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

#' @param dataVector A vector of data, which can be reasonably assumed to
#' represent independent, identically-distributed random variables, such as
#' estimates of a parameter from a Bayesian MCMC.

#' @param alpha The threshold used for defining the highest density frequency
#' cut-off. If the highest density interval is applied to a Bayesian MCMC
#' posterior sample, then the interval is effectively calculated for this
#' value as a posterior probability density.

#' @param coda Default is \code{FALSE}. If \code{TRUE}, unimodal highest density regions will
#' instead be calculated using \code{\link{HPDinterval}} from package \code{coda}, which is
#' similar to the function \code{quantile} in that it calculates only a single interval.

#' @param ... Additional arguments passed to \code{\link{density}}. 
#' A user may want to mess with this to adjust bandwidth, et cetera.

#' @param verboseMultimodal If \code{TRUE}, the function will print a message
#' indicating when the inferred highest density interval is discontinuous, and
#' thus likely reflects that the supplied data is multimodal.

#' @references
#' Hyndman, R. J. 1996. Computing and Graphing Highest
#' Density Regions. \emph{The American Statistician} 50(2):120-126.

#' @seealso \code{\link{HPDinterval}} in package \code{coda} for an
#' alternative used by older versions of \code{TreEvo}
#' which cannot properly handle multimodal distributions.

#' @examples
#' set.seed(444)
#' 
#' # let's imagine we have some variable with
#' 	# an extreme bimodal distribution
#' z <- sample(c(rnorm(50, 1, 2), rnorm(100, 50, 3)))
#' hist(z)
#' 
#' # now let's say we want to know the what sort of values
#' # are reasonably consistent with this distribution
#' 
#' # for example, let's say we wanted the ranges within
#' # which 80% of our data occurs
#' 
#' # one way to do this would be a quantile
#' 	# two tailed 80% quantiles
#' quantile(z, probs = c(0.1, 0.9))
#' 
#' # that seems overly broad - there's essentially no density
#' # in the central valley - but we want to exclude values found there!
#' # A value of 30 doesn't match this data sample, right??
#' 
#' # the problem is that all quantile methods are essentially based on
#' # the empirical cumulative distribution function - which is monotonic
#' # (as any cumulutative function should be), and thus
#' # quantiles can only be a single interval
#' 
#' # A different approach is to use density from stats
#' density(z)
#' # we could then take the density estimates for
#' # particular intervals, rank-order them, and
#' # then cumulative sample until we reach
#' # our desired probability density (alpha)
#' 
#' # let's try that
#' alpha <- 0.8
#' zDensOut <- density(z)
#' zDensity <- zDensOut$y/sum(zDensOut$y)
#' inHPD <- cumsum(-sort(-zDensity)) <= alpha
#' # now reorder
#' inHPD <- inHPD[order(order(-zDensity))]
#' colDens <- rep(1, length(zDensity))
#' colDens[inHPD] <- 2
#' # and let's plot it, with colors
#' plot(zDensOut$x, zDensity, col = colDens)
#' 
#' # and we can do all that (except the plotting)
#' 	# with highestDensityInterval
#' highestDensityInterval(z, alpha = 0.8)
#' 
#' #############################
#' # example with output from doRun_prc
#' 
#' data(simRunExample)
#' 
#' # we'll use summarizePosterior, which just automates picking
#'   # the last generation, and freely-varying parameters for HDRs 
#'
#' # alpha = 0.95
#' summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.95)
#' 
#' # you might be tempted to use alphas like 95%, but with bayesian statistics
#' # we often don't sample the distribution well enough to know
#' # its shape to exceeding detail. alpha = 0.8 may be more reasonable.
#' summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.8)
#' 
#' 





#' @name highestDensityInterval
#' @rdname highestDensityInterval
#' @export
highestDensityInterval <- function(dataVector, alpha,
		coda = FALSE, verboseMultimodal = TRUE,...){
	#
	# test that its a vector
	dataVector <- as.numeric(dataVector)
	if(!is.vector(dataVector)){
		stop("Unable to coerce dataVector to being a numeric vector")
		}
	# test alpha is of length 1
	if(length(alpha)>1){
		stop("alpha should be of length 1")
		}
	if(alpha >= 1 | alpha <= 0){
		stop("alpha should be numeric and be between zero and one")
		}
	#
	if(coda){
		codaHPD <- getHPDcoda(dataVector, alpha)
		resMatrix <- matrix(codaHPD[1:2], 1, 2)
	}else{
		#
		densOut <- density(dataVector, ...)
		densityScaled <- densOut$y/sum(densOut$y)
		#
		# count max number of ties
		maxTies <- max(table(densityScaled))
		# stop if more than half the dataset is tied
		if(maxTies>(length(dataVector)/2)){
			stop("Values of distribution are more than half tied with each other, may be flat")}
		#
		inHPD <- cumsum(-sort(-densityScaled)) <= alpha
		# now reorder
		inHPD <- inHPD[order(order(-densityScaled))]
		# get breaks
		startInt <- densOut$x[c(inHPD[1], 
			sapply(2:length(inHPD), function(x) inHPD[x] & !inHPD[x-1])
			)]
		endInt <- densOut$x[c(
			sapply(1:(length(inHPD)-1), function(x) inHPD[x] & !inHPD[x+1])
			, inHPD[length(inHPD)]
			)]
		resMatrix <- cbind(startInt, endInt)
		}
	# name the columns and rows
		# paste("LowerHPD_", alpha, sep = ""), paste("UpperHPD_", alpha, sep = "")
	colLabels <- c(paste0(c("LowerBound_alpha=", "UpperBound__alpha="), alpha))
	colnames(resMatrix) <- colLabels
	#
	if(nrow(resMatrix)>1 & verboseMultimodal){
		message("Inferred highest density intervals are discontinuous, suggesting data may have multiple modes")
		}
	#
	return(resMatrix)
	}



# internal function for getting unimodal HDR using coda
getHPDcoda <- function(data, alpha){
	codaHPD <- coda::HPDinterval(coda::as.mcmc(data), prob = alpha)
	return(codaHPD)
	}
