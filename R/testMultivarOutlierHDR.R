#' Test If an Outlier Is Within a Highest Density Region of a Multivariate Cloud of Data Points
#'
#' This function takes a matrix, consisting of multiple observations
#' for a set of variables, with observations assumed to be independent
#' and identically distributed, calculates a highest density interval
#' for each of those variables (using
#' \link{\code{highestDensityInterval}}), and then tests if some
#' potential 'outlier' (a particular observation for that set of
#' variables), is within that highest density interval.
#' Typically, users may want to account for possible non-independence
#' of the variables, which could lead to inflated false-positive rates
#' with this test of coverage. To account for this, this function
#' allows for the option of applying principal component analysis to
#' rotate the variables, and then calculate the highest density
#' intervals from the orthogonal principal component scores.
#' 


#' @details
#' Quantiles, or highest density intervals are
#' essentially a univariate concept. They are calculated
#' independently on each variable, as if each variable was
#' independent. They are often used to try to piece together 
#' 'regions' of confidence or credibility in statistics.
#' Unfortunately, this means that that if variables are correlated, 
#' the box-like region they describe in the multivariate space
#' may not cover enough of the actual data scatter, while covering
#' lots of empty space.
#' 
#' There are several solutions. There are a number of approaches
#' for fitting ellipsoids to bivariate data. But what about
#' multivariate datasets, such as sets of parameter estimates from
#' the posterior of a Bayesian analysis, which may have an
#' arbitrary number of variables (and thus dimensions)?
#' 
#' A different solution, applied here when \code{pca = TRUE}, 
#' is to remove the non-independence among variables
#' by rotating the dataset with principal components analysis.
#' While this approach has its own flaws, such as assuming that the data
#' reflects a multivariate normal. However, this is absolutely an
#' improvement on not-rotating, particularly if being
#' done for the purpose of this functon, which is essentially to measure
#' coverage--to measure whether some particular set of values falls within a highest
#' density region (a multivariate highest density interval).
#' 

#' @inheritParams highestDensityInterval

#' @param pca If \code{TRUE} (the default), a principal components analysis is
#' applied to the provided input data, to reorient the data along orthogonal axes, 
#' with the purpose of potentially reducing covariation among the variables. If
#' your variables are known to have little covariation among them, and/or have
#' strange distributions that may violate the multivariate normal assumptions
#' of a principal components analysis, then it may be advisable to set {pca = FALSE}.

#' @param dataMatrix A matrix consisting of rows of data observations, 
#' with one or more variables as the columns, such that each multivariate
#' observation can be reasonably assumed to represent independent, 
#' identically-distributed random variables. For example, a matrix of
#' sampled parameter estimates from the post-burnin posterior of a Bayesian MCMC.

#' @param outlier A vector of consisting of a single observation of
#' one or more variables, with the same length as the number
#' of columns in /code{dataMatrix}. Not necessarily
#' a true 'outlier' \emph{per se} but
#' some data point of interest that you wish to test whether it is
#' inside some given probability density estimated from a sample of points.


#' @return
#' A logical, indicating whether the given data point (\code{outlier})
#' is within the multivariate data cloud at the specified threshold
#' probability density.

#' @seealso
#' This function is essentially a wrapper for applying \code{highestDensityInterval} to
#' multivariate data, along with \code{\link{princomp}} in package \code{stats}.

#' @author David W. Bapst

#' @examples
#' # First, you should understand the problems
#' # with dealing with correlated variables and looking at quantiles
#' 
#' # simulate two correlated variables
#' set.seed(444)
#' x <- rnorm(100, 1, 1)
#' y <- (x*1.5)+rnorm(100)
#' 
#' # find the highest density intervals for each variable
#' pIntX <- highestDensityInterval(x, alpha = 0.8)
#' pIntY <- highestDensityInterval(y, alpha = 0.8)
#' 
#' # These define a box-like region that poorly
#' # describes the actual distribution of
#' # the data in multivariate space.
#' 
#' # Let's see this ourselves...
#' xLims <- c(min(c(x, pIntX)), max(c(x, pIntX)))
#' yLims <- c(min(c(y, pIntY)), max(c(y, pIntY)))
#' plot(x, y, xlim = xLims, ylim = yLims)
#' rect(pIntX[1], pIntY[1], pIntX[2], pIntY[2])
#' 
#' # So, that doesn't look good.
#' # Let's imagine we wanted to test if some outlier
#' # was within that box:
#' 
#' outlier <- c(2, -1)
#' points(x = outlier[1], y = outlier[2], col = 2)
#' 
#' # we can use testMultivarOutlierHDR with pca = FALSE
#' # to do all of the above without visually checking
#' testMultivarOutlierHDR(dataMatrix = c(x, y), 
#' 	outlier = outlier, alpha = 0.8, pca = FALSE)
#' 	
#' # Should that really be considered to be within
#' # the 80% density region of this data cloud?
#' 
#' #####
#' 
#' # let's try it with PCA
#' 
#' pcaRes <- princomp(cbind(x, y))
#' PC <- pcaRes$scores
#' 
#' pIntPC1 <- highestDensityInterval(PC[, 1], alpha = 0.8)
#' pIntPC2 <- highestDensityInterval(PC[, 2], alpha = 0.8)
#' 
#' # plot it
#' xLims <- c(min(c(PC[, 1], pIntPC1)), max(c(PC[, 1], pIntPC1)))
#' yLims <- c(min(c(PC[, 2], pIntPC2)), max(c(PC[, 2], pIntPC2)))
#' plot(PC[, 1], PC[, 2], xlim = xLims, ylim = yLims)
#' rect(pIntPC1[1], pIntPC2[1], pIntPC1[2], pIntPC2[2])
#' 
#' # That looks a lot better, isnt' filled with lots of
#' # white space not supported by the observed data.
#' 
#' # And now let's apply testMultivarOutlierHDR, with pca = TRUE
#' testMultivarOutlierHDR(dataMatrix = c(x, y), 
#' 	outlier = outlier, alpha = 0.8, pca = TRUE)
#' 
#' #####################
#' 
#' # Example with four complex variables
#' 	# including correlated and multimodal variables
#' 
#' x <- rnorm(100, 1, 1)
#' y <- (x*0.8)+rnorm(100)
#' z <- sample(c(rnorm(50, 3, 2), 
#'   rnorm(50, 30, 3)))
#' a <- sample(c(rnorm(50, 3, 2), 
#'   rnorm(50, 10, 3)))+x^2
#' 
#' #plot(x, y)
#' #plot(x, z)
#' #plot(x, a)
#' data <- cbind(x, y, z, a)
#' 
#' # actual outlier, but maybe not obvious if PCA isn't applied
#' outlier <- c(2, 0.6, 10, 8)
#' # this one should appear to be an outlier (but only if PCA is applied)
#' testMultivarOutlierHDR(dataMatrix = data,
#'   outlier = outlier, alpha = 0.8)
#' testMultivarOutlierHDR(dataMatrix = data, 
#'   outlier = outlier, alpha = 0.8, pca = FALSE)
#' 
#' # this one should be within the 80% area
#' outlier <- c(1, 0, 30, 5)
#' testMultivarOutlierHDR(dataMatrix = data, 
#'   outlier = outlier, alpha = 0.8)
#' testMultivarOutlierHDR(dataMatrix = data, 
#'   outlier = outlier, alpha = 0.8, pca = FALSE)
#' 
#' # this one should be an obvious outlier no matter what
#' outlier <- c(3, -2, 20, 18)
#' # this one should be outside the 80% area
#' testMultivarOutlierHDR(dataMatrix = data, 
#'   outlier = outlier, alpha = 0.8)
#' testMultivarOutlierHDR(dataMatrix = data, 
#'   outlier = outlier, alpha = 0.8, pca = FALSE)
#' 





#' @name testMultivarOutlierHDR
#' @rdname testMultivarOutlierHDR
#' @export
testMultivarOutlierHDR <- function(dataMatrix, outlier, alpha, coda = FALSE, pca = TRUE, ...){
	if(!is.matrix(dataMatrix)){
		stop("dataMatrix is apparently not a matrix")
		}
	if(nrow(dataMatrix)<3){
		stop("Two or less observations in dataMatrix - insufficient data for generating a highest density region")
		}
	if(nrow(dataMatrix)>10){
		warning("Less than 10 observations are given as input - estimates of the highest density region may be very imprecise")
		}
	#
	dataAll <- rbind(data, outlier)
	# use PCA or not
	if(pca){
		pcaRes <- stats::princomp(dataAll)
		varia <- pcaRes$scores
	}else{
		# use raw
		varia <- dataAll
		}
	# separate out the outlier to be tested
	variaOut <- varia[nrow(varia), ]
	# now remove outlier from variable sample to get HPD from
	varia <- varia[-nrow(varia), ]
	#
	# get HPDs
	HDR <- lapply(1:ncol(varia), function(i)
		 highestDensityInterval(varia[, i], alpha = alpha, coda = coda, ...))
	#
	# test if outlier is within intervals listed for each variable
	withinHDR <- sapply(1:ncol(varia), function(x) 
		any(sapply(1:nrow(HPD[[x]]), function(y) 
			HPD[[x]][y, 1] <= variaOut[x] & HPD[[x]][y, 2] >= variaOut[x]
			))
		)
	res <- all(withinHPD)
	return(res)
	}
