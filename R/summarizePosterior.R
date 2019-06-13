#' Summarize Posterior Distribution for a Free Parameter
#' 
#' This function summarizes the distribution of parameter
#' estimates from the posterior of an ABC analysis using the
#' \code{doRun} functions in \code{TreEvo}, for all
#' freely varying parameters. Only the final
#' generation is considered. This summary includes
#' the mean, standard deviation and Highest Posterior
#' Density (at a 0.8 alpha) for each parameter.
#' 

#' @param particleDataFrame A \code{particleDataFrame} object, as
#' found among the output from \code{\link{doRun}} functions.

#' @inheritParams highestDensityInterval

#' @return 
#' Returns a list, wherein each element of the list is secondary list containing
#' the weighted mean, standard deviation, and a matrix giving the highest density
#' intervals (e.g. the highest posterior density intervals). Because posterior
#' estimates of parameter values may be multimodal, multiple sets of bounds
#' may be reported for complex posterior distributions, which each constitute one 
#' row of the output matrix. See \code{\link{highestDensityInterval}} for details.

#' @author David W Bapst

#' @seealso
#' This function is essentially a wrapper for independently applying 
#' a few summary statistics and applying
#' \code{\link{highestDensityInterval}} to multiple parameter estimates, taken from
#' the last generation of an ABC analysis in \code{TreEvo}. As each parameter
#' is handled independently, the returned HPD intervals may not properly account
#' for covariation among parameter estimates from the posterior. If testing
#' whether a given observation is within a given density of the posterior or
#' not, please look at function \code{\link{testMultivarOutlierHDR}}.

#' @examples
#' # example with output from doRun_prc
#' data(simRunExample)
#'
#' # alpha = 0.95
#' summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.95)
#' 
#' # you might be tempted to use alphas like 95%, but with bayesian statistics
#' # we often don't sample the distribution well enough to know
#' # its shape to exceeding detail. alpha = 0.8 may be more reasonable.
#' summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.8)
#' 
#' # or even better, for coverage purposes, maybe 0.5
#' summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.5)
#' 


#' @name summarizePosterior
#' @rdname summarizePosterior
#' @export
summarizePosterior <- function(particleDataFrame, 
        alpha = 0.8, 
        coda = FALSE, 
        verboseMultimodal = TRUE, 
		stopIfFlat = TRUE,
        ...){
    ########################################################
    #    
    generation <- particleDataFrame$generation
    # get the max generation
    maxGen <- max(particleDataFrame$generation)
    # subset on particles from the last generation
    subpDF <- particleDataFrame[
		particleDataFrame$generation == maxGen & particleDataFrame$weight>0,
		]    
    #
    # now only get the parameter estimates
    subpDF <- as.data.frame(subpDF[, 7:dim(particleDataFrame)[2]])
    
    #
    # old code - for removing parameters with bad values??
    # maybe applicable for fixed params?
    #
    for(i in 1:dim(subpDF)[2]){
        if(length(subpDF[, i])>1){
            if(sd(subpDF[, i], na.rm = TRUE)  ==  0) {
                subpDF <- subpDF[, -i]
                }
            }
        }
    #
    res<-apply(
		subpDF,2,
		getSummary,
        alpha = alpha, 
		coda = coda, 
        verboseMultimodal=verboseMultimodal, 
		stopIfFlat = stopIfFlat,
		...)
    # name the elements of the list
    names(res) <- colnames(subpDF)
    return(res)
    }
                
getSummary <- function(
		param, 
		alpha, 
		coda, 
        verboseMultimodal, 
		stopIfFlat,
		...
		){
    ######################################
    #
    HPD <- highestDensityInterval(
				dataVector = param,
                alpha = alpha, 
				coda = coda, 
                verboseMultimodal=verboseMultimodal, 
				stopIfFlat = stopIfFlat,
				...)
    output<-list(
        mean = mean(param),
        sd = sd(param),
        HPD = HPD
        )
    return(output)
    }
    
    
##################################################################################################    
# OLD
#    
#        generation <- particleDataFrame$generation
#    #
#    alphaTail <- (1-alpha)/2
#    Ints <- matrix(nrow = (dim(particleDataFrame)[2]-6), ncol = 4)
#    colnames(Ints) <- c("mean", "sd", "LowerCI", "UpperCI")
#    rownames(Ints) <- names(particleDataFrame[7: dim(particleDataFrame)[2]])
#    subpDF <- subset(particleDataFrame[which(particleDataFrame$weight>0), ], 
#        generation == max(particleDataFrame$generation), drop = FALSE)[7:dim(particleDataFrame)[2]]
#    for(i in 1:dim(subpDF)[2]){
#        if(length(subpDF[, i])>1){
#            if(sd(subpDF[, i])  !=  0) {
#                Ints[i, 1] <- mean(subpDF[, i]) #not weighted
#                Ints[i, 2] <- sd(subpDF[, i])
#                Ints[i, 3] <- quantile(subpDF[, i], probs = 0+alphaTail)
#                Ints[i, 4] <- quantile(subpDF[, i], probs = 1-alphaTail)
#                }
#            }
#        }    
#    

# @param returnData Option to return data that falls within HPD interval.

