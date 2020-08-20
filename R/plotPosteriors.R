#' Plot posterior distributions
#' 
#' For each free parameter in the posterior, this function creates a plot of the distribution of values estimated in the last generation.
#' This function can also be used to visually compare against true (generating) parameter values in a simulation.

#' @note
#' \code{realParam} and \code{realParamValues} should only be used with simulated data, where
#' the true values are known. 
#' 

#' @param particleDataFrame A \code{particleDataFrame} object returned by TreEvo ABC analyses, can be a
#' single data frame or a list of data frames.
#' If the \code{particleDataFrame} is a list of separate TreEvo runs, posteriors will
#' be layered over each other to check for repeatability.  The relative gray-scale of
#' posterior distributions in the plot depends on their total number of runs.

#' @param priorsList A \code{priorList} list object, returned by TreEvo \code{\link{doRun}}
#' functions, where each element is itself a list of length two, composed of the name for the
#' prior function (e.g. \code{'uniform'}, and a vector of the parameters for that prior distribution). Can be a
#' single such list, or a list of \code{priorList} lists from a series of analyses (i.e. a list of depth 3).
#' If the \code{priorList} is a list of lists (if it is found to be a list of depth 3, only the first list will be
#' evaluated, and all other lists will be ignored.
#' In other words, priors are expected to be identical for all runs considered by this function, such
#' the choice of using the list of priors from the first analysis in the list is acceptable for all comparisons.
	

#' @param realParam If \code{TRUE}, this function will plot line segments where real
#' parameter values are known. (Usually only true when simulated data is analyzed.)

#' @param realParamValues Values for real parameters, include a value for each
#' parameter (including fixed values). Otherwise should be \code{NA}.

#' @return Returns a plot for each free parameter.

#' @author Barb Banbury and Brian O'Meara

# @references O'Meara and Banbury, unpublished
# @keywords plotPosteriors

#' @seealso
#' \code{\link{plotPrior}} for a set of functions for visualizing and summarizing
#' prior and posterior distributions, including a visual comparison for single parameters.

#' @examples
#' data(simRunExample)
#' 
#' # make a list of particleDataFrames to plot multiple runs
#' resultsPDFlist <- lapply(resultsBMExample, function(x) x$particleDataFrame)
#' 
#' plotPosteriors(
#'    resultsPDFlist, 
#'    priorsList = resultsBMExample[[1]]$priorList, 
#'    realParam = TRUE, 
#'    realParamValues = c(ancStateExample, genRateExample)
#'    )
#' 


#' @name plotPosteriors
#' @rdname plotPosteriors
#' @export
plotPosteriors <- function(particleDataFrame, priorsList, realParam = FALSE, realParamValues = NA) {
    # particleDataFrame can be single or a list of particleDataFrames (1:n)
    #
    # priorsList can also be single list or a list of list (of lists) 
	#
	#	 Note that priorsList have to be the
    #     same to make comparison across runs, therefore if a list of priorsList is
    #     given, this function will extract only the first list
    #
    # priorsList should be $priorList from TreEvo output
    #
    # realParamValues should include a real value for every prior (fixed or not).
    #
    # fix this function so it works when "fixed" params are not last (ie the second param)

	if(class(particleDataFrame) == "data.frame"){
         # ugh ugh
        #generation <- NULL #to appease R CMD CHECK
        # yes??? I think this is right, not sure
        generation <- particleDataFrame$generation
    
        #make generation and other names by column so it works for partial and complete
        data1 <- subset(particleDataFrame[which(particleDataFrame$weight>0), ],
			generation == max(particleDataFrame$generation)
			)
        run <- rep(1, dim(data1)[1])
        all <- cbind(run, data1)
        }
        
    if(class(particleDataFrame) == "list"){
        all <- data.frame()
        for (list in 1:length(particleDataFrame)) {
            #make generation and other names by column so it works for partial and complete
            data1 <- subset(particleDataFrame[[list]][which(particleDataFrame[[list]]$weight>0), ], 
                generation == max(particleDataFrame[[list]]$generation))
            run <- rep(list, dim(data1)[1])
            particleDataFrame[[list]] <- cbind(run, data1)
            all <- rbind(all, particleDataFrame[[list]])        
        }
    }
	#
	# if priors list is a list of lists of lists, use first list only
	if(is.list(priorsList[[1]][[1]])){
		priorsList<-priorsList[[1]]
		}
	#
	areFreeParam <- sapply(priorsList,function(x) 
		x[[1]]  !=  "fixed")
    nFreeParams <- length(which(areFreeParam))
    #dev.new(width = 2.5 * nFreeParams, height = 3)
    nf <- layout(
		matrix(1 : nFreeParams, nrow = 1, byrow = TRUE),
		respect = TRUE)
    #layout.show(nf)
	#
    #alternatively, we can plotPrior above the posterior density to check if it is moving
    #nf <- layout(matrix(1:(2*paramsToPlot), nrow = 2, byrow = FALSE), respect = TRUE)
    ##layout.show(nf)
	#
    v <- vector("list", max(all$run))
    nParticles <- dim(subset(all[which(all$weight>0), ], run == max(all$run)))[1]
	#
    for (param in 1:length(priorsList)) {
        #message(param)
        #v <- vector("list", max(all$run))
        which.param <- param+7
        r <- c()
        q <- c()
		typePrior <- priorsList[[param]][[1]]
		paramPriorParams <- priorsList[[param]][[2]]
		#
        if (areFreeParam[param]) {
			#
            if (typePrior  ==  "uniform" || typePrior  ==  "exponential") {
                for (dens in 1:max(all$run)) {
                    v[[dens]] <- density(subset(all[which(all$run == dens), ], )[, which.param], 
                        weights = nParticles*subset(all[which(all$weight>0), ], run == dens)[, 7]/sum(
                            nParticles*subset(all[which(all$weight>0), ], run == dens)[, 7]), 
                        from = min(subset(all[which(all$run == dens), ], )[, which.param]), 
                        to = max(subset(all[which(all$run == dens), ], )[, which.param]))
                    v[[dens]]$x <- c(min(v[[dens]]$x), v[[dens]]$x, max(v[[dens]]$x))
                    v[[dens]]$y <- c(0, v[[dens]]$y, 0)
                    q <- c(q, v[[dens]]$x)
                    r <- c(r, v[[dens]]$y)
                }
            }else{
                if (typePrior  ==  "normal" ||
                        typePrior  ==  "lognormal" ||
                        typePrior  ==  "gamma"){
                    ####################################
                    for (dens in 1:max(all$run)) {
                        v[[dens]] <- density(subset(all[which(all$run == dens), ], )[, which.param], 
                            weights = nParticles*subset(all[which(all$weight>0), ], run == dens)[, 7]/sum(
                                nParticles*subset(all[which(all$weight>0), ], run == dens)[, 7]))
                        q <- c(q, v[[dens]]$x)
                        r <- c(r, v[[dens]]$y)
                        }
                    }
                }

            if (typePrior == "uniform") {
                #message("made it to uniform priorFn")
                min = as.numeric(min(paramPriorParams[1], paramPriorParams[2]))
                max = as.numeric(max(paramPriorParams[1], paramPriorParams[2]))
                x <- runif(1000, min, max)
                poly <- curve(dunif(x), 
                    xlim = c(min(min-(.3*(max-min)),
                        min(q)),max(max+(.3*(max-min)), max(q))), 
                    ylim = c(0, max(1.5, max(r))), 
                    type = "n", 
                    xlab = names(all)[which.param], 
                    ylab = "Density", bty = "n")
                rect(min, 0, max, 1, border = "blue", lwd = 1.5)
                if (realParam) {
                    segments(realParamValues[param], 0, realParamValues[param],
                        max(1.5, max(r)), col = "red", lwd = 1.5)
                }
            }
            else if (typePrior == "normal") {
                #message("made it to normal priorFn")
                mean = as.numeric(paramPriorParams[1])
                stdev = as.numeric(paramPriorParams[2])
                x <- rnorm(1000, mean, stdev)
                w <- density(x)
                poly <- curve(dnorm(x, mean, stdev), 
                    from = min(x), to = max(x), 
                    xlim = c(min(min(w$x), min(q)), max(max(w$x), max(q))),
                    ylim = c(0, max(max(w$y), max(r))), 
                    xlab = names(all)[which.param], 
                    ylab = "Density", col = "blue", 
                    lwd = 1.5, bty = "n")
                if (realParam) {
                    segments(realParamValues[param], 0, realParamValues[param],
                        max(max(w$y), max(r)), col = "red", lwd = 1.5)
                }
            }
            else if (typePrior == "lognormal") {
                mean = as.numeric(paramPriorParams[1])
                stdev = as.numeric(paramPriorParams[2])
                x <- rlnorm(1000, mean, stdev)
                w <- density(x)
                poly <- curve(dlnorm(x, mean, stdev), 
                    from = 0, to = qlnorm(0.99, mean, stdev), 
                    xlim = c(min(min(w$x), min(q)), max(max(w$x), max(q))), 
                    ylim = c(0, max(max(w$y), max(r))), 
                    xlab = names(all)[which.param], 
                    ylab = "Density", col = "blue", 
                    lwd = 1.5, bty = "n")
                if (realParam) {
                    segments(realParamValues[param], 0, realParamValues[param],
                        max(max(w$y), max(r)), col = "red", lwd = 1.5)
                }
            }
            else if (typePrior == "gamma") {
                shape = as.numeric(paramPriorParams[1])
                scale = as.numeric(paramPriorParams[2])
                x <- rgamma(1000, shape, scale)
                w <- density(x)
                poly <- curve(dgamma(x, shape, scale), 
                    from = 0, to = qgamma(0.99, shape, scale), 
                    xlim = c(min(min(w$x), min(q)), max(max(w$x), max(q))), 
                    ylim = c(0, max(max(w$y), max(r))), 
                    xlab = names(all)[which.param], ylab = "Density", col = "blue", lwd = 1.5, bty = "n")
                if (realParam) {
                    segments(realParamValues[param], 0, realParamValues[param], 
                        max(max(w$y), max(r)), col = "red", lwd = 1.5)
                }
            }
            else if (typePrior == "exponential") {
                #message("made it to exponential priorFn")
                rate = as.numeric(paramPriorParams[1])
                x <- rexp(1000, rate)
                w <- density(x)
                poly <- curve(dexp(x, rate), from = 0, to = qexp(0.99, rate), 
                    xlim = c(min(min(w$x), min(q)), max(max(w$x), max(q))),
                    ylim = c(0, max(max(w$y), max(r))), 
                    xlab = names(all)[which.param],
                    ylab = "Density", 
                    col = "blue", lwd = 1.5, bty = "n"
                    ) #, sub = names(data)[which.param]
                if (realParam) {
                    segments(realParamValues[param], 0, realParamValues[param],
                        max(max(w$y), max(r)), col = "red", lwd = 1.5)
                }
            }
            
        ## plotting    
            #plot(q, r, type = "n", ylab = "density", xlim = c(min(q), max(q)),
            #ylim = c(min(r), max(r)), sub = names(all)[which.param])        
            #
            
            for (lines in 1:length(v)) {
                polygon(v[[lines]]$x, v[[lines]]$y,
                    border = NA,
                    col = rgb(0, 0, 0, (.7/max(all$run))))
            }
        } #if (typePrior  !=  "fixed")
    } #for
    layout(1)
    }
