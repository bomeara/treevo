#' Pairwise Effective Sample Size
#' 
#' This function calculates Effective Sample Size (ESS) on results.  
#' The ESS algorithm performs best when results are
#' from multiple runs.
#' 
#' 

#' @param inputData The input dataset can be a single data frame or a
#' list of data frames.

#' @return Returns a matrix with ESS values of all pairwise runs.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#' 
#' data(simRunExample)
#' 
#' # this will give a warning
#' pairwiseESS(results[[1]]$particleDataFrame)
#' 
#' # ESS should be calculated over multiple runs
#' pairwiseESS(results)
#'
#' # you can also manually assemble a list of particleDataFrame tables
#'     # and use this as the input
#' inputList <- list(results[[2]]$particleDataFrame, results[[1]]$particleDataFrame)
#' pairwiseESS(inputList)

#' @name pairwiseESS
#' @rdname pairwiseESS
#' @export
pairwiseESS <- function(inputData) {
    #Combine doRun$particleDataFrame results and test effective sample size    
    # inputData can be single or a list of particleDataFrames (1:n)
    if ( inherits(inputData, "multiRun_doRun_prc") ) {
        particleInput <- lapply(inputData, function(x) x$particleDataFrame)
		if(length(particleInput) < 2){
				message(
					"A single run of doRun_prc found, extracting particleDataFrame object"
					)
				particleInput <- inputData$particleDataFrame
			}else{
				message(
					"Multiple runs of doRun_prc found, extracting particleDataFrame objects"
					)
				}		
    }else{
        particleInput <- inputData
        }
    #
    ESSmatrix <- NA
    #
    if(class(particleInput) == "data.frame"){
        warning(
			"ESS on a single run should be high; consider combining several runs"
			)
		#
        data1 <- subset(particleInput, 
			particleInput$generation == max(particleInput[, 1])
			)
        data1 <- data1[which(data1$weight>0), ]
        # get ESS
        ESS <- coda::effectiveSize(
			data1[, 7:dim(data1)[2]]
			)
        ESSmatrix <- matrix(c(as.numeric(ESS)))
        rownames(ESSmatrix) <- names(ESS)
        colnames(ESSmatrix) <- "ESS"
        }
    #
    if(class(particleInput) == "list"){
		if(length(particleInput) < 2){
			stop(
				"The number of runs provided is less than two - no pairwise comparisons possible."
				)
			}
        pairs <- as.matrix(pairings(length(particleInput)))
        ESS <- vector("list")
        ESSmatrix <- matrix(
			nrow = dim(particleInput[[1]])[2]-6, 
			ncol = dim(pairs)[2]
			)
        rownames(ESSmatrix) <- names(particleInput[[1]])[
			7:length(names(particleInput[[1]]))
			]
        colnames(ESSmatrix) <- seq(1:dim(ESSmatrix)[2])
        data1 <- vector("list")
        for (list in 1:length(particleInput)) {
            data1[[list]] <- subset(
				particleInput[[list]][which(particleInput[[list]][, 6]>0), ], 
                particleInput[[list]]$generation == max(
					particleInput[[list]]$generation
					)
				)
            }
        for (combination in 1:dim(pairs)[2]){
            runNames <- c()
            data2 <- c()
            runNames <- c()
            for (run in 1:dim(pairs)[1]){
                if (pairs[run, combination] == 1){
                    data2 <- rbind(data2, data1[[run]])
                    runNames <- append(runNames, run)
                    }
                }
            # get ESS
            ESS[[combination]] <- coda::effectiveSize(
				data2[, 7:dim(data2)[2]]
				)
			# name by combining runNames
            names(ESS)[[combination]] <- paste0("Runs_", 
				paste0(runNames, collapse = "_")
				)
			#
			# define ESS matrix
            ESSmatrix[, combination] <- c(as.numeric(ESS[[combination]]))
			# give column names to ESSmatrix
            colnames(ESSmatrix)[combination] <- paste0("Runs_", 
				paste0(runNames, collapse = "_")
				)
            }
        }
    return(ESSmatrix)
    }

    
    
    
    
#   Find unique pairings
#
#   This function finds all possible pairwise combinations for multiple runs.
#
#   Used in several internal functions for analysis of results
#
#   @param nRuns Numeric value representing haw many runs to compare
#   @return Returns a matrix of combinations
#   @author Brian O'Meara and Barb Banbury
#	@references O'Meara and Banbury, unpublished


pairings <- function (nRuns) {
	#library(partitions)
	# check
	if(nRuns < 2){
		stop("Inadequate sample to draw pairwise comparisons from.")
		}
	#
	#each output column has a 1 in the rows corresponding
		#to the runs you will combine for the ESS
	#
    possibilities <- blockparts(1:nRuns, nRuns, include.fewer = TRUE)
    possibilities <- possibilities[, which(apply(possibilities, 2, max) == 1)]
    possibilities <- possibilities[, which(apply(possibilities, 2, sum)>1)]
    return(possibilities)
    }
