#   Find unique pairings
#   
#   This function finds all possible pairwise combinations for multiple runs.
#   
#   Used in several internal functions for analysis of results
#   
#   @param nRuns Numeric value representing haw many runs to compare
#   @return Returns a matrix of combinations
#   @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished

#' @name intrinsicModels
#' @rdname intrinsicModels
#' @export
pairings<-function (nRuns) {
#library(partitions)
#each output colum has a 1 in the rows corresponding to the runs you will combine for the ESS
#	library(partitions)
	possibilities<-blockparts(1:nRuns,nRuns,include.fewer=TRUE)
	possibilities<-possibilities[,which(apply(possibilities, 2, max)==1)]
	possibilities<-possibilities[,which(apply(possibilities, 2, sum)>1)]
	return(possibilities)
}
