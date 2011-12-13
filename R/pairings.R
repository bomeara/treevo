pairings<-function (nRuns) {
library(partitions)
#each output colum has a 1 in the rows corresponding to the runs you will combine for the ESS
	library(partitions)
	possibilities<-blockparts(1:nRuns,nRuns,include.fewer=TRUE)
	possibilities<-possibilities[,which(apply(possibilities, 2, max)==1)]
	possibilities<-possibilities[,which(apply(possibilities, 2, sum)>1)]
	return(possibilities)
}