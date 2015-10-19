#' Naming function for summary stats
#' 
#' This function creates a vector of summary stat names
#' 
#' 
#' @param phy Tree (Phylogenetic tree in phylo format)
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
sumStatNames<-function(phy) {
  ntips=Ntip(phy)
  return(c("BML", "BMbeta", "BMaic", "lambdaL", "lambdabeta", "lambdalambda", "lambdaaic", "deltaL", "deltabeta", "deltadelta", "deltaaic", "ouL", "oubeta", "oualpha", "ouaic", "whiteL", "whiteaic", "rawmean", "rawmax", "rawmin", "rawvar", "rawmedian", phy$tip.label, paste("pic", sequence(ntips-1), sep=""), paste("ancState", sequence(ntips-1), sep=""), paste("ancCIrange", sequence(ntips-1), sep="")))
}
