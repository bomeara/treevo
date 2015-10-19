#create abctaxon class object


#' abctaxon
#' 
#' Internal TreEvo function to create a list of objects
#' 
#' 
#' @param id Simulation ID
#' @param name Taxon of interest
#' @param timeSinceSpeciation Time since speciation
#' @param states States of taxon
#' @param nextstates States of taxon after mutation
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
abctaxon <- function(id=0, name="taxon0", timeSinceSpeciation=0, states=NA, nextstates=NA) {
	taxon <- list(id=id, name=name, timeSinceSpeciation=timeSinceSpeciation, states=states, 
		nextstates=nextstates)
	class(taxon) <- "abctaxon"
	return(taxon)
}
