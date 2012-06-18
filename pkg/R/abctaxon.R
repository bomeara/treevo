#create abctaxon class object
abctaxon <- function(id=0, name="taxon0", timeSinceSpeciation=0, states=NA, nextstates=NA) {
	taxon <- list(id=id, name=name, timeSinceSpeciation=timeSinceSpeciation, states=states, 
		nextstates=nextstates)
	class(taxon) <- "abctaxon"
	return(taxon)
}
