#create abctaxonS3 class object
abctaxonS3 <- function(id=0, name="taxon0", timeSinceSpeciation=0, states=NA, nextstates=NA) {
	taxon <- list(id=id, name=name, timeSinceSpeciation=timeSinceSpeciation, state=states, 
		nextstates=nextstates)
	class(taxon) <- "abctaxonS3"
	return(taxon)
}
