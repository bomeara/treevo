# Get data.frame of taxa
#
# This function retrieves branch lengths and ancestor-decendant relationships
# from a tree.
#
# This function is used by other TreEvo functions for internal calculations.
#

# @inheritParams doSimulation

# @return A list of abctaxon objects

# @author Brian O'Meara

# @references O'Meara and Banbury, unpublished

# @name getTaxonDFWithPossibleExtinction
# @rdname getTaxonDFWithPossibleExtinction
# @export
getTaxonDFWithPossibleExtinction<-function(phy) {
	heights <- data.frame(cbind(phy$edge, phy$edge.length, phytools::nodeHeights(phy)))
	colnames(heights) <- c("rootward.id", "tipward.id", "edge.length", "rootward.height", "tipward.height")
	heights$name <- paste("internal", heights$tipward.id, sep=".")
	for (taxon.index in sequence(ape::Ntip(phy))) {
		heights$name[which(heights$tipward.id==taxon.index)] <- phy$tip.label[taxon.index]
	}
	result <- data.frame(t(sapply(split(heights, seq(nrow(heights))), createAbcTaxonFromHeightsRow)), stringsAsFactors=FALSE)
	for (col.id in sequence(dim(result)[2])) {
		result[,col.id] <- unlist(result[,col.id])
	}
	return(result)
}


# Create abctaxon object
#
# This creates an abctaxon object from a row in the heights df
# created within getAbcTaxonListWithPossibleExtinction
#
# @param x A row of the heights df
# @return An abctaxon object
# @author Brian O'Meara
createAbcTaxonFromHeightsRow <- function(x) {
	return(abctaxon(id=x$tipward.id, name=x$name, timeSinceSpeciation=0, states=NA, nextstates=NA, 
	startTime=x$rootward.height, endTime=x$tipward.height, terminal=!grepl("internal.", x$name), ancestorId=x$rootward.id))
}
