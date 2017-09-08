# Simulation Splits
# 
# This function retrieves branch lengths and ancestor-decendant relationships
# from a tree
# 
# This function is used by other TreEvo functions for internal calculations.
# 

# @param phy A phylogenetic tree, in package \code{ape}'s \code{phylo} format.

# @return A data frame of branching times, ancestor and descendant vectors

# @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

# @name intrinsicModels
# @rdname intrinsicModels
# @export
getSimulationSplits<-function(phy) {
	phy$node.label<-NULL
	branchingTimes<-sort(branching.times(phy), decreasing=TRUE)
	branchingTimesNames<-names(branchingTimes)
	ancestorVector<-c()
	descendant1Vector<-c()
	descendant2Vector<-c()
	for (i in 1:length(branchingTimes)) {
		relationshipVector<-phy$edge[phy$edge[, 1]==branchingTimesNames[i]]
		ancestorVector<-c(ancestorVector, relationshipVector[2])
		descendant1Vector<-c(descendant1Vector, relationshipVector[3])
		descendant2Vector<-c(descendant2Vector, relationshipVector[4])
	}
	
	return(data.frame(branchingTimes, ancestorVector, descendant1Vector, descendant2Vector))
}
