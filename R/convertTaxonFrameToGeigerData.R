# Converting Simulated Data to Geiger Formatted Data
# 
# This function takes simulated data and converts it to geigers formatting.
# 
# 

# @param taxonframe Data frame with species values

# @InheritParams doSimulation

# @return Returns a data frame in geiger format

# @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished


# @name convertTaxonFrameToGeigerData
# @rdname convertTaxonFrameToGeigerData
# @export
# convertTaxonFrameToGeigerData<-function(taxonframe, phy) {
# 	ntax<-dim(taxonframe)[1]
# 	newmatrix<-matrix(data=taxonframe[, 4:(dim(taxonframe)[2])], nrow= ntax)
# 	newrownames<-c(rep(0, ntax))
# 	for (i in 1:ntax) {
# 		newrownames[i]<-phy$tip.label[(taxonframe$taxonid[i])]
# 	}
# 	geigerframe<-data.frame(newmatrix, row.names= newrownames, stringsAsFactors=F)
# 	geigerframe
# }
