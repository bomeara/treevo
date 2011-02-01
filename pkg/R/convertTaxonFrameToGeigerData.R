

convertTaxonFrameToGeigerData<-function(taxonframe, phy) {
	ntax<-dim(taxonframe)[1]
	newmatrix<-matrix(data=taxonframe[, 4:(dim(taxonframe)[2])], nrow= ntax)
	newrownames<-c(rep(0, ntax))
	for (i in 1:ntax) {
		newrownames[i]<-phy$tip.label[(taxonframe$taxonid[i])]
	}
	geigerframe<-data.frame(newmatrix, row.names= newrownames, stringsAsFactors=F)
	geigerframe
}
