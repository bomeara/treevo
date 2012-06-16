summarizeTaxonStates<-function(taxa) {
#print("in summarizeTaxonStates")
	statesvector<-c()
	taxonid<-c()
	taxonname<-c()
	taxontimesincespeciation<-c()
	for (i in 1:length(taxa)) {
		statesvector<-c(statesvector, taxa[[i]]$states)
		taxonid<-c(taxonid, taxa[[i]]$id )
		taxonname<-c(taxonname, taxa[[i]]$name )
		taxontimesincespeciation<-c(taxontimesincespeciation, taxa[[i]]$timeSinceSpeciation)
	}
	statesmatrix<-matrix(statesvector, ncol=length(taxa[[1]]$states), byrow=TRUE) #each row represents one taxon
	taxonframe<-data.frame(taxonid, taxonname, taxontimesincespeciation, statesmatrix)
	return(taxonframe)
}
