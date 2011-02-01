summarizeTaxonStates<-function(taxa) {
#print("in summarizeTaxonStates")
	statesvector<-c()
	taxonid<-c()
	taxonname<-c()
	taxontimesincespeciation<-c()
	for (i in 1:length(taxa)) {
#print(i)
#		print(taxa)
		statesvector<-c(statesvector, states(taxa[[i]]))
#print("statesvector")
		taxonid<-c(taxonid, id(taxa[[i]]) )
#print("taxonid")
		taxonname<-c(taxonname, name(taxa[[i]]) )
#print("taxonname")
		taxontimesincespeciation<-c(taxontimesincespeciation, timeSinceSpeciation(taxa[[i]]))
#print("finished ", i)
	}
	statesmatrix<-matrix(statesvector, ncol=length(states(taxa[[1]])), byrow=TRUE) #each row represents one taxon
	taxonframe<-data.frame(taxonid, taxonname, taxontimesincespeciation, statesmatrix)
	return(taxonframe)
}
