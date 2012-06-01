
loadSimulations<-function(filenames) {
	results<-c() #check type
	for (i in sequence(length(filenames))) {
		load(filenames[i]) 
		results<-rbind(results,trueFreeValuesANDSummaryValues)
	}
	return(results)
}