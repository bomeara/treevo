
presentABCOutput<-function(ABCOutput, plot=FALSE, priors=c(), truth=c()) {
	library(Hmisc)
	lastGen<-ABCOutput[which(ABCOutput$generation==max(ABCOutput$generation)), ]
	finalParticles<-lastGen[which(lastGen$weight>0), ]
	nParams<-dim(finalParticles)[2]-6
	nParticles<-dim(finalParticles)[1]
	resultsMatrix<-matrix(nrow=13, ncol=0)
	
	for (variable in 1:nParams) {
		resultsMatrix<-cbind(resultsMatrix, wtd.quantile(finalParticles[, 6+variable], weights= nParticles*finalParticles[, 6], probs=c(0, 0.001, 0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 0.999, 1)))
	}
	
	if (plot==TRUE) {
		par(mfcol=c(1, nParams))
		for (variable in 1: nParams) {
			cat("prior index from ", 1+(variable-1)*2, " to ", 2+(variable-1)*2, "\n")
			densityResults<-density(finalParticles[, 6+variable], weights= nParticles*finalParticles[, 6]/sum(nParticles*finalParticles[, 6]), from=min(c(priors[1+(variable-1)*2], priors[2+(variable-1)*2])), to=max(c(priors[1+(variable-1)*2], priors[2+(variable-1)*2])))
			plot(x=c(priors[1+(variable-1)*2], priors[2+(variable-1)*2]), y=c(0, max(densityResults$y)), yaxt="n", xlab=(names(finalParticles)[6+variable]), type="n", ylab="", main="", open )
			lines(densityResults)
				lines(x=c(truth[variable], truth[variable]), y=c(0, max(densityResults$y)), lty=3)
		}
	}
	results<-data.frame(resultsMatrix)
	return(results)
}
