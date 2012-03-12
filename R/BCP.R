BCP<-function(RealParam, FinalParamPredictions, verbose=F){
#BCP = Bayesian Coverage Probability
#RealParam can be "RealParams" from doSimulation or a vector c(x1, x2, ...)
#Calculates what percent of the time the real parameter falls into the 95% credible interval
if(length(RealParam) != dim(FinalParamPredictions[[1]])[1]){ warning("RealParams and FinalParamPredictions do not match")}
Covered<-matrix(nrow=length(FinalParamPredictions), ncol=length(RealParam))
colnames(Covered)<-rownames(FinalParamPredictions[[1]])
	for(i in 1:length(FinalParamPredictions)){
		for(j in 1:length(RealParam)){
			if(FinalParamPredictions[[i]][j,4]-FinalParamPredictions[[i]][j,3] != 0) { #keep only parameters that vary (ie, not fixed)
				if(FinalParamPredictions[[i]][j,4] > RealParam[j] && RealParam[j] > FinalParamPredictions[[i]][j,3]) {
					Covered[i,j]<-1
				}
				else{
					Covered[i,j]<-0
				}
			}
		}
	}
	CoverProb<-apply(Covered, 2, mean)
	if(verbose){
		CoverProb<-vector("list")
		CoverProb$byRun<-Covered
		CoverProb$BCP<-apply(Covered, 2, mean)
	}
	return(CoverProb)
}