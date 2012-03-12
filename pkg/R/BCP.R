BCP<-function(RealParam, HPD, verbose=F){
#BCP = Bayesian Coverage Probability
#RealParam can be "RealParams" from doSimulation or a vector c(x1, x2, ...)
#Calculates what percent of the time the real parameter falls into the 95% credible interval
if(length(RealParam) != dim(HPD[[1]])[1]){ warning("RealParams and HPD do not match")}
Covered<-matrix(nrow=length(HPD), ncol=length(RealParam))
colnames(Covered)<-rownames(HPD[[1]])
	for(i in 1:length(HPD)){
		for(j in 1:length(RealParam)){
			if(!is.na(HPD[[i]][j,4])) { #keep only parameters that vary (ie, not fixed)
				if(HPD[[i]][j,4] > RealParam[j] && RealParam[j] > HPD[[i]][j,3]) {
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