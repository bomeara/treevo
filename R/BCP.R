#BCP = Bayesian Coverage Probability
#RealParam can be "RealParams$vector" from doSimulation c(x1, x2, ...) or a list
#HPD should be a list of HPD output from different runs
#Calculates what percent of the time the real parameter falls into the HPD
BCP<-function(RealParam, HPD, verbose=F){  #alter so it can be a single or a list of realParmas
	if(class(RealParam)=="numeric"){
		rps<-vector("list", length=length(HPD))
		for (i in 1: length(HPD)){
			rps[[i]]<-RealParam
		}
	}
	rps<-RealParam
	#if(length(RealParam) != dim(HPD[[1]])[1]){ warning("RealParams and HPD do not match")}  #need something like this, but it will have to be after changing the RealParam to take a list 
	Covered<-matrix(nrow=length(HPD), ncol=length(rps[[i]]))
	colnames(Covered)<-rownames(HPD[[1]])
	for(i in 1:length(HPD)){
		for(j in 1:length(rps[[i]])){
			if(!is.na(HPD[[i]][j,4])) { #keep only parameters that vary (ie, not fixed)
				if(HPD[[i]][j,4] > rps[[i]][j] && rps[[i]][j] > HPD[[i]][j,3]) {
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