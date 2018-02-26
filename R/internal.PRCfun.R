# internal side functions for use in doRun_prc

getlnTransitionProb<-function(newvalue,meantouse,Fn,priorValues,stdFactor){
		#newvalue = newparticleList[[1]]$startingValues[j],
		#meantouse = prevGenParticleList[[i]]$startingValues[j],
		#Fn=startingPriorsFns[j],
		#priorValues= startingPriorsValues[,j],
		#stdFactor = standardDevFactor
										
	if (Fn=="uniform") {
		sdtouse<-stdFactor*((max(priorValues)-min(priorValues))/sqrt(12))
		#message(paste0("Fn is uniform and sdtouse = ", sdtouse))
	}
	else if (Fn=="exponential") {
		sdtouse<-stdFactor*(1/priorValues[1])
		#message(paste0("Fn is exponential and sdtouse = ", sdtouse))
	}
	else {
		sdtouse<-stdFactor*(priorValues[2])
	}
	#
	lnlocalTransitionProb<-dnorm(newvalue, mean=meantouse, sd=sdtouse,log=TRUE
		) - ((log(1)/pnorm(min(priorValues), mean=meantouse, sd=sdtouse, lower.tail=TRUE, log.p=TRUE))
			* pnorm(max(priorValues), mean=meantouse , sd=sdtouse, lower.tail=FALSE, log.p=TRUE))
	if(length(lnlocalTransitionProb)!=1){
		#message(lnlocalTransitionProb)
		stop("Somehow, multiple lnlocalTransitionProb values produced")
		}
	if (is.nan(lnlocalTransitionProb)) {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
		lnlocalTransitionProb<-.Machine$double.xmin
	}
	if (min(priorValues)==max(priorValues)) {
		lnlocalTransitionProb=log(1)
	}
	if(!is.finite(lnlocalTransitionProb) || is.na(lnlocalTransitionProb)) {
		message(paste0("issue with lnlocalTransitionProb = ",lnlocalTransitionProb))
		}
	return(lnlocalTransitionProb)
	}	

sumLogTranProb<-function(prevGenParticleList
	,newStartingValues, newIntrinsicValues, newExtrinsicValues
	,startingPriorsFns,intrinsicPriorsFns,extrinsicPriorsFns
	,startingPriorsValues,intrinsicPriorsValues,extrinsicPriorsValues
	,standardDevFactor){
	#now get weights, using correction in Sisson et al. 2007
	#
	newWeight=0
	#
	for (i in 1:length(prevGenParticleList)) {
		#
		LLTPstart<-sapply(length(newStartingValues),
			function(j) getlnTransitionProb(newvalue = newStartingValues[j],
				meantouse = prevGenParticleList[[i]]$startingValues[j],
				Fn=startingPriorsFns[j],
				priorValues= startingPriorsValues[,j],
				stdFactor = standardDevFactor))
		LLTPintr<-sapply(length(newIntrinsicValues),
			function(j) getlnTransitionProb(newvalue = newIntrinsicValues[j],
				meantouse = prevGenParticleList[[i]]$intrinsicValues[j],
				Fn=intrinsicPriorsFns[j],
				priorValues= intrinsicPriorsValues[,j],
				stdFactor = standardDevFactor))
		LLTPextr<-sapply(length(newExtrinsicValues),
			function(j) getlnTransitionProb(newvalue = newExtrinsicValues[j],
				meantouse = prevGenParticleList[[i]]$extrinsicValues[j],
				Fn=extrinsicPriorsFns[j],
				priorValues= extrinsicPriorsValues[,j],
				stdFactor = standardDevFactor))
		#
		lnTransitionProb=log(1)
		#
		lnTransitionProb<-lnTransitionProb+sum(LLTPstart)+sum(LLTPintr)+sum(LLTPextr)
		#
		#if(!is.finite(lnTransitionProb) || is.na(lnTransitionProb)) {
		#	warning(paste0("Issue with lnTransitionProb: ",
		#		" lnTransitionProb = ",lnTransitionProb))
		#	}
		#
		newWeight<-newWeight+(prevGenParticleList[[i]]$weight)*exp(lnTransitionProb)
		} #for (i in 1:length(prevGenParticleList)) bracket
	#
	if (!is.finite(newWeight) || is.na(newWeight)) {
		warning(paste0("Possible error; newWeight is ",newWeight))
		}
	return(newWeight)
	}

	