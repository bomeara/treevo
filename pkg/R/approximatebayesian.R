#minimal ABC classes, pulled from  abctaxon.R, abc.R, abcparticle.R

library(geiger)
library(car)
library(mixOmics)
library(Hmisc)


#class for particles in ABC-PRC
setClass(
	Class="abcparticle", 
	representation=representation(
		id="numeric", 
		generation="numeric", 
		weight="numeric", 
		distance="numeric", 
		startingStates="numeric", 
		intrinsicValues="numeric", 
		extrinsicValues="numeric", 
		kernel="numeric"
	)
)

setGeneric("id", function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id", "abcparticle", function(object) {
	object@id
})

setReplaceMethod("id", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@id <- value
	object
}
)

setGeneric("setId", function(object, value) {
	standardGeneric("setId")
})

setMethod("setId", "abcparticle", function(object, value) {
	object@id<-value
	object
})

setGeneric("getId", function(object) {
	standardGeneric("getId")
})

setMethod("getId", "abcparticle", function(object) {
	object@id
})

setGeneric("setWeight", function(object, value) {
	standardGeneric("setWeight")
})

setMethod("setWeight", "abcparticle", function(object, value) {
	object@weight<-value
	object
})

setGeneric("getWeight", function(object) {
	standardGeneric("getWeight")
})

setMethod("getWeight", "abcparticle", function(object) {
	object@weight
})

setGeneric("generation", function(object) {
	standardGeneric("generation")
})

setGeneric("generation<-", function(object, value) {
	standardGeneric("generation<-")
})


setMethod("generation", "abcparticle", function(object) {
	object@generation
})

setReplaceMethod("generation", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@generation <- value
	object
}
)

setGeneric("weight", function(object) {
	standardGeneric("weight")
})

setGeneric("weight<-", function(object, value) {
	standardGeneric("weight<-")
})


setMethod("weight", "abcparticle", function(object) {
	object@weight
})

setReplaceMethod("weight", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@weight <- value
	object
}
)

setGeneric("distance", function(object) {
	standardGeneric("distance")
})

setGeneric("distance<-", function(object, value) {
	standardGeneric("distance<-")
})


setMethod("distance", "abcparticle", function(object) {
	object@distance
})

setReplaceMethod("distance", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@distance <- value
	object
}
)

setGeneric("setDistance", function(object, value) {
	standardGeneric("setDistance")
})

setMethod("setDistance", "abcparticle", function(object, value) {
	object@distance<-value
	object
})

setGeneric("kernel", function(object) {
	standardGeneric("kernel")
})

setGeneric("kernel<-", function(object, value) {
	standardGeneric("kernel<-")
})


setMethod("kernel", "abcparticle", function(object) {
	object@kernel
})

setReplaceMethod("kernel", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@kernel <- value
	object
}
)


setGeneric("startingStates", function(object) {
	standardGeneric("startingStates")
})

setGeneric("startingStates<-", function(object, value) {
	standardGeneric("startingStates<-")
})


setMethod("startingStates", "abcparticle", function(object) {
	object@startingStates
})

setReplaceMethod("startingStates", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@startingStates <- value
	object
}
)

setGeneric("intrinsicValues", function(object) {
	standardGeneric("intrinsicValues")
})

setGeneric("intrinsicValues<-", function(object, value) {
	standardGeneric("intrinsicValues<-")
})


setMethod("intrinsicValues", "abcparticle", function(object) {
	object@intrinsicValues
})

setReplaceMethod("intrinsicValues", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@intrinsicValues <- value
	object
}
)

setGeneric("extrinsicValues", function(object) {
	standardGeneric("extrinsicValues")
})

setGeneric("extrinsicValues<-", function(object, value) {
	standardGeneric("extrinsicValues<-")
})


setMethod("extrinsicValues", "abcparticle", function(object) {
	object@extrinsicValues
})

setReplaceMethod("extrinsicValues", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@extrinsicValues <- value
	object
}
)

setGeneric("initializeStatesFromMatrices", function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns){standardGeneric("initializeStatesFromMatrices")})

setMethod("initializeStatesFromMatrices", signature="abcparticle", definition=function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns) {
		#cat("startingPriorsValues\n")
		#print(startingPriorsValues)
		#cat("\nintrinsicPriorsValues\n")
		#print(intrinsicPriorsValues)
		#cat("\nextrinsicPriorsValues\n")
		#print(extrinsicPriorsValues)
		#cat("\ndim(startingPriorsValues)[2]=", dim(startingPriorsValues)[2], " dim(intrinsicPriorsValues)[2]=", dim(intrinsicPriorsValues)[2], " dim(extrinsicPriorsValues)[2])=", dim(extrinsicPriorsValues)[2], "\n")
		x@startingStates<-vector(mode="numeric", length=dim(startingPriorsValues)[2])
		x@intrinsicValues<-vector(mode="numeric", length=dim(intrinsicPriorsValues)[2])
		x@extrinsicValues <-vector(mode="numeric", length=dim(extrinsicPriorsValues)[2])
		#cat("\n initializeStatesFromMatrices\n")
		for (j in 1:dim(startingPriorsValues)[2]) {
			x@startingStates[j]=pullFromPrior(priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
			#cat("starting states j=", j, " x@startingStates[j] =", x@startingStates[j], "\n")
		}
		for (j in 1:dim(intrinsicPriorsValues)[2]) {
			x@intrinsicStates[j]=pullFromPrior(priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
			#cat("intrinsicValues j=", j, " x@ intrinsicValues[j] =", x@intrinsicValues[j], "\n")

		}
		for (j in 1:dim(extrinsicPriorsValues)[2]) {
			x@extrinsicStates[j]=pullFromPrior(priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
			#cat("extrinsicValues j=", j, " x@ extrinsicValues[j] =", x@extrinsicValues[j], "\n")
		}
		#print(x)
		x
		}
)

pullFromPrior<-function(priorValues, priorFn) {
	#fixed, uniform, normal, lognormal, gamma
	x<-NA
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed") {
		x<-priorValues[1]
	}
	else if (priorFn=="uniform") {
		x<-runif(n=1,min=min(priorValues), max=max(priorValues))
	}
	else if (priorFn=="normal") {
		x<-rnorm(n=1, mean=priorValues[1], sd=priorValues[2])
	}
	else if (priorFn=="lognormal") {
		x<-rlnorm(n=1, meanlog = priorValues[1], sdlog = priorValues[2])
	}
	else if (priorFn=="gamma") {
		x<-rgamma(n=1, shape=priorValues[1], scale = priorValues[2])
	}
	else if (priorFn=="exponential") {
		x<-rexp(n=1, rate=priorValues[1])
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}
	x
}

meanFromPrior<-function(priorValues, priorFn) {
	#fixed, uniform, normal, lognormal, gamma
	x<-NA
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed") {
		x<-priorValues[1]
	}
	else if (priorFn=="uniform") {
		x<-mean(priorValues)
	}
	else if (priorFn=="normal") {
		x<-priorValues[1]
	}
	else if (priorFn=="lognormal") {
		x<- priorValues[1]
	}
	else if (priorFn=="gamma") {
		x<-priorValues[1]*priorValues[2]
	}
	else if (priorFn=="exponential") {
		x<-1.0/priorValues[1]
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}
	x
}


setGeneric("mutateStates", function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor){standardGeneric("mutateStates")})

setMethod("mutateStates", signature="abcparticle", definition=function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor) {
	#dput(x)
	#dput(x@startingStates)
	#typeof(x)
	#typeof(x@startingStates)
		replacementVector<-rep(NA, length(x@startingStates))
		for (j in 1:length(x@startingStates)) {
			replacementVector[j]<-mutateState(startingState=x@startingStates[j], standardDevFactor=standardDevFactor, priorValues=x@startingPriorsValues[, j], priorFn=x@startingPriorsFns[j])
		}
		x@startingStates<-replacementVector
		
		replacementVector<-rep(NA, length(x@intrinsicStates))
		for (j in 1:length(x@intrinsicStates)) {
			replacementVector[j]<-mutateState(intrinsicState=x@intrinsicStates[j], standardDevFactor=standardDevFactor, priorValues=x@intrinsicPriorsValues[, j], priorFn=x@intrinsicPriorsFns[j])
		}
		x@intrinsicStates<-replacementVector

		replacementVector<-rep(NA, length(x@extrinsicStates))
		for (j in 1:length(x@extrinsicStates)) {
			replacementVector[j]<-mutateState(extrinsicState=x@extrinsicStates[j], standardDevFactor=standardDevFactor, priorValues=x@extrinsicPriorsValues[, j], priorFn=x@extrinsicPriorsFns[j])
		}
		x@extrinsicStates<-replacementVector
		
		#kernel(x)<-lnTransitionProb
		x
		}
)

mutateState<-function(startingState, standardDevFactor, priorValues, priorFn) {
	newState<-NA
	minBound=-Inf
	maxBound=Inf
	validNewState<-False
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed" || priorFn=="uniform") {
		minBound<-min(priorValues)
		maxBound<-max(priorValues)
	}
	else if (priorFn=="lognormal" || priorFn=="gamma" || priorFn=="exponential") {
		minBound<-0
	}
	
	sdToUse<-standardDevFactor
	if (priorFn=="fixed") {
		sdToUse<-0
	}
	else if (priorFn=="uniform") {
		sdToUse<-standardDevFactor*range(priorValues)
	}
	else if (priorFn=="normal") {
		sdToUse<-standardDevFactor*priorValues[2]
	}
	else if (priorFn=="lognormal") {
		sdToUse<-standardDevFactor*priorValues[2]
	}
	else if (priorFn=="gamma") {
		sdToUse<-standardDevFactor*sqrt(priorValues[1]*priorValues[2]*priorValues[2])
	}
	else if (priorFn=="exponential") {
		sdToUse<-standardDevFactor/priorValues[1]
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}
	
	while(!validNewState) {
		newState<-rnorm(n=1, mean=startingState, sd=sdToUse)
	}
	newState
}
	

setGeneric("simulateTips", function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep){standardGeneric("simulateTips")})
	
setMethod("simulateTips", signature="abcparticle", definition=function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, x@startingStates, x@intrinsicValues, x@extrinsicValues, timeStep), phy)
	return(newtips)
	})	
	
setClass(
Class="abctaxon", 
representation=representation(
id="numeric", 
name="character", 
states="numeric", 
nextstates="numeric", 
timeSinceSpeciation="numeric"
), 
prototype(id=0, name="taxon0", timeSinceSpeciation=0)
)

setGeneric("id", function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id", "abctaxon", function(object) {
	object@id
})

setReplaceMethod("id", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@id <- value
	object
}
)

setGeneric("name", function(object) {
	standardGeneric("name")
})

setGeneric("name<-", function(object, value) {
	standardGeneric("name<-")
})


setMethod("name", "abctaxon", function(object) {
	object@name
})

setReplaceMethod("name", signature(object="abctaxon", value="character") , 
function(object, value) {
	object@name <- value
	object
}
)



setGeneric("states", function(object) {
	standardGeneric("states")
})

setGeneric("states<-", function(object, value) {
	standardGeneric("states<-")
})


setMethod("states", "abctaxon", function(object) {
	object@states
})

setReplaceMethod("states", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@states <- value
	object
}
)



setGeneric("nextstates", function(object) {
	standardGeneric("nextstates")
})

setGeneric("nextstates<-", function(object, value) {
	standardGeneric("nextstates<-")
})


setMethod("nextstates", "abctaxon", function(object) {
	object@nextstates
})

setReplaceMethod("nextstates", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@nextstates <- value
	object
}
)


setGeneric("timeSinceSpeciation", function(object) {
	standardGeneric("timeSinceSpeciation")
})

setGeneric("timeSinceSpeciation<-", function(object, value) {
	standardGeneric("timeSinceSpeciation<-")
})


setMethod("timeSinceSpeciation", "abctaxon", function(object) {
	object@timeSinceSpeciation
})

setReplaceMethod("timeSinceSpeciation", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@timeSinceSpeciation <- value
	object
}
)

getSimulationSplits<-function(phy) {
	branchingTimes<-sort(branching.times(phy), decreasing=T)
	branchingTimesNames<-names(branchingTimes)
	ancestorVector<-c()
	descendant1Vector<-c()
	descendant2Vector<-c()
	for (i in 1:length(branchingTimes)) {
		relationshipVector<-phy$edge[phy$edge[, 1]==branchingTimesNames[i]]
		ancestorVector<-c(ancestorVector, relationshipVector[2])
		descendant1Vector<-c(descendant1Vector, relationshipVector[3])
		descendant2Vector<-c(descendant2Vector, relationshipVector[4])
	}
	
	return(data.frame(branchingTimes, ancestorVector, descendant1Vector, descendant2Vector))
}

doSimulation<-function(splits, intrinsicFn, extrinsicFn, startingStates, intrinsicValues, extrinsicValues, timeStep) {
	numberofsteps<-floor(splits[1, 1]/timeStep)
	mininterval<-min(splits[1:(dim(splits)[1]-1), 1]-splits[2:(dim(splits)[1]), 1])
	if (numberofsteps<1000) {
#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", splits[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
#warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more. I would suggest decreasing timeStep to no more than ", mininterval/50))
	}
#initial setup
	timefrompresent=splits[1, 1]
	taxa<-c(new("abctaxon", id=splits[1, 3], states=startingStates), new("abctaxon", id=splits[1, 4], states=startingStates))
	splits<-splits[2:dim(splits)[1], ] #pop off top value
	
#start running
	while(timefrompresent>0) {
#print(timefrompresent)
#speciation if needed
		while ((timefrompresent-timeStep)<=splits[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
			originallength<-length(taxa)
			taxontodelete<-Inf
			for (i in 1:originallength) {
				if (id(taxa[[i]])==splits[1, 2]) {
					taxontodelete<-i
					taxa<-c(taxa, taxa[[i]], taxa[[i]])
					id(taxa[[originallength+1]])<-splits[1, 3]
					timeSinceSpeciation(taxa[[originallength+1]])<-0
					id(taxa[[originallength+2]])<-splits[1, 4]
					timeSinceSpeciation(taxa[[originallength+2]])<-0
				}
			}
#cat("taxontodelete = ", taxontodelete)
			taxa<-taxa[-1*taxontodelete]
			if(dim(splits)[1]>1) {
				splits<-splits[2:(dim(splits)[1]), ] #pop off top value
			}
			else {
				splits[1, ]<-c(-1, 0, 0, 0)
			}
#print("------------------- speciation -------------------")
#print(taxa)
#summarizeTaxonStates(taxa)
		}
#trait evolution step
		for (i in 1:length(taxa)) {
			otherstatesvector<-c()
			for (j in 1:length(taxa)) {
				if (j!=i) {
					otherstatesvector<-c(otherstatesvector, states(taxa[[j]]))
				}
			}
#print(taxa)
#print(length(otherstatesvector))
			otherstatesmatrix<-matrix(otherstatesvector, ncol=length(states(taxa[[i]])), byrow=T) #each row represents one taxon
			newvalues<-states(taxa[[i]])+intrinsicFn(params=intrinsicValues, states=states(taxa[[i]]), timefrompresent =timefrompresent)+extrinsicFn(params=extrinsicValues, selfstates=states(taxa[[i]]), otherstates=otherstatesmatrix, timefrompresent =timefrompresent)
			nextstates(taxa[[i]])<-newvalues
		}
		for (i in 1:length(taxa)) {
#print("\nbefore\n")
#print(taxa[[i]])
			states(taxa[[i]])<-nextstates(taxa[[i]])
#print("\nafter\n")
#print(taxa[[i]])
			
		}
#print("------------------- step -------------------")
#print(taxa)
#summarizeTaxonStates(taxa)
		
		timefrompresent<-timefrompresent-timeStep
		for (i in 1:length(taxa)) {
			timeSinceSpeciation(taxa[[i]])<-timeSinceSpeciation(taxa[[i]])+timeStep
		}
	}
	return(summarizeTaxonStates(taxa))
}

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
	statesmatrix<-matrix(statesvector, ncol=length(states(taxa[[1]])), byrow=T) #each row represents one taxon
	taxonframe<-data.frame(taxonid, taxonname, taxontimesincespeciation, statesmatrix)
	return(taxonframe)
}

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


summaryStatsLong<-function(phy, data, todo=c()) {
	sink(file="/dev/null")
	if (length(todo)==0) {
		todo=rep(1, 22+dim(data)[1]) #by default, include everything -- the 22 summary stats and the raw tip data
	}
	
	#thought here: include brown<-try(fitContinouous()), store brown, then do try(brown$lnl), etc. Faster than calling each fn
	brown<-try(fitContinuous(phy=phy, data=data[max(todo[1:3])], model="BM")[[1]]) #will only run if want to do brownian summary stats
	#brown.lnl<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[1]], model="BM")[[1]]$lnl)) #if todo[i]==0, will cause an error right away, saving on computation time
	brown.lnl<-as.numeric(try(brown$lnl/todo[1])) #divide by zero so we get Inf if we don't want that summary stat
	#brown.beta<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[2]], model="BM")[[1]]$beta))
	brown.beta <-as.numeric(try(brown$beta/todo[2]))
	#brown.aic<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[3]], model="BM")[[1]]$aic))
	brown.aic <-as.numeric(try(brown$aic/todo[3]))

	lambda<-try(fitContinuous(phy=phy, data=data[max(todo[4:7])], model="lambda")[[1]])
	#lambda.lnl<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[4]], model="lambda")[[1]]$lnl))
	lambda.lnl <-as.numeric(try(lambda$lnl/todo[4]))
	#lambda.beta<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[5]], model="lambda")[[1]]$beta))
	lambda.beta <-as.numeric(try(lambda$beta/todo[5]))
	#lambda.lambda<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[6]], model="lambda")[[1]]$lambda))
	lambda.lambda <-as.numeric(try(lambda$lambda/todo[6]))
	#lambda.aic<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[7]], model="lambda")[[1]]$aic))
	lambda.aic <-as.numeric(try(lambda$aic/todo[7]))

	delta<-try(fitContinuous(phy=phy, data=data[max(todo[8:11])], model="delta")[[1]])
	#delta.lnl<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[8]], model="delta")[[1]]$lnl))
	delta.lnl <-as.numeric(try(delta$lnl/todo[8]))
	#delta.beta<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[9]], model="delta")[[1]]$beta))
	delta.beta <-as.numeric(try(delta$beta/todo[9]))
	#delta.delta<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[10]], model="delta")[[1]]$delta))
	delta.delta <-as.numeric(try(delta$delta/todo[10]))
	#delta.aic<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[11]], model="delta")[[1]]$aic))
	delta.aic <-as.numeric(try(delta$aic/todo[11]))
	
	ou<-try(fitContinuous(phy=phy, data=data[max(todo[12:15])], model="OU")[[1]])
	#ou.lnl<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[12]], model="OU")[[1]]$lnl))
	ou.lnl <-as.numeric(try(ou$lnl/todo[12]))
	#ou.beta<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[13]], model="OU")[[1]]$beta))
	ou.beta <-as.numeric(try(ou$beta/todo[13]))
	#ou.alpha<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[14]], model="OU")[[1]]$alpha))
	ou.alpha <-as.numeric(try(ou$alpha/todo[14]))
	#ou.aic<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[15]], model="OU")[[1]]$aic))
	ou.aic <-as.numeric(try(ou$aic/todo[15]))
	
	white<-try(fitContinuous(phy=phy, data=data[max(todo[16:17])], model="white")[[1]])
	#white.lnl<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[16]], model="white")[[1]]$lnl))
	white.lnl<-as.numeric(try(white$lnl/todo[16]))
	#white.aic<-as.numeric(try(fitContinuous(phy=phy, data=data[todo[17]], model="white")[[1]]$aic))
	white.aic<-as.numeric(try(white$aic/todo[17]))
	
	
	raw.mean<-as.numeric(try(mean(data[todo[18]])))
	raw.max<-as.numeric(try(max(data[todo[19]])))
	raw.min<-as.numeric(try(min(data[todo[20]])))
	raw.var<-as.numeric(try(var(data[todo[21]])))
	raw.median<-as.numeric(try(median(data[todo[22], ])))
	#cat("summaryStatsLong")
	summarystats<-c(brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic, raw.mean, raw.max, raw.min, raw.var, raw.median, data[[1]] )
	#cat("\n summaryStatsLong summarystats1\n")
	#print(summarystats)
	summarystats[which(todo==0)]<-NA
	#cat("\n summaryStatsLong summarystats2\n")
	#print(summarystats)
	summarystats[which(is.finite(summarystats)==FALSE)]<-NA
	#cat("\n summaryStatsLong summarystats3\n")
	#print(summarystats)
	while(sink.number()>0) {sink()}
	summarystats
	
	
}

boxcoxplsSummary<-function(todo, summaryValues, prunedPlsResult, boxcoxLambda, boxcoxAddition) {
	for (summaryValueIndex in 1:length(summaryValues)) {
		summaryValues[summaryValueIndex]<-(summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
	}
	prunedSummaryValues<-summaryValues[which(todo>0)]
	predictResult<-(predict(prunedPlsResult, matrix(prunedSummaryValues, nrow=1))$predict[, , 1])
	predictResult
}


#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), timeStep, toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, nrepSim=100, plot=T, vipthresh=0.8, epsilonProportion=0.2, epsilonMultiplier=0.5, nStepsPRC=4) {
	#print("in do run")
	splits<-getSimulationSplits(phy) #initialize this info


	#figure out number of free params
	numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	numberParametersFree<-numberParametersTotal
	numberParametersStarting<-0
	numberParametersIntrinsic<-0
	numberParametersExtrinsic<-0
	freevariables<-matrix(data=NA, nrow=2, ncol=0)
	titlevector<-c()
	freevector<-c()
	
		
	for (i in 1:dim(startingPriorsValues)[2]) {
		priorFn<-match.arg(arg=startingPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersStarting<-numberParametersStarting+1
			freevariables<-cbind(freevariables, startingPriorsValues[, i])
			titlevector <-c(titlevector, paste("Starting", numberParametersStarting))
			freevector<-c(freevector, TRUE)
		}
	}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=intrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersIntrinsic<-numberParametersIntrinsic+1
			freevariables<-cbind(freevariables, intrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Intrinsic", numberParametersIntrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=extrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersExtrinsic<-numberParametersExtrinsic+1
			freevariables<-cbind(freevariables, extrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Extrinsic", numberParametersExtrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	

	#initialize guesses, if needed
	if (length(startingStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		startingStatesGuess<-rep(NA,length(startingPriorsFns))
		for (i in 1:length(startingPriorsFns)) {
			startingStatesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
		}
	}
	if (length(intrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		intrinsicStatesGuess<-rep(NA,length(intrinsicPriorsFns))
		for (i in 1:length(intrinsicPriorsFns)) {
			intrinsicStatesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
		}
	}
	if (length(extrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		extrinsicStatesGuess<-rep(NA,length(extrinsicPriorsFns))
		for (i in 1:length(extrinsicPriorsFns)) {
			extrinsicStatesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
		}
	}
	
	
	#----------------- Find best set of summary stats to use for this problem. (start) -----------------
	#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data
	for (simIndex in 1:nrepSim) {
		trueStarting<-rep(NaN, dim(startingPriorsValues)[2])
		trueIntrinsic<-rep(NaN, dim(intrinsicPriorsValues)[2])
		trueExtrinsic<-rep(NaN, dim(extrinsicPriorsValues)[2])
		for (j in 1:dim(startingPriorsValues)[2]) {
			trueStarting[j]=pullFromPrior(startingPriorsValues[,j],startingPriorsFns[j])
		}
		for (j in 1:dim(intrinsicPriorsValues)[2]) {
			trueIntrinsic[j]=pullFromPrior(intrinsicPriorsValues[,j],intrinsicPriorsFns[j])
		}
		for (j in 1:dim(extrinsicPriorsValues)[2]) {
			trueExtrinsic[j]=pullFromPrior(extrinsicPriorsValues[,j],extrinsicPriorsFns[j])
		}
		trueInitial<-c(trueStarting, trueIntrinsic, trueExtrinsic)
		trueFreeValues<-rbind(trueFreeValues, trueInitial[freevector])
		#cat("summaryValues\n")
		#print(summaryValues)
		summaryValues<-rbind(summaryValues, summaryStatsLong(phy, convertTaxonFrameToGeigerData (doSimulation(splits=splits, intrinsicFn= intrinsicFn, extrinsicFn= extrinsicFn, startingStates= trueStarting, intrinsicValues= trueIntrinsic, extrinsicValues= trueExtrinsic, timeStep=timeStep), phy)))
		while(sink.number()>0) {sink()}
	}
	
	library("car")
	#now put this into the boxcox function to get best lambda for each summary stat
	boxcoxLambda<-rep(NA, dim(summaryValues)[2])
	boxcoxAddition<-rep(NA, dim(summaryValues)[2])
	for (summaryValueIndex in 1:dim(summaryValues)[2]) {
		boxcoxAddition[summaryValueIndex]<-0
		lowValue<-min(summaryValues[, summaryValueIndex])-4*sd(summaryValues[, summaryValueIndex])
		if (lowValue<=0) {
			boxcoxAddition[summaryValueIndex]<-4*abs(lowValue) #just for some protection against low values, since box.cox needs non-negative values
		}
		cat("\nsummary values[", summaryValueIndex, "] = ")
		print(summaryValues[, summaryValueIndex])
		summary<-summaryValues[, summaryValueIndex]+boxcoxAddition[summaryValueIndex]
		boxcoxLambda[summaryValueIndex]<-1
		if(sd(summaryValues[, summaryValueIndex])>0) { #box.cox fails if all values are identical
			newLambda<-as.numeric(try(box.cox.powers.hacked(summary)$lambda))
			if (!is.na(newLambda)) {
				boxcoxLambda[summaryValueIndex]<-newLambda
			}
		}
		summaryValues[, summaryValueIndex]<-summary^boxcoxLambda[summaryValueIndex]
	}
	
	
	#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
	library("mixOmics")
	plsResult<-pls(Y=trueFreeValues, X=summaryValues)
	vipResult<-vip(plsResult)
	todo<-rep(1, dim(summaryValues)[2])#initialize the vector that indicates which summary stats to include
	for (summaryIndex in 1:dim(summaryValues)[2]) {
		if (max(vipResult[summaryIndex, ])<vipthresh) {
				todo[summaryIndex]<-0 #exclude this summary stat, because it's too unimportant
		}	
	}
	while(sink.number()>0) {sink()}
	#cat("todo", "\n")
	#print(todo)
	
	prunedSummaryValues<-summaryValues[, which(todo>0)]
	#cat("prunedSummaryValues", "\n")
	#print(prunedSummaryValues)

	prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)
	#cat("prunedPlsResult", "\n")
	#print(prunedPlsResult)

	#----------------- Find best set of summary stats to use for this problem. (end) -----------------
	
	#----------------- Find distribution of distances (start) ----------------------
	predictResult<-(predict(prunedPlsResult, prunedSummaryValues)$predict[, , 1])
	#cat("predictResult", "\n")
	#print(predictResult)

	distanceVector<-rep(NA, dim(predictResult)[1])
	for (simulationIndex in 1:dim(predictResult)[1]) {
			distanceVector[simulationIndex]<-dist(matrix(c(trueFreeValues[simulationIndex, ], predictResult[simulationIndex, ]), nrow=2, byrow=T))[1]
	}
	#cat("distanceVector", "\n")
	#print(distanceVector)
	densityDistanceVector<-density(distanceVector)
	#plot(densityDistanceVector)
	epsilonDistance<-quantile(distanceVector, probs=epsilonProportion) #this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected 
	#lines(x=c(epsilonDistance, epsilonDistance), y=c(0, max(densityDistanceVector$y)), lty="dotted")
	toleranceVector<-rep(epsilonDistance, nStepsPRC)
	for (step in 2:nStepsPRC) {
		toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
	}
	#----------------- Find distribution of distances (end) ---------------------
	
	#------------------ ABC-PRC (start) ------------------
	#do not forget to use boxcoxLambda, and prunedPlsResult when computing distances
	nameVector<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	if (plot) {
		plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 5*max(toleranceVector)), type="n")
	}
	for (i in 1:dim(startingPriorsValues)[2]) {
		nameVector<-append(nameVector, paste("StartingStates", i, sep=""))
	}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		nameVector<-append(nameVector, paste("IntrinsicValue", i, sep=""))
	}
	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		nameVector<-append(nameVector, paste("ExtrinsicValue", i, sep=""))
	}
	particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
	particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
	particleDistance=rep(NA, numParticles)
	particle<-1
	attempts<-0
	particleDataFrame<-data.frame()
	cat("successes", "attempts", "expected number of attempts required\n")
	particleVector<-c()
	originalSummaryStats<-boxcoxplsSummary(todo, summaryStatsLong(phy, traits, todo), prunedPlsResult, boxcoxLambda, boxcoxAddition)
	#cat("originalSummaryStats\n")
	#print(originalSummaryStats)
	while (particle<=numParticles) {
		attempts<-attempts+1
		
		newparticleVector<-c(new("abcparticle", id=particle, generation=1, weight=0))
		newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]], startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues)
		#cat("\nextrinsicVector\n")
		#print(extrinsicValues(newparticleVector[[1]]))
		#cat("\nintrinsicVector\n")
		#print(intrinsicValues(newparticleVector[[1]]))

		newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=T))[1])
		if (is.na(distance(newparticleVector[[1]]))) {
			while(sink.number()>0) {sink()}
			warning("distance(newparticleVector[[1]]) = NA")
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
		else if (is.na(toleranceVector[1])) {
			while(sink.number()>0) {sink()}
			warning("toleranceVector[1] = NA")
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
		else if (distance(newparticleVector[[1]])<toleranceVector[1]) {
			newparticleVector[[1]]<-setId(newparticleVector[[1]], particle)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 1/numParticles)
			particleWeights[particle]<-1/numParticles
			particle<-particle+1
			particleVector<-append(particleVector, newparticleVector)
		}
		else {
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
		while(sink.number()>0) {sink()}
		#print(newparticleVector)
		vectorForDataFrame<-c(1, attempts, getId(newparticleVector[[1]]), 0, distance(newparticleVector[[1]]), getWeight(newparticleVector[[1]]), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]))
#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingStates = ", length(startingStates), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", distance(newparticleVector[[1]]), "\nweight = ", getWeight(newparticleVector[[1]]), "\n", vectorForDataFrame, "\n")
		particleDataFrame<-rbind(particleDataFrame, data.frame(rbind(vectorForDataFrame)))
		cat(particle-1, attempts, floor(numParticles*attempts/particle), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), distance(newparticleVector[[1]]), "\n")
		
	}
	
	for (dataGenerationStep in 2:length(toleranceVector)) {
		cat("\n\n\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n\n\n")
		particleWeights<-particleWeights/(sum(particleWeights)) #normalize to one
		cat("particleWeights\n", particleWeights, "\n\n")
		oldParticleVector<-particleVector
		oldParticleWeights<-particleWeights
		particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
		particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
		particleDistance=rep(NA, numParticles)
		particle<-1
		attempts<-0
		cat("successes", "attempts", "expected number of attempts required\n")
		particleVector<-c()
		weightScaling=0;
		while (particle<=numParticles) {
			attempts<-attempts+1
			particleToSelect<-which.max(as.vector(rmultinom(1, size = 1, prob=oldParticleWeights)))
			#cat("particle to select = ", particleToSelect, "\n")
			#cat("dput(oldParticleVector)\n")
			#dput(oldParticleVector)
			#cat("dput(oldParticleVector[particleToSelect])\n")
			#dput(oldParticleVector[particleToSelect])
			#cat("dput(oldParticleVector[[particleToSelect]])\n")
			#dput(oldParticleVector[[particleToSelect]])
			newparticleVector<-oldParticleVector[particleToSelect]
			#cat("dput(newparticleVector[[1]])\n")
			#dput(newparticleVector[[1]])
			#cat("mutateStates\n")
			newparticleVector[[1]]<-mutateStates(newparticleVector[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor)
			#cat("dput(newparticleVector[[1]]) AFTER MUTATE STATES\n")
			#dput(newparticleVector[[1]])
			newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=T))[1])
			if (plot) {
				plotcol="grey"
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				plotcol="black"
				if (dataGenerationStep==length(toleranceVector)) 
				{
					plotcol="red"
					}
				}
				text(x=intrinsicValues(newparticleVector[[1]]), y=distance(newparticleVector[[1]]), labels= dataGenerationStep, col=plotcol) 
			}
			#dput(convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleVector[[1]]@startingStates, newparticleVector[[1]]@intrinsicValues, newparticleVector[[1]]@extrinsicValues, timeStep), phy))
			#cat("dput(newparticleVector[[1]]) AFTER computeABCDistance\n")
			#dput(newparticleVector[[1]])
			
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				newparticleVector[[1]]<-setId(newparticleVector[[1]], particle)
				particle<-particle+1
				particleVector<-append(particleVector, newparticleVector)
				#now get weights, using correction in sisson et al. 2007
				newWeight=0
				for (i in 1:length(oldParticleVector)) {
					lnTransitionProb=log(1)
					for (j in 1:length(newparticleVector[[1]]@startingStates)) {
						newvalue<-newparticleVector[[1]]@startingStates[j]
						meantouse= oldParticleVector[[i]]@startingStates[j]
						sdtouse=standardDevFactor*(max(startingPriorsValues[, j])-min(startingPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(startingPriorsValues[, j])==max(startingPriorsValues[, j])) {
							localTransitionProb=1
						} 
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@intrinsicValues)) {
						newvalue<-newparticleVector[[1]]@intrinsicValues[j]
						meantouse= oldParticleVector[[i]]@intrinsicValues[j]
						sdtouse=standardDevFactor*(max(intrinsicPriorsValues[, j])-min(intrinsicPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(intrinsicPriorsValues[, j])==max(intrinsicPriorsValues[, j])) {
							localTransitionProb=1
						} 						
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@extrinsicValues)) {
						newvalue<-newparticleVector[[1]]@extrinsicValues[j]
						meantouse= oldParticleVector[[i]]@extrinsicValues[j]
						sdtouse=standardDevFactor*(max(extrinsicPriorsValues[, j])-min(extrinsicPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(extrinsicPriorsValues[, j])==max(extrinsicPriorsValues[, j])) {
							localTransitionProb=1
						} 						
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}					
					newWeight<-newWeight+getWeight(oldParticleVector[[i]])*exp(lnTransitionProb)
				
				
				}
				newparticleVector[[1]]<-setWeight(newparticleVector[[1]], newWeight)
				particleWeights[particle-1]<-newWeight
				weightScaling<-weightScaling+getWeight(newparticleVector[[1]])
		}
		else {
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
			while(sink.number()>0) {sink()}
			#print(newparticleVector)
			vectorForDataFrame<-c(dataGenerationStep, attempts, getId(newparticleVector[[1]]), particleToSelect, distance(newparticleVector[[1]]), getWeight(newparticleVector[[1]]), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]))
#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingStates = ", length(startingStates), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", distance(newparticleVector[[1]]), "\nweight = ", getWeight(newparticleVector[[1]]), "\n", vectorForDataFrame, "\n")
			particleDataFrame<-rbind(particleDataFrame, data.frame(rbind(vectorForDataFrame))) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
			cat(particle-1, attempts, floor(numParticles*attempts/particle), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), distance(newparticleVector[[1]]), "\n")
			
		}
		particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))
	}
	names(particleDataFrame)<-nameVector
	if(plot) {
		quartz()
		plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 1), type="n")
		for (i in 1:(length(toleranceVector)-1)) {
			graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
			lines(density(subset(particleDataFrame, generation==i)[, 8]), col= graycolor)
		}
		lines(density(subset(particleDataFrame, generation==length(toleranceVector))[, 8]), col= "red")
	}
	return(particleDataFrame)

	
	
	#------------------ ABC-PRC (end) ------------------
	
	

}


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

box.cox.powers.hacked<-function(X, start=NULL, hypotheses=NULL, ...){
	#modified from the car package to change the optim fn. The rest of the code is from J. Fox, author of car.
    modified.power<-function(x, lambda, gm){
        if (lambda == 0) log(x)*gm
        else (gm^(1-lambda))*((x^lambda)-1)/lambda
        }
    neg.kernel.profile.logL<-function(X, lambda, gm){
        for (j in 1:ncol(X)){
            X[, j]<-modified.power(X[, j], lambda[j], gm[j])
            }
        (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
        }
    univ.neg.kernel.logL <- function(x, lambda, gm){
        x <- modified.power(x, lambda, gm)
        (length(x)/2)*log(((length(x)-1)/length(x))*var(x))
        }
    X<-as.matrix(X)
    nc <- ncol(X)
    if(any(X<=0)) stop("All values must be > 0")
    gm<-apply(X, 2, function(x) exp(mean(log(x))))
    if (is.null(start)) {
        start <- rep(1, nc)
        for (j in 1:nc){
            res<- optimize(
                f = function(lambda) univ.neg.kernel.logL(x=X[, j], lambda=lambda, gm=gm[j]), 
                lower=-50, upper=+50)
            start[j] <- res$minimum
            }
        }
    res<-optim(start, neg.kernel.profile.logL, hessian=TRUE, method="BFGS", X=X, gm=gm, ...)
    result<-list()
    result$start<-start
    result$criterion<-res$value
    result$names<-colnames(X)
    result$lambda<-res$par
    result$stderr<-sqrt(diag(inv(res$hessian)))
    result$LR0<-2*(neg.kernel.profile.logL(X, rep(0, nc), gm)-res$value)
    result$LR1<-2*(neg.kernel.profile.logL(X, rep(1, nc), gm)-res$value)
    if (!is.null(hypotheses)) {
        for (i in 1:length(hypotheses)){
            if (length(hypotheses[[i]]) != nc) 
                stop(paste("hypothesis", i, "that powers =", hypotheses[[i]], "does not have", nc, "values"))
            hypotheses[[i]] <- list(test=2*(neg.kernel.profile.logL(X, hypotheses[[i]], gm)-res$value), 
                hypothesis=hypotheses[[i]])
            }
        result$hypotheses <- hypotheses
        }
    result$return.code<-res$convergence
    if(result$return.code != 0) 
        warning(paste("Convergence failure: return code =", 
            result$return.code))
    class(result)<-"box.cox.powers"
    result
    }

doSimulationForPlotting<-function(splits, intrinsicFn, extrinsicFn, startingStates, intrinsicValues, extrinsicValues, timeStep) {
	startVector<-c()
	endVector<-c()
	startTime<-c()
	endTime<-c()
	numberofsteps<-floor(splits[1, 1]/timeStep)
	mininterval<-min(splits[1:(dim(splits)[1]-1), 1]-splits[2:(dim(splits)[1]), 1])
	if (numberofsteps<1000) {
#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", splits[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
#warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more. I would suggest decreasing timeStep to no more than ", mininterval/50))
	}
#initial setup
	timefrompresent=splits[1, 1]
	taxa<-c(new("abctaxon", id=splits[1, 3], states=startingStates), new("abctaxon", id=splits[1, 4], states=startingStates))
	splits<-splits[2:dim(splits)[1], ] #pop off top value
	
#start running
	while(timefrompresent>0) {
#print(timefrompresent)
#speciation if needed
		while ((timefrompresent-timeStep)<=splits[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
			originallength<-length(taxa)
			taxontodelete<-Inf
			for (i in 1:originallength) {
				if (id(taxa[[i]])==splits[1, 2]) {
					taxontodelete<-i
					taxa<-c(taxa, taxa[[i]], taxa[[i]])
					id(taxa[[originallength+1]])<-splits[1, 3]
					timeSinceSpeciation(taxa[[originallength+1]])<-0
					id(taxa[[originallength+2]])<-splits[1, 4]
					timeSinceSpeciation(taxa[[originallength+2]])<-0
				}
			}
#cat("taxontodelete = ", taxontodelete)
			taxa<-taxa[-1*taxontodelete]
			if(dim(splits)[1]>1) {
				splits<-splits[2:(dim(splits)[1]), ] #pop off top value
			}
			else {
				splits[1, ]<-c(-1, 0, 0, 0)
			}
#print("------------------- speciation -------------------")
#print(taxa)
#summarizeTaxonStates(taxa)
		}
#trait evolution step
		for (i in 1:length(taxa)) {
			otherstatesvector<-c()
			for (j in 1:length(taxa)) {
				if (j!=i) {
					otherstatesvector<-c(otherstatesvector, states(taxa[[j]]))
				}
			}
#print(taxa)
#print(length(otherstatesvector))
			otherstatesmatrix<-matrix(otherstatesvector, ncol=length(states(taxa[[i]])), byrow=T) #each row represents one taxon
			newvalues<-states(taxa[[i]])+intrinsicFn(params=intrinsicValues, states=states(taxa[[i]]), timefrompresent =timefrompresent)+extrinsicFn(params=extrinsicValues, selfstates=states(taxa[[i]]), otherstates=otherstatesmatrix, timefrompresent =timefrompresent)
			nextstates(taxa[[i]])<-newvalues
			startVector<-append(startVector, states(taxa[[i]]))
			endVector <-append(endVector, newvalues)
			startTime <-append(startTime, timefrompresent+timeStep)
			endTime <-append(endTime, timefrompresent)
			
			
		}
		for (i in 1:length(taxa)) {

			states(taxa[[i]])<-nextstates(taxa[[i]])

			
		}
#print("------------------- step -------------------")
#print(taxa)
#summarizeTaxonStates(taxa)
		
		timefrompresent<-timefrompresent-timeStep
		for (i in 1:length(taxa)) {
			timeSinceSpeciation(taxa[[i]])<-timeSinceSpeciation(taxa[[i]])+timeStep
		}
	}
	plot(x=c(min(c(startVector, endVector)), max(c(startVector, endVector))), y=c(0, max(c(startTime, endTime))), type="n", ylab="Time", xlab="Trait value", main="", bty="n")
	for (i in 1:length(startVector)) {
		lines(x=c(startVector[i], endVector[i]), y=max(c(startTime, endTime)) - c(startTime[i], endTime[i]))	
	}
	return(data.frame(startVector, endVector, startTime, endTime))
}