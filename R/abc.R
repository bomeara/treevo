
setClass(
Class="abctaxon",
representation=representation(
id="numeric",
name="character",
states="numeric",
nextstates="numeric",
timeSinceSpeciation="numeric"
),
prototype(id=0,name="taxon0",timeSinceSpeciation=0)
)

setGeneric("id", function(object) {
	standardGeneric("id")
	})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
	})


setMethod("id","abctaxon",function(object) {
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


setMethod("name","abctaxon",function(object) {
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


setMethod("states","abctaxon",function(object) {
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


setMethod("nextstates","abctaxon",function(object) {
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


setMethod("timeSinceSpeciation","abctaxon",function(object) {
	object@timeSinceSpeciation
	})

setReplaceMethod("timeSinceSpeciation", signature(object="abctaxon", value="numeric") ,
	function(object, value) {
		object@timeSinceSpeciation <- value
		object
	}
)

setGeneric("updateStates", function(object) {
	standardGeneric("updateStates")
	})

	
setMethod("updateStates", signature(object='abctaxon'),
	function(object) {
		object@states<-object@nextstates
		object
		})

getSimulationSplits<-function(phy) {
	branchingTimes<-sort(branching.times(phy),decreasing=T)
	branchingTimesNames<-names(branchingTimes)
	ancestorVector<-c()
	descendant1Vector<-c()
	descendant2Vector<-c()
	for (i in 1:length(branchingTimes)) {
		relationshipVector<-phy$edge[phy$edge[,1]==branchingTimesNames[i]]
		ancestorVector<-c(ancestorVector, relationshipVector[2])
		descendant1Vector<-c(descendant1Vector, relationshipVector[3])
		descendant2Vector<-c(descendant2Vector, relationshipVector[4])
	}
	
	return(data.frame(branchingTimes,ancestorVector,descendant1Vector,descendant2Vector))
}

doSimulation<-function(splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep) {
	numberofsteps<-floor(splits[1,1]/timeStep)
	mininterval<-min(splits[1:(dim(splits)[1]-1),1]-splits[2:(dim(splits)[1]),1])
	if (numberofsteps<1000) {
#warning(paste("You have only ",numberofsteps," but should probably have a lot more. I would suggest decreasing timeStep to no more than ",splits[1,1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
#warning(paste("You have only ",floor(mininterval/timeStep)," on the shortest interval but should probably have a lot more. I would suggest decreasing timeStep to no more than ",mininterval/50))
	}
#initial setup
	timefrompresent=splits[1,1]
	taxa<-c(new("abctaxon",id=splits[1,3],states=startingStates),new("abctaxon",id=splits[1,4],states=startingStates))
	splits<-splits[2:dim(splits)[1],] #pop off top value
	
#start running
	while(timefrompresent>0) {
#print(timefrompresent)
#speciation if needed
		if ((timefrompresent-timeStep)<=splits[1,1]) { #do speciation
			originallength<-length(taxa)
			taxontodelete<-Inf
			for (i in 1:originallength) {
				if (id(taxa[[i]])==splits[1,2]) {
					taxontodelete<-i
					taxa<-c(taxa,taxa[[i]],taxa[[i]])
					id(taxa[[originallength+1]])<-splits[1,3]
					timeSinceSpeciation(taxa[[originallength+1]])<-0
					id(taxa[[originallength+2]])<-splits[1,4]
					timeSinceSpeciation(taxa[[originallength+2]])<-0
				}
			}
			#cat("taxontodelete = ",taxontodelete)
			taxa<-taxa[-1*taxontodelete]
			if(dim(splits)[1]>1) {
				splits<-splits[2:(dim(splits)[1]),] #pop off top value
			}
			else {
				splits[1,]<-c(-1,0,0,0)
			}
			#print("------------------- speciation -------------------")
#print(taxa)
			#summarizeTaxonStates(taxa)
		}
#trait evolution step
		for (i in 1:length(taxa)) {
			otherstatesvector<-c()))
			for (j in 1:length(taxa)) {
				if (j!=i) {
					otherstatesvector<-c(otherstatesvector,states(taxa[[j]]))
				}
			}
			otherstatesmatrix<-matrix(otherstatesvector,ncol=length(states(taxa[[i]])),byrow=T) #each row represents one taxon
			newvalues<-intrinsicFn(params=intrinsicValues,states=states(taxa[[i]]))+extrinsicFn(params=extrinsicValues,selfstates=states(taxa[[i]]),otherstates=otherstatesmatrix)
			nextstates(taxa[[i]])<-newvalues
		}
		for (i in 1:length(taxa)) {
			updateStates(taxa[[i]])
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
		#print("finished ",i)
	}
	statesmatrix<-matrix(statesvector,ncol=length(states(taxa[[1]])),byrow=T) #each row represents one taxon
	taxonframe<-data.frame(taxonid,taxonname,taxontimesincespeciation,statesmatrix)
	return(taxonframe)
}

brownianIntrinsic<-function(params,states) {
	newstates<-rnorm(n=length(states),mean=states,sd=params) #note use of standard deviation
#print("newstates =")
#	print(newstates)
	return(newstates)
}

brownianExtrinsic<-function(params,selfstates,otherstates) {
	newstates<-0*selfstates
#	print("extrinsic=")
#	print(newstates)
	return(newstates)
}

computeDistance<-function(simulationOutput,originalValues) {
	#both must have taxa on rows and chars in cols
	simulatedTraits<-unlist(simulationOutput[,4:dim(simulationOutput)[2]])
	names(simulatedTraits)<-NULL
	originalTraits<-unlist(originalValues)
	names(originalTraits)<-NULL
	euclideandistance<-dist(matrix(c(simulatedTraits, originalTraits),nrow=2,byrow=T))[1]
	return(euclideandistance)
	}
	
plotdistance1D<-function(xmin=0,xmax=1,numpoints=10,numreplicates=50) {
	xvalues<-seq(from=xmin,to=xmax,length.out=numpoints)
	yvalues<-seq(from=0,to=0,length.out=numpoints)
	for (i in 1:numpoints) { 
yvalues[i]<-median(replicate(numreplicates,computeDistance(doSimulation(splits=physplits,intrinsicFn=brownianIntrinsic,extrinsicFn=brownianExtrinsic,startingStates=c(xvalues[i]),intrinsicValues=0.5,extrinsicValues=0,timeStep=0.001),true)))
cat(xvalues[i]," ",yvalues[i])
		}
		plot(x=xvalues,y=yvalues)
	}
	
	convertTaxonFrameToGeigerData<-function(taxonframe,phy) {
		ntax<-dim(taxonframe)[1]
		newmatrix<-matrix(data=taxonframe[,4:(dim(taxonframe)[2])], nrow= ntax)
		newrownames<-c(rep(0, ntax))
		for (i in 1:ntax) {
			newrownames[i]<-phy$tip.label[(taxonframe$taxonid[i])]
			}
		geigerframe<-data.frame(newmatrix,row.names= newrownames, stringsAsFactors=F)
		geigerframe
		}
	
	geigerUnivariateSummaryStats<-function(phy, data) {
		#uses a bunch of stats from geiger. Only works for one character right now
		brownian<-fitContinuous(phy=phy,data= data,model="BM")
		white<-fitContinuous(phy=phy,data= data,model="white")
		lambda<-fitContinuous(phy=phy,data= data,model="lambda")
		kappa<-fitContinuous(phy=phy,data= data,model="kappa")
		delta<-fitContinuous(phy=phy,data= data,model="delta")
		EB<-fitContinuous(phy=phy,data= data,model="EB")
		summarystats<-c(brownian$Trait1$lnl, brownian$Trait1$beta, white$Trait1$lnl, white$Trait1$mean, lambda$Trait1$lambda, kappa$Trait1$lambda, delta$Trait1$delta, EB$Trait1$a)
		summarystats
		}
		
		