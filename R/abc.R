
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
									#cat("taxa[[i]]@id = ",taxa[[i]]@id, " splits[1,2] = ",splits[1,2],"\n")
				if (taxa[[i]]@id==splits[1,2]) {
					taxontodelete<-i
					taxa<-c(taxa,taxa[[i]],taxa[[i]])
					taxa[[originallength+1]]@id<-splits[1,3]
					taxa[[originallength+1]]@timeSinceSpeciation<-0
					taxa[[originallength+2]]@id<-splits[1,4]
					taxa[[originallength+2]]@timeSinceSpeciation<-0
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
			otherstatesvector<-c()
			for (j in 1:length(taxa)) {
				if (j!=i) {
					otherstatesvector<-c(otherstatesvector,taxa[[j]]@states)
				}
			}
			otherstatesmatrix<-matrix(otherstatesvector,ncol=length(taxa[[i]]@states),byrow=T) #each row represents one taxon
			newvalues<-intrinsicFn(params=intrinsicValues,states=taxa[[i]]@states)+extrinsicFn(params=extrinsicValues,selfstates=taxa[[i]]@states,otherstates=otherstatesmatrix)
			taxa[[i]]@nextstates<-newvalues
		}
		for (i in 1:length(taxa)) {
			taxa[[i]]@states<-taxa[[i]]@nextstates
		}
		#print("------------------- step -------------------")
#print(taxa)
		#summarizeTaxonStates(taxa)
		
		timefrompresent<-timefrompresent-timeStep
		for (i in 1:length(taxa)) {
			taxa[[i]]@timeSinceSpeciation<-taxa[[i]]@timeSinceSpeciation+timeStep
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
		statesvector<-c(statesvector, taxa[[i]]@states)
		#print("statesvector")
		taxonid<-c(taxonid, taxa[[i]]@id )
		#print("taxonid")
		taxonname<-c(taxonname, taxa[[i]]@name )
		#print("taxonname")
		taxontimesincespeciation<-c(taxontimesincespeciation, taxa[[i]]@timeSinceSpeciation)
		#print("finished ",i)
	}
	statesmatrix<-matrix(statesvector,ncol=length(taxa[[1]]@states),byrow=T) #each row represents one taxon
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
	euclideandistance<-dist(matrix(c(simulatedTraits, originalTraits),nrow=2))[1]
	return(euclideandistance)
	}