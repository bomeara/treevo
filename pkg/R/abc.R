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
			otherstatesvector<-c()
			for (j in 1:length(taxa)) {
				if (j!=i) {
					otherstatesvector<-c(otherstatesvector,states(taxa[[j]]))
				}
			}
			#print(taxa)
			#print(length(otherstatesvector))
			otherstatesmatrix<-matrix(otherstatesvector,ncol=length(states(taxa[[i]])),byrow=T) #each row represents one taxon
			newvalues<-intrinsicFn(params=intrinsicValues,states=states(taxa[[i]]))+extrinsicFn(params=extrinsicValues,selfstates=states(taxa[[i]]),otherstates=otherstatesmatrix)
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
	#sink(file="/dev/null") #because I really don't need output of which model I'm fitting. I already know.
#uses a bunch of stats from geiger. Only works for one character right now
	brownian<-fitContinuous(phy=phy,data= data,model="BM")
	white<-fitContinuous(phy=phy,data= data,model="white")
	lambda<-fitContinuous(phy=phy,data= data,model="lambda")
	kappa<-fitContinuous(phy=phy,data= data,model="kappa")
	delta<-fitContinuous(phy=phy,data= data,model="delta")
	EB<-fitContinuous(phy=phy,data= data,model="EB")
	summarystats<-c(brownian$Trait1$lnl, brownian$Trait1$beta, white$Trait1$lnl, white$Trait1$mean, lambda$Trait1$lambda, kappa$Trait1$lambda, delta$Trait1$delta, EB$Trait1$a)
	#sink() #turns output back on
	summarystats
}

geigerUnivariateSummaryStats2<-function(phy,data) {
	#sink(file="/dev/null") #because I really don't need output of which model I'm fitting. I already know.
#uses a bunch of stats from geiger. Only works for one character right now
	brownian<-fitContinuous(phy=phy,data= data,model="BM")
	white<-fitContinuous(phy=phy,data= data,model="white")
	lambda<-fitContinuous(phy=phy,data= data,model="lambda")
	kappa<-fitContinuous(phy=phy,data= data,model="kappa")
	delta<-fitContinuous(phy=phy,data= data,model="delta")
	EB<-fitContinuous(phy=phy,data= data,model="EB")
	summarystats<-c(brownian$Trait1$beta, white$Trait1$mean, lambda$Trait1$lambda, kappa$Trait1$lambda, delta$Trait1$delta, EB$Trait1$a)
	#sink() #turns output back on
	summarystats
}

abcprc2<-function(phy,originalData,intrinsicFn,extrinsicFn,summaryFn,startingMatrix,intrinsicMatrix,extrinsicMatrix,timeStep,toleranceVector,numParticles=10) {
#*matrix entries have lower (first row) and upper (second row) values for each parameter. Actual values are drawn from a uniform distribution for each parameter. If you want a parameter set, have max and min values the same number
#other things as in Sisson et al. Sequential Monte Carlo without likelihoods. P Natl Acad Sci Usa (2007) vol. 104 (6) pp. 1760-1765
#toleranceVector=vector of epsilons
#this implementation assumes flat priors
	nameVector<-c("generation","attempt","id","parentid","distance","weight")
	for (i in 1:dim(startingMatrix)[2]) {
		nameVector<-append(nameVector,paste("StartingStates",i,sep=""))
	}
	for (i in 1:dim(intrinsicMatrix)[2]) {
		nameVector<-append(nameVector,paste("IntrinsicValue",i,sep=""))
	}
	for (i in 1:dim(extrinsicMatrix)[2]) {
		nameVector<-append(nameVector,paste("ExtrinsicValue",i,sep=""))
	}
	#for (i in 1:Ntip(phy)) {
	#	nameVector<-append(nameVector,paste("Tip_",i,sep=""))
	#}
	
	splits<-getSimulationSplits(phy)
	originalSummary <-summaryFn(phy,originalData)
	particleWeights=rep(0,numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
	particleParameters<-matrix(nrow=numParticles,ncol=dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]) #stores parameters in model for each particle
	particleDistance=rep(NA,numParticles)
	particle<-1
	attempts<-0
	sink()
	cat("successes","attempts","expected number of attempts required")
	while (particle<=numParticles) {
		attempts<-attempts+1

		newparticleVector<-c(new("abcparticle",id=particle,generation=1))
		#print("new particle 1 = ")
		#print(newparticleVector[[1]])
		newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]],startingMatrix, intrinsicMatrix, extrinsicMatrix)
		#print("new particle 2 = ")
		#print(newparticleVector[[1]])
		#print("new particle 3 = ")
		#print(newparticleVector[[1]])
		newparticleVector[[1]]<-computeABCDistance(newparticleVector[[1]], summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep)
		#print("new particle 4 = ")
		#print(newparticleVector[[1]])
		if (distance(newparticleVector[[1]])<toleranceVector[1]) {
			id(newparticleVector[[1]])<-particle
			weight(newparticleVector[[1]])<-1/numParticles
			particleWeights[particle]<-1/numParticles
			particle<-particle+1
		}
		vectorForDataFrame<-c(1,attempts,id(newparticleVector[[1]]),0,distance(newparticleVector[[1]]),weight(newparticleVector[[1]]),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]))
		
		STICK THE ABOVE IN A DATA FRAME USING THE NAMEVECTOR. NOTE THAT PARTICLES NOT SELECTED WILL HAVE ID==0
		
		
		
		cat(particle-1,attempts,floor(numParticles*attempts/particle),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),distance(newparticleVector[[1]]),"\n")

	}

	for (dataGenerationStep in 2:length(toleranceVector)) {
	}
}


abcprc<-function(phy,originalData,intrinsicFn,extrinsicFn,summaryFn,startingMatrix,intrinsicMatrix,extrinsicMatrix,timeStep,toleranceVector,numParticles=50) {
#*matrix entries have lower (first row) and upper (second row) values for each parameter. Actual values are drawn from a uniform distribution for each parameter. If you want a parameter set, have max and min values the same number
#other things as in Sisson et al. Sequential Monte Carlo without likelihoods. P Natl Acad Sci Usa (2007) vol. 104 (6) pp. 1760-1765
#toleranceVector=vector of epsilons
#this implementation assumes flat priors
	splits<-getSimulationSplits(phy)
	originalSummaryStats<-summaryFn(phy,originalData)
	
	particleWeights=rep(0,numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
	particleParameters<-matrix(nrow=numParticles,ncol=dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]) #stores parameters in model for each particle
	particleDistance=rep(NA,numParticles)
	particle<-1
	attempts<-0
	cat("successes","attempts","expected number of attempts required")
	while (particle<=numParticles) {
		attempts<-attempts+1)))
		startingStates<-vector(mode="numeric",length=dim(startingMatrix)[2])
		intrinsicValues<-vector(mode="numeric",length=dim(intrinsicMatrix)[2])
		extrinsicValues <-vector(mode="numeric",length=dim(extrinsicMatrix)[2])
		#print(intrinsicValues)
		#print(extrinsicValues)
		placementCounter=1;
		for (j in 1:length(startingStates)) {
			startingStates[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
			particleParameters[particle,placementCounter]<-startingStates[j]
			placementCounter<-placementCounter+1
		}
		for (j in 1:length(intrinsicValues)) {
			intrinsicValues[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
			particleParameters[particle,placementCounter]<-intrinsicValues[j]
			placementCounter<-placementCounter+1
			
		}
		for (j in 1:length(extrinsicValues)) {
			extrinsicValues[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
			particleParameters[particle,placementCounter]<-extrinsicValues[j]
			placementCounter<-placementCounter+1
		}
		particleDistance[particle]<-dist(matrix(c(summaryFn(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep),phy)), originalSummaryStats),nrow=2,byrow=T))[1]
		if (particleDistance[particle]<toleranceVector[1]) {
			particleWeights[particle]<-1/numParticles
			particle<-particle+1
		}
		cat(particle-1,attempts,floor(numParticles*attempts/particle))
	}

	for (dataGenerationStep in 2:length(toleranceVector)) {
	}
}

generateRandomSampleFromMatrices<-function(startingMatrix,intrinsicMatrix,extrinsicMatrix) {
		startingStates<-vector(mode="numeric",length=dim(startingMatrix)[2])
		intrinsicValues<-vector(mode="numeric",length=dim(intrinsicMatrix)[2])
		extrinsicValues <-vector(mode="numeric",length=dim(extrinsicMatrix)[2])
		#startingStates=rep(NA,dim(startingMatrix)[2])
		#intrinsicValues=rep(NA,dim(intrinsicMatrix)[2])
		#extrinsicValues=rep(NA,dim(extrinsicMatrix)[2])
		for (j in 1:length(startingStates)) {
			startingStates[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
		}
		for (j in 1:length(intrinsicValues)) {
			#print(intrinsicValues)
			intrinsicValues[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
			
		}
		for (j in 1:length(extrinsicValues)) {
			#print(extrinsicValues)
			extrinsicValues[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
		}

	
	}






profileAcrossUniform<-function(phy,originalData,intrinsicFn,extrinsicFn,startingMatrix,intrinsicMatrix,extrinsicMatrix,timeStep,numreps=100) {
#*matrix entries have lower (first row) and upper (second row) values for each parameter. Actual values are drawn from a uniform distribution for each parameter. If you want a parameter set, have max and min values the same number
	print("cow")
	splits<-getSimulationSplits(phy)
	print(splits)
	originalSummaryStats<-geigerUnivariateSummaryStats(phy,originalData)
	print(originalSummaryStats)
	resultsMatrix<-matrix(nrow=numreps,ncol=length(originalSummaryStats) + dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]+1)
	namesVector<-rep(NA,length(originalSummaryStats) + dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]+1)
	for (i in 1:numreps) {
		startingStates<-vector(mode="numeric",length=dim(startingMatrix)[2])
		intrinsicValues<-vector(mode="numeric",length=dim(intrinsicMatrix)[2])
		extrinsicValues <-vector(mode="numeric",length=dim(extrinsicMatrix)[2])
		placementCounter=1;
		for (j in 1:length(startingStates)) {
			startingStates[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
			resultsMatrix[i,placementCounter]<-startingStates[j]
			if (i==1) {
				namesVector[placementCounter]<-paste("startingValue",j,sep="")
			}
			placementCounter<-placementCounter+1
		}
		for (j in 1:length(intrinsicValues)) {
			intrinsicValues[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
			resultsMatrix[i,placementCounter]<-intrinsicValues[j]
			if (i==1) {
				namesVector[placementCounter]<-paste("intrinsicValue",j,sep="")
			}
			placementCounter<-placementCounter+1
			
		}
		for (j in 1:length(extrinsicValues)) {
			extrinsicValues[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
			resultsMatrix[i,placementCounter]<-extrinsicValues[j]
			if (i==1) {
				namesVector[placementCounter]<-paste("extrinsicValue",j,sep="")
			}
			placementCounter<-placementCounter+1
		}
		individualResult<-geigerUnivariateSummaryStats(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep),phy))
		for (j in 1:length(individualResult)) {
			resultsMatrix[i,placementCounter]<-individualResult[j]
			if (i==1) {
				namesVector[placementCounter]<-paste("SummaryStat",j,sep="")
			}
			placementCounter<-placementCounter+1
		}
		resultsMatrix[i,placementCounter]<-dist(matrix(c(individualResult, originalSummaryStats),nrow=2,byrow=T))[1]
		if (i==1) {
			namesVector[placementCounter]<-"distance"
		}
		placementCounter<-placementCounter+1
		print(resultsMatrix[i,])
	}
	resultsDataFrame<-data.frame(resultsMatrix)
	names(resultsDataFrame)<-namesVector
	return(resultsDataFrame)
}


#test code
library(geiger)
phy<-rcoal(9)
char<-data.frame(5+sim.char(phy,model.matrix=matrix(10),1))
Rprof()
profile<-profileAcrossUniform(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001,numreps=10)
Rprof(NULL)
summaryRprof()
profile
plot(profile$startingValue1,profile$distance)

#test code2
library(geiger)
phy<-rcoal(9)
char<-data.frame(5+sim.char(phy,model.matrix=matrix(10),1))
abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(50,2),summaryFn= geigerUnivariateSummaryStats2)
