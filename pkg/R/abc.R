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
	sink(file="/dev/null") #because I really don't need output of which model I'm fitting. I already know.
#uses a bunch of stats from geiger. Only works for one character right now
	brownian<-fitContinuous(phy=phy,data= data,model="BM")
	white<-fitContinuous(phy=phy,data= data,model="white")
	lambda<-fitContinuous(phy=phy,data= data,model="lambda")
	kappa<-fitContinuous(phy=phy,data= data,model="kappa")
	delta<-fitContinuous(phy=phy,data= data,model="delta")
	EB<-fitContinuous(phy=phy,data= data,model="EB")
	summarystats<-c(brownian$Trait1$beta, white$Trait1$mean, lambda$Trait1$lambda, kappa$Trait1$lambda, delta$Trait1$delta, EB$Trait1$a)
	sink() #turns output back on
	summarystats
}

rawValuesSummaryStats<-function(phy,data) {
	nrow=dim(data[1])[1]
	summarystats<-rep(NA,nrow)
	for (i in 1:nrow) {
		summarystats[i]<-data[i,1]
	}
	summarystats
}

abcprc2<-function(phy,originalData,intrinsicFn,extrinsicFn,summaryFn,startingMatrix,intrinsicMatrix,extrinsicMatrix,timeStep,toleranceVector,numParticles=10, standardDevFactor=0.05, plot=F) {
#*matrix entries have lower (first row) and upper (second row) values for each parameter. Actual values are drawn from a uniform distribution for each parameter. If you want a parameter set, have max and min values the same number
#other things as in Sisson et al. Sequential Monte Carlo without likelihoods. P Natl Acad Sci Usa (2007) vol. 104 (6) pp. 1760-1765
#toleranceVector=vector of epsilons
#this implementation assumes flat priors
if (plot) {
	plot(x=c(min(intrinsicMatrix),max(intrinsicMatrix)),y=c(0,5*max(toleranceVector)),type="n")
	}
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
	particleDataFrame<-data.frame()
	cat("successes","attempts","expected number of attempts required\n")
	particleVector<-c()
	while (particle<=numParticles) {
		attempts<-attempts+1
		
		newparticleVector<-c(new("abcparticle",id=particle,generation=1,weight=0))
		newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]],startingMatrix, intrinsicMatrix, extrinsicMatrix)
		newparticleVector[[1]]<-computeABCDistance(newparticleVector[[1]], summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep)
		if (distance(newparticleVector[[1]])<toleranceVector[1]) {
			newparticleVector[[1]]<-setId(newparticleVector[[1]],particle)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]],1/numParticles)
			particleWeights[particle]<-1/numParticles
			particle<-particle+1
			particleVector<-append(particleVector,newparticleVector)
		}
		if (plot) {
			plotcol="grey"
			if (distance(newparticleVector[[1]])<toleranceVector[1]) {
				plotcol="black"
				}
			text(x=intrinsicValues(newparticleVector[[1]]),y=distance(newparticleVector[[1]]),labels=1,col=plotcol) 
			}
		sink()
		#print(newparticleVector)
		vectorForDataFrame<-c(1,attempts,getId(newparticleVector[[1]]),0,distance(newparticleVector[[1]]),getWeight(newparticleVector[[1]]),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]))
#cat("\n\nlength of vectorForDataFrame = ",length(vectorForDataFrame),"\n","length of startingStates = ",length(startingStates),"\nlength of intrinsicValues = ",length(intrinsicValues),"\nlength of extrinsicValues = ",length(extrinsicValues),"\ndistance = ",distance(newparticleVector[[1]]),"\nweight = ",getWeight(newparticleVector[[1]]),"\n",vectorForDataFrame,"\n")
		particleDataFrame<-rbind(particleDataFrame,data.frame(rbind(vectorForDataFrame)))
		cat(particle-1,attempts,floor(numParticles*attempts/particle),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),distance(newparticleVector[[1]]),"\n")
		
	}
	
	for (dataGenerationStep in 2:length(toleranceVector)) {
		cat("\n\n\n", "STARTING DATA GENERATION STEP ", dataGenerationStep,  "\n\n\n")
		particleWeights<-particleWeights/(sum(particleWeights)) #normalize to one
		cat("particleWeights\n", particleWeights,"\n\n")
		oldParticleVector<-particleVector
		oldParticleWeights<-particleWeights
		particleWeights=rep(0,numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
		particleParameters<-matrix(nrow=numParticles,ncol=dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]) #stores parameters in model for each particle
		particleDistance=rep(NA,numParticles)
		particle<-1
		attempts<-0
		cat("successes","attempts","expected number of attempts required\n")
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
			newparticleVector[[1]]<-mutateStates(newparticleVector[[1]],startingMatrix, intrinsicMatrix, extrinsicMatrix, standardDevFactor)
			#cat("dput(newparticleVector[[1]]) AFTER MUTATE STATES\n")
			#dput(newparticleVector[[1]])
			newparticleVector[[1]]<-computeABCDistance(newparticleVector[[1]], summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep)
			if (plot) {
				plotcol="grey"
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				plotcol="black"
				if (dataGenerationStep==length(toleranceVector)) 
				{
					plotcol="red"
					}
				}
				text(x=intrinsicValues(newparticleVector[[1]]),y=distance(newparticleVector[[1]]),labels= dataGenerationStep,col=plotcol) 
			}
			#dput(convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,newparticleVector[[1]]@startingStates,newparticleVector[[1]]@intrinsicValues,newparticleVector[[1]]@extrinsicValues,timeStep),phy))
			#cat("dput(newparticleVector[[1]]) AFTER computeABCDistance\n")
			#dput(newparticleVector[[1]])
			
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				newparticleVector[[1]]<-setId(newparticleVector[[1]],particle)
				particle<-particle+1
				particleVector<-append(particleVector,newparticleVector)
				#now get weights, using correction in sisson et al. 2007
				newWeight=0
				for (i in 1:length(oldParticleVector)) {
					lnTransitionProb=log(1)
					for (j in 1:length(newparticleVector[[1]]@startingStates)) {
						newvalue<-newparticleVector[[1]]@startingStates[j]
						meantouse= oldParticleVector[[i]]@startingStates[j]
						sdtouse=standardDevFactor*(max(startingMatrix[,j])-min(startingMatrix[,j]))
						localTransitionProb=dnorm(newvalue,mean= meantouse,sd= sdtouse)/(1-pnorm(min(startingMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=T)+pnorm(max(startingMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=F))
						if (min(startingMatrix[,j])==max(startingMatrix[,j])) {
							localTransitionProb=1
						} 
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@intrinsicValues)) {
						newvalue<-newparticleVector[[1]]@intrinsicValues[j]
						meantouse= oldParticleVector[[i]]@intrinsicValues[j]
						sdtouse=standardDevFactor*(max(intrinsicMatrix[,j])-min(intrinsicMatrix[,j]))
						localTransitionProb=dnorm(newvalue,mean= meantouse,sd= sdtouse)/(1-pnorm(min(intrinsicMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=T)+pnorm(max(intrinsicMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=F))
						if (min(intrinsicMatrix[,j])==max(intrinsicMatrix[,j])) {
							localTransitionProb=1
						} 						
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@extrinsicValues)) {
						newvalue<-newparticleVector[[1]]@extrinsicValues[j]
						meantouse= oldParticleVector[[i]]@extrinsicValues[j]
						sdtouse=standardDevFactor*(max(extrinsicMatrix[,j])-min(extrinsicMatrix[,j]))
						localTransitionProb=dnorm(newvalue,mean= meantouse,sd= sdtouse)/(1-pnorm(min(extrinsicMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=T)+pnorm(max(extrinsicMatrix[,j]),mean= meantouse ,sd= sdtouse,lower.tail=F))
						if (min(extrinsicMatrix[,j])==max(extrinsicMatrix[,j])) {
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
			sink()
			#print(newparticleVector)
			vectorForDataFrame<-c(dataGenerationStep,attempts,getId(newparticleVector[[1]]), particleToSelect,distance(newparticleVector[[1]]),getWeight(newparticleVector[[1]]),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]))
#cat("\n\nlength of vectorForDataFrame = ",length(vectorForDataFrame),"\n","length of startingStates = ",length(startingStates),"\nlength of intrinsicValues = ",length(intrinsicValues),"\nlength of extrinsicValues = ",length(extrinsicValues),"\ndistance = ",distance(newparticleVector[[1]]),"\nweight = ",getWeight(newparticleVector[[1]]),"\n",vectorForDataFrame,"\n")
			particleDataFrame<-rbind(particleDataFrame,data.frame(rbind(vectorForDataFrame))) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
			cat(particle-1,attempts,floor(numParticles*attempts/particle),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),distance(newparticleVector[[1]]),"\n")
			
		}
		
		
		
		
		
	}
	names(particleDataFrame)<-nameVector
	if(plot) {
		quartz()
		plot(x=c(min(intrinsicMatrix),max(intrinsicMatrix)),y=c(0,1),type="n")
		for (i in 1:(length(toleranceVector)-1)) {
			graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
			lines(density(subset(particledata,generation==i)[,8]),col= graycolor)
		}
		lines(density(subset(particledata,generation==length(toleranceVector))[,8]),col= "red")
	}
	particleDataFrame
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
		sink()
		print(newparticleVector)
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


#test code2
library(geiger)
phy<-rcoal(9)
char<-data.frame(5+sim.char(phy,model.matrix=matrix(20),1))
#Rprof()
#particledata<-abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(500,2),summaryFn= geigerUnivariateSummaryStats2)

particledata<-abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(500,400,300, 200, 100),standardDevFactor=0.1, summaryFn= rawValuesSummaryStats,plot=T,numParticles=30)

#Rprof(NULL)
#summaryRprof()

#test code
#library(geiger)
#phy<-rcoal(9)
#char<-data.frame(5+sim.char(phy,model.matrix=matrix(10),1))
#Rprof()
#profile<-profileAcrossUniform(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15,0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001,numreps=10)
#Rprof(NULL)
#summaryRprof()
#profile
#plot(profile$startingValue1,profile$distance)


