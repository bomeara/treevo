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

#uses try to deal with rare errors (singular matrices, for example) in Geiger. Rather than aborting the entire run, just gives back an NA, and dist fn can deal with NAs (though this introduces a bias in creating small distances when something fails)
	summarystats<-as.numeric(
		c(
			try(fitContinuous(phy=phy,data= data,model="BM")$Trait1$beta,silent=T),  
			try(fitContinuous(phy=phy,data= data,model="white")$Trait1$mean,silent=T), 
			try(fitContinuous(phy=phy,data= data,model="lambda")$Trait1$lambda,silent=T), 
			try(fitContinuous(phy=phy,data= data,model="kappa")$Trait1$lambda,silent=T), 
			try(fitContinuous(phy=phy,data= data,model="delta")$Trait1$delta,silent=T), 
			try(fitContinuous(phy=phy,data= data,model="EB")$Trait1$a,silent=T)
		)
	)
	sink() #turns output back on
	sink() #do it again, as it seems to require more than once
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

summaryStatsLong<-function(phy,data,todo=c()) {
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
	raw.median<-as.numeric(try(median(data[todo[22],])))
	#cat("summaryStatsLong")
	summarystats<-c(brown.lnl, brown.beta, brown.aic,  lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic,  raw.mean, raw.max, raw.min, raw.var, raw.median, data[[1]] )
	#cat("\n summaryStatsLong summarystats1\n")
	#print(summarystats)
	summarystats[which(todo==0)]<-NA
	#cat("\n summaryStatsLong summarystats2\n")
	#print(summarystats)
	summarystats[which(is.finite(summarystats)==FALSE)]<-NA
	#cat("\n summaryStatsLong summarystats3\n")
	#print(summarystats)
	sink()
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

estimateMeanAndSigma<-function(splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep,nrep=10,phy,returnTraits=F) {
	#calculates mean and sample covariance
	simulatedValues<-replicate(n=nrep,expr=convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn= intrinsicFn,extrinsicFn= extrinsicFn,startingStates= startingStates,intrinsicValues= intrinsicValues,extrinsicValues= extrinsicValues,timeStep=timeStep),phy))
	simulatedMatrix<-as.matrix(as.data.frame(simulatedValues,row.names=phy$tip.label))
	ntax<-dim(simulatedMatrix)[1]
	means<-rep(NA, ntax)
	for (i in 1:ntax) {
		means[i]<-mean(simulatedMatrix[i,])
	}
	names(means)<-phy$tip.label
	sigma<-matrix(data=NA,nrow=ntax,ncol=ntax,dimnames=list(phy$tip.label, phy$tip.label))
	for (i in 1:ntax) {
		for (j in i:ntax) {
			if (i!=j) {
				sigma[i,j]<-cov(x=simulatedMatrix[i,],y=simulatedMatrix[j,])
				sigma[j,i]<-sigma[i,j]
			}	
			else {
				sigma[i,j]<-var(x=simulatedMatrix[i,])
			}
		}	
	}
	estimatedMeanAndSigma<-list(mean=means,sigma=sigma)
	if (returnTraits) {
		estimatedMeanAndSigma<-list(mean=means, sigma=sigma, traits=simulatedMatrix)
	}
	return(estimatedMeanAndSigma)
}

GeigerEstimateMeanAndSigma<-function(nrep=10,phy,returnTraits=F) {
	#calculates mean and sample covariance
	simulatedValues<-replicate(n=nrep,expr=data.frame(5+sim.char(phy,model.matrix=matrix(10),1)))
	simulatedMatrix<-as.matrix(as.data.frame(simulatedValues,row.names=phy$tip.label))
	ntax<-dim(simulatedMatrix)[1]
	means<-rep(NA, ntax)
	for (i in 1:ntax) {
		means[i]<-mean(simulatedMatrix[i,])
	}
	names(means)<-phy$tip.label
	sigma<-matrix(data=NA,nrow=ntax,ncol=ntax,dimnames=list(phy$tip.label, phy$tip.label))
	for (i in 1:ntax) {
		for (j in i:ntax) {
			if (i!=j) {
				sigma[i,j]<-cov(x=simulatedMatrix[i,],y=simulatedMatrix[j,])
				sigma[j,i]<-sigma[i,j]
			}	
			else {
				sigma[i,j]<-var(x=simulatedMatrix[i,])
			}
		}	
	}
	estimatedMeanAndSigma<-list(mean=means,sigma=sigma)
	if (returnTraits) {
		estimatedMeanAndSigma<-list(mean=means, sigma=sigma, traits=simulatedMatrix)
	}
	return(estimatedMeanAndSigma)
}

estimateLnLikelihoodDmvnorm <-function(x, mean,sigma) {
	if (!is.vector(x)) {
		ntax<-dim(x)[1]
		newvector<-rep(NA,ntax)
		for (i in 1:ntax) {
			newvector[i]<-x[i,1]	
		}	
		x<-newvector
	}
	if (!is.vector(mean)) {
		ntax<-dim(mean)[1]
		newvector<-rep(NA,ntax)
		for (i in 1:ntax) {
			newvector[i]<-mean[i,1]	
		}	
		mean <-newvector
	}
	if (is.vector(mean)) {
       mean <- matrix(mean, ncol = length(mean))
    }
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
 	return(dmvnorm(x=x,mean=mean,sigma=sigma, log=T))
}



estimateLnLikelihood<-function(x, mean,sigma,forcePositiveDet=F) {
	ntax<-dim(sigma)[1]
	
	#this rigamarole is to deal with weird input
	diff<-rep(NA,ntax)
		if (!is.vector(x)) {
		ntax<-dim(x)[1]
		newvector<-rep(NA,ntax)
		for (i in 1:ntax) {
			newvector[i]<-x[i,1]	
		}	
		x<-newvector
	}
	if (!is.vector(mean)) {
		ntax<-dim(mean)[1]
		newvector<-rep(NA,ntax)
		for (i in 1:ntax) {
			newvector[i]<-mean[i,1]	
		}	
		mean <-newvector
	}
	for (i in 1:ntax) {
		diff[i]<-x[i]-mean[i]	
	}
	
	#now to business
	#invV<-chol2inv(chol(sigma))
	invV<-ginv(sigma) #gotta use this as have problems with non-invertible matrices
	detV<-det(sigma)
	if(forcePositiveDet && detV<0) {
		detV<-abs(detV)	
		warningtext<-paste("detV = ",detV,", used in estimateLnLikelihood, is negative, which will make taking the squareroot difficult. abs(detV) taken instead.",sep="")
		warning(warningtext)
	}
	lnL<-((-0.5*(diff)%*%invV%*%(diff))-log(sqrt(((2*pi)^ntax)*detV)))
	as.vector(lnL)

}

simulateForLnL<-function(data,splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep,nrep=10,phy) {
	estimatedMeanAndSigma<-estimateMeanAndSigma(splits,intrinsicFn,extrinsicFn,startingStates,intrinsicValues,extrinsicValues,timeStep,nrep,phy,returnTraits=F)
	lnL<-estimateLnLikelihood(x=data,sigma=estimatedMeanAndSigma$sigma,mean=estimatedMeanAndSigma$mean,forcePositiveDet=T)
	lnL
}

testApproach<-function(phy,originalData,intrinsicFn,extrinsicFn,summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2),startingMatrix,intrinsicMatrix,extrinsicMatrix,startingStatesGuess,intrinsicValuesGuess,extrinsicValuesGuess,timeStep,toleranceVector,nrepSim=100, nrepEst=1000, standardDevFactor=0.05, plot=T) {
	#here, startingStatesGuess, intrinsicValuesGuess, and extrinsicValuesGuess are set to a hypothesized true value
		numberParametersTotal<-dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]
		numberParametersFree<-numberParametersTotal
		numberParametersStarting<-0
		numberParametersIntrinsic<-0
		numberParametersExtrinsic<-0
		freevariables<-matrix(data=NA, nrow=2,ncol=0)
		truthvector<-cbind(startingStatesGuess, intrinsicValuesGuess, extrinsicValuesGuess) #includes free and fixed parameters
		titlevector<-c()
		freevector<-c()
		
		for (i in 1:dim(startingMatrix)[2]) {
			if (startingMatrix[1,i]== startingMatrix[2,i]) {
				numberParametersFree<-numberParametersFree-1
				freevector<-c(freevector,FALSE)
			}	
			else {
				numberParametersStarting<-numberParametersStarting+1
				freevariables<-cbind(freevariables,startingMatrix[,i])
				titlevector <-c(titlevector,paste("Starting", numberParametersStarting))
				freevector<-c(freevector,TRUE)
			}
		}
		for (i in 1:dim(intrinsicMatrix)[2]) {
			if (intrinsicMatrix[1,i]== intrinsicMatrix[2,i]) {
				numberParametersFree<-numberParametersFree-1
				freevector<-c(freevector,FALSE)
			}	
			else {
				numberParametersIntrinsic <-numberParametersIntrinsic +1
				freevariables<-cbind(freevariables, intrinsicMatrix[,i])
				titlevector <-c(titlevector,paste("Intrinsic", numberParametersIntrinsic))
				freevector<-c(freevector,TRUE)
			}

		}
		for (i in 1:dim(extrinsicMatrix)[2]) {
			if (extrinsicMatrix[1,i]== extrinsicMatrix[2,i]) {
				numberParametersFree<-numberParametersFree-1
				freevector<-c(freevector,FALSE)
			}	
			else {
				numberParametersExtrinsic <-numberParametersExtrinsic +1
				freevariables<-cbind(freevariables, extrinsicMatrix[,i])
				titlevector <-c(titlevector,paste("Extrinsic", numberParametersExtrinsic))
				freevector<-c(freevector,TRUE)
			}

		}
		truthvector<-truthvector[freevector] #so it now includes free params only
		#numberValuesPerParameter<-floor(nrep^(1/numberParametersFree))
		#cat("\n\nThis will try ", numberValuesPerParameter," unique values for each of ", numberParametersFree, " free parameters in a grid\n\n")
		particleDataFrame<-data.frame()
		splits<-getSimulationSplits(phy)
		simulatedValues<-as.data.frame(replicate(n=nrepSim+1,expr=convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn= intrinsicFn,extrinsicFn= extrinsicFn,startingStates= startingStatesGuess,intrinsicValues= intrinsicValuesGuess,extrinsicValues= extrinsicValuesGuess,timeStep=timeStep),phy)),row.names=phy$tip.label)
		#print(simulatedValues)
		guessSummaries<-matrix(nrow=nrepSim,ncol=2)
		toOrigSummaries<-matrix(nrow=nrepSim,ncol=2)
		#print(guessSummaries)
		for (i in 1:length(summaryFns)) {
			for (j in 1:nrepSim) {
				guessSummaries[j,i]<-dist(matrix(c(summaryFns[i][[1]](phy, simulatedValues[j]), summaryFns[i][[1]](phy, simulatedValues[j+1])),nrow=2,byrow=T))[1]
				toOrigSummaries[j,i]<-dist(matrix(c(summaryFns[i][[1]](phy, simulatedValues[j]), summaryFns[i][[1]](phy, originalData)),nrow=2,byrow=T))[1]

			}	
			#print(guessSummaries)
		}
		quartz()
		par(mfcol=c(1,length(summaryFns)))
		for(i in 1:length(summaryFns)) {
			dguess<-density(guessSummaries[,i])
			dorig<-density(toOrigSummaries[,i])
			plot(x=c(min(dguess$x,dorig$x),max(dguess$x,dorig$x)),y=c(min(dguess$y,dorig$y),max(dguess$y,dorig$y)),type="n")
			lines(dguess)
			lines(dorig,col="red")
			
			#title(main=as.character(summmaryFns[i]))
		}
		quartz()
		par(mfrow=c(length(summaryFns),numberParametersFree))
		infomatrix<-matrix(nrow=0,ncol=3+numberParametersFree) #will store c(summaryindex, summaryX_newguessdistance, summaryX_neworiginaldistance, summaryX_freeparam1, summaryS_freeparam2,... summaryS_newguessdistance, ....)
		for (currentrep in 1:nrepEst) {
			cat(currentrep,"\n")
			newparticleVector<-c(new("abcparticle",id=1,generation=0,weight=0))
			newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]],startingMatrix, intrinsicMatrix, extrinsicMatrix)
			simulatedTips<-simulateTips(newparticleVector[[1]], splits, phy, intrinsicFn, extrinsicFn, timeStep)
			
			for(i in 1:length(summaryFns)) {
				newGuessDistance<-dist(matrix(c(summaryFns[i][[1]](phy, simulatedValues[1]), summaryFns[i][[1]](phy, simulatedTips)),nrow=2,byrow=T))[1]
				newOriginalDistance<-dist(matrix(c(summaryFns[i][[1]](phy, originalData), summaryFns[i][[1]](phy, simulatedTips)),nrow=2,byrow=T))[1]
				values<-c(startingStates (newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]))
				values<-values[freevector]
				infomatrix<-rbind(infomatrix,c(i, newGuessDistance, newOriginalDistance, values))
				#for (j in 1:numberParametersFree) {
				#	par(mfg=c(i,j))
				#	sink()
				#	sink()
				#	cat(i," ",j, " ",values[j],"\n")
				#	points(x=c(values[j]),y=c(newGuessDistance),pch=1,col=rgb(0,0,0,,0.7))
				#	points(x=c(values[j]),y=c(newOriginalDistance),pch=".",col="red")
				#}
			}
			if(currentrep>2) {
				sink()
				sink()
				sink()
				#cat("infomatrix\n")
				#print(infomatrix)
				for(i in 1:length(summaryFns)) {
					
					infomatrixA<-subset(infomatrix,infomatrix[,1]==i) #pulls out rows with correct summary fn
					#cat("infomatrixA\n")
					#print(infomatrixA)
					
					cutoffs<-c(0.05,0.1,.25,.50,1)
					for (j in 1:numberParametersFree) {
						#par(mfg=c(i,j))
						par(bg="white")
					plot(x=c(min(freevariables[,j]),max(freevariables[,j])),y=c(0,1.25*max(infomatrixA[,2:3])),type="n",main= paste("Summary fn",i,titlevector[j]),xlab="parameter value",ylab="distance",frame.plot=FALSE)

						
						densityheight=0.15*max(infomatrixA[,2:3])
						
						lines(x=c(truthvector[j],truthvector[j]),y=c(-0.2*max(infomatrixA[,2:3]),1.25*max(infomatrix[,2:3])),col="red")
						for (cutoffid in 1:length(cutoffs)) {
							cutoff<-cutoffs[cutoffid]
							if (1/cutoff < currentrep) {
								mindistance<-quantile(infomatrixA[,2],probs=cutoff)
								infomatrixB<-subset(infomatrixA,infomatrixA[,2]<mindistance)
								#cat("infomatrixB\n")
								#print(infomatrixB)
								if (dim(infomatrixB)[1]>=2) {
									densityvalues<-density(infomatrixB[,j+3],from=min(freevariables[,j]),to=max(freevariables[,j]))
									densityvalues$y<-densityvalues$y*densityheight/max(densityvalues$y) + mindistance
									
									lines(densityvalues)
									
									points(x=infomatrixB[,j+3] , y=rep(mindistance, dim(infomatrixB)[1]), pch=".", col=rgb(0,0,0,.5))
									
								}
							}
						}	
					}
				}	
			}
		}
}


boxcoxplsSummary<-function(todo,rawSummaryValues,prunedPlsResult, boxcoxLambda, boxcoxAddition) {
	#cat("\n\n boxcoxplsSummary \ntodo\n")
	#print(todo)
	#cat("\n rawSummaryValues\n")
	#print(rawSummaryValues)
	#cat("\n prunedPlsResult\n")
	#print(prunedPlsResult)
	#cat("\n boxcoxLambda\n")
	#print(boxcoxLambda)
	#cat("\n boxcoxAddition\n")
	#print(boxcoxAddition)
	summaryValues<-boxcoxAddition+rawSummaryValues #for boxcox, need positive numbers
	#cat("\n summaryValues\n")
	#print(summaryValues)

	for (summaryValueIndex in 1:length(summaryValues)) {
		summaryValues[summaryValueIndex]<-box.cox(x=summaryValues[summaryValueIndex], p=boxcoxLambda[summaryValueIndex])
	}
	#cat("\n summaryValues\n")
	#print(summaryValues)
	prunedSummaryValues<-summaryValues[which(todo>0)]
	#cat("\n prunedSummaryValues\n")
	#print(prunedSummaryValues)
	predictResult<-(predict(prunedPlsResult, matrix(prunedSummaryValues,nrow=1))$predict[, , 1])
	#cat("\n predictResult\n")
	#print(predictResult)

	predictResult
}

#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun<-function(phy,traits,intrinsicFn,extrinsicFn,summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2),startingMatrix,intrinsicMatrix,extrinsicMatrix,startingStatesGuess=c(),intrinsicValuesGuess=c(),extrinsicValuesGuess=c(),timeStep,toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, nrepSim=100, plot=T,vipthresh=0.8,epsilonProportion=0.2,epsilonMultiplier=0.5,nStepsPRC=4) {
	#print("in do run")
	splits<-getSimulationSplits(phy) #initialize this info


	#figure out number of free params
	numberParametersTotal<-dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]
	numberParametersFree<-numberParametersTotal
	numberParametersStarting<-0
	numberParametersIntrinsic<-0
	numberParametersExtrinsic<-0
	freevariables<-matrix(data=NA, nrow=2,ncol=0)
	titlevector<-c()
	freevector<-c()
		
	for (i in 1:dim(startingMatrix)[2]) {
		if (startingMatrix[1,i]== startingMatrix[2,i]) {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector,FALSE)
		}	
		else {
			numberParametersStarting<-numberParametersStarting+1
			freevariables<-cbind(freevariables,startingMatrix[,i])
			titlevector <-c(titlevector,paste("Starting", numberParametersStarting))
			freevector<-c(freevector,TRUE)
		}
	}
	for (i in 1:dim(intrinsicMatrix)[2]) {
		if (intrinsicMatrix[1,i]== intrinsicMatrix[2,i]) {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector,FALSE)
		}	
		else {
			numberParametersIntrinsic <-numberParametersIntrinsic +1
			freevariables<-cbind(freevariables, intrinsicMatrix[,i])
			titlevector <-c(titlevector,paste("Intrinsic", numberParametersIntrinsic))
			freevector<-c(freevector,TRUE)
		}
	}
	for (i in 1:dim(extrinsicMatrix)[2]) {
		if (extrinsicMatrix[1,i]== extrinsicMatrix[2,i]) {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector,FALSE)
		}	
		else {
			numberParametersExtrinsic <-numberParametersExtrinsic +1
			freevariables<-cbind(freevariables, extrinsicMatrix[,i])
			titlevector <-c(titlevector,paste("Extrinsic", numberParametersExtrinsic))
			freevector<-c(freevector,TRUE)
		}
	}

	

	#initialize guesses, if needed
	if (length(startingStatesGuess)==0) { #if no user guesses, try midpoint of the flat prior
		startingStatesGuess<-colMeans(startingMatrix)
	}
	if (length(intrinsicValuesGuess)==0) { #if no user guesses, try midpoint of the flat prior
		intrinsicValuesGuess<-colMeans(intrinsicMatrix)
	}
	if (length(extrinsicValuesGuess)==0) { #if no user guesses, try midpoint of the flat prior
		extrinsicValuesGuess<-colMeans(extrinsicMatrix)
	}
	
	
	#----------------- Find best set of summary stats to use for this problem. (start) -----------------
	#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data
	for (simIndex in 1:nrepSim) {
		trueStarting<-rep(NaN,dim(startingMatrix)[2])
		trueIntrinsic<-rep(NaN,dim(intrinsicMatrix)[2])
		trueExtrinsic<-rep(NaN,dim(extrinsicMatrix)[2])
		for (j in 1:dim(startingMatrix)[2]) {
			trueStarting[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
		}
		for (j in 1:dim(intrinsicMatrix)[2]) {
			trueIntrinsic[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
		}
		for (j in 1:dim(extrinsicMatrix)[2]) {
			trueExtrinsic[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
		}
		trueInitial<-c(trueStarting, trueIntrinsic, trueExtrinsic)
		trueFreeValues<-rbind(trueFreeValues, trueInitial[freevector])
		#cat("summaryValues\n")
		#print(summaryValues)
		summaryValues<-rbind(summaryValues,summaryStatsLong(phy,convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn= intrinsicFn,extrinsicFn= extrinsicFn,startingStates= trueStarting,intrinsicValues= trueIntrinsic,extrinsicValues= trueExtrinsic,timeStep=timeStep),phy)))
		rep(sink(),50)
		#cat("summaryValues+ summaryStatsLong\n")
		#print(summaryValues)
		#cat("rep ",simIndex," ", summaryValues[dim(summaryValues)[1],], "\n")
	}
	
	#cat("1\n")
	library("car")
	#now put this into the boxcox function to get best lambda for each summary stat
	boxcoxAddition<-10*(abs(min(summaryValues)))
	summaryValues<-boxcoxAddition+summaryValues #for boxcox, need positive numbers
	boxcoxLambda<-rep(NA,dim(summaryValues)[2])
	#cat("2 ",dim(summaryValues)[2], "\n")
	for (summaryValueIndex in 1:dim(summaryValues)[2]) {
		#cat ("3 ", summaryValueIndex, "\n")
		summary<-summaryValues[,summaryValueIndex]
		data<-list(trueFreeValues,summary)
		#cat ("4 \n")
		#print(data)
		boxcoxResult<-boxcox(trueFreeValues ~ summary, data = data, plotit=FALSE)
		#cat ("5","\n")
		#print(boxcoxResult)
		boxcoxLambda[summaryValueIndex]<-boxcoxResult$x[which(boxcoxResult$y==max(boxcoxResult$y))]
		#cat ("6 ","\n")
		#print(boxcoxLambda)
		#Transform each summary stat according to its best lambda
		summaryValues[,summaryValueIndex]<-box.cox(x=summaryValues[,summaryValueIndex], p=boxcoxLambda[summaryValueIndex])
		#cat ("7 \n")
		#print(summaryValues)
	
	}
	
	
	#Use integrOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (integromics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
	library("integrOmics")
	plsResult<-pls(Y=trueFreeValues, X=summaryValues)
	#cat("plsResult\n")
	#print(plsResult)
	vipResult<-vip(plsResult)
	#cat("vipResult\n")
	#print(vipResult)
	todo<-rep(1,dim(summaryValues)[2])#initialize the vector that indicates which summary stats to include
	for (summaryIndex in 1:dim(summaryValues)[2]) {
		if (max(vipResult[summaryIndex,])<vipthresh) {
				todo[summaryIndex]<-0 #exclude this summary stat, because it's too unimportant
		}	
	}
rep(sink(),100)
	#cat("todo","\n")
	#print(todo)
	
	prunedSummaryValues<-summaryValues[,which(todo>0)]
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

	distanceVector<-rep(NA,dim(predictResult)[1])
	for (simulationIndex in 1:dim(predictResult)[1]) {
			distanceVector[simulationIndex]<-dist(matrix(c(trueFreeValues[simulationIndex,], predictResult[simulationIndex,]),nrow=2,byrow=T))[1]
	}
	#cat("distanceVector", "\n")
	#print(distanceVector)
	densityDistanceVector<-density(distanceVector)
	plot(densityDistanceVector)
	epsilonDistance<-quantile(distanceVector,probs=epsilonProportion) #this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected 
	lines(x=c(epsilonDistance, epsilonDistance),y=c(0,max(densityDistanceVector$y)),lty="dotted")
	toleranceVector<-rep(epsilonDistance,nStepsPRC)
	for (step in 2:nStepsPRC) {
		toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
	}
	#----------------- Find distribution of distances (end) ---------------------
	
	#------------------ ABC-PRC (start) ------------------
	#do not forget to use boxcoxLambda, boxcoxAddition, and prunedPlsResult when computing distances
	nameVector<-c("generation","attempt","id","parentid","distance","weight")
	if (plot) {
		plot(x=c(min(intrinsicMatrix),max(intrinsicMatrix)),y=c(0,5*max(toleranceVector)),type="n")
	}
	for (i in 1:dim(startingMatrix)[2]) {
		nameVector<-append(nameVector,paste("StartingStates",i,sep=""))
	}
	for (i in 1:dim(intrinsicMatrix)[2]) {
		nameVector<-append(nameVector,paste("IntrinsicValue",i,sep=""))
	}
	for (i in 1:dim(extrinsicMatrix)[2]) {
		nameVector<-append(nameVector,paste("ExtrinsicValue",i,sep=""))
	}
	particleWeights=rep(0,numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
	particleParameters<-matrix(nrow=numParticles,ncol=dim(startingMatrix)[2] +  dim(intrinsicMatrix)[2] + dim(extrinsicMatrix)[2]) #stores parameters in model for each particle
	particleDistance=rep(NA,numParticles)
	particle<-1
	attempts<-0
	particleDataFrame<-data.frame()
	cat("successes","attempts","expected number of attempts required\n")
	particleVector<-c()
	originalSummaryStats<-boxcoxplsSummary(todo,summaryStatsLong(phy,traits,todo),prunedPlsResult, boxcoxLambda, boxcoxAddition)
	#cat("originalSummaryStats\n")
	#print(originalSummaryStats)
	while (particle<=numParticles) {
		attempts<-attempts+1
		
		newparticleVector<-c(new("abcparticle",id=particle,generation=1,weight=0))
		newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]],startingMatrix, intrinsicMatrix, extrinsicMatrix)
		newparticleVector[[1]]<-setDistance(newparticleVector[[1]],dist(matrix(c(boxcoxplsSummary(todo,summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),timeStep),phy),todo),prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats),nrow=2,byrow=T))[1])
		if (distance(newparticleVector[[1]])<toleranceVector[1]) {
			newparticleVector[[1]]<-setId(newparticleVector[[1]],particle)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]],1/numParticles)
			particleWeights[particle]<-1/numParticles
			particle<-particle+1
			particleVector<-append(particleVector,newparticleVector)
		}
		else {
			newparticleVector[[1]]<-setId(newparticleVector[[1]],-1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
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
			newparticleVector[[1]]<-setDistance(newparticleVector[[1]],dist(matrix(c(boxcoxplsSummary(todo,summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),timeStep),phy),todo),prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats),nrow=2,byrow=T))[1])
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
		else {
			newparticleVector[[1]]<-setId(newparticleVector[[1]],-1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
			sink()
			#print(newparticleVector)
			vectorForDataFrame<-c(dataGenerationStep,attempts,getId(newparticleVector[[1]]), particleToSelect,distance(newparticleVector[[1]]),getWeight(newparticleVector[[1]]),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]))
#cat("\n\nlength of vectorForDataFrame = ",length(vectorForDataFrame),"\n","length of startingStates = ",length(startingStates),"\nlength of intrinsicValues = ",length(intrinsicValues),"\nlength of extrinsicValues = ",length(extrinsicValues),"\ndistance = ",distance(newparticleVector[[1]]),"\nweight = ",getWeight(newparticleVector[[1]]),"\n",vectorForDataFrame,"\n")
			particleDataFrame<-rbind(particleDataFrame,data.frame(rbind(vectorForDataFrame))) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
			cat(particle-1,attempts,floor(numParticles*attempts/particle),startingStates(newparticleVector[[1]]),intrinsicValues(newparticleVector[[1]]),extrinsicValues(newparticleVector[[1]]),distance(newparticleVector[[1]]),"\n")
			
		}
		particleDataFrame[which(particleDataFrame$generation==dataGenerationStep),]$weight<-doRunOutput[which(particleDataFrame$generation==dataGenerationStep),]$weight/(sum(doRunOutput[which(particleDataFrame$generation==dataGenerationStep),]$weight))
	}
	names(particleDataFrame)<-nameVector
	if(plot) {
		quartz()
		plot(x=c(min(intrinsicMatrix),max(intrinsicMatrix)),y=c(0,1),type="n")
		for (i in 1:(length(toleranceVector)-1)) {
			graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
			lines(density(subset(particleDataFrame,generation==i)[,8]),col= graycolor)
		}
		lines(density(subset(particleDataFrame,generation==length(toleranceVector))[,8]),col= "red")
	}
	particleDataFrame

	
	
	#------------------ ABC-PRC (end) ------------------
	
	

}

presentABCOutput<-function(ABCOutput,plot=FALSE) {
	library(Hmisc)
	lastGen<-ABCOutput[which(ABCOutput$generation==max(ABCOutput$generation)),]
	finalParticles<-lastGen[which(lastGen$weight>0),]
	nParams<-dim(finalParticles)[2]-6
	nParticles<-dim(finalParticles)[1]
	resultsMatrix<-matrix(nrow=13,ncol=0)
	
	for (variable in 1:nParams) {
		resultsMatrix<-cbind(resultsMatrix,wtd.quantile(finalParticles[,6+variable],weights= nParticles*finalParticles[,6],probs=c(0,0.001, 0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 0.999,1)))
	}
	results<-data.frame(resultsMatrix)
	results
}



#test code2
#library(geiger)
#phy<-rcoal(9)
#char<-data.frame(5+sim.char(phy,model.matrix=matrix(20),1))
#Rprof()
#particledata<-abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(500,2),summaryFn= geigerUnivariateSummaryStats2)

#particledata<-abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,10),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(500,400,300, 200, 100),standardDevFactor=0.1, summaryFn= rawValuesSummaryStats,plot=T,numParticles=10)

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

#test code3
library(geiger)

#phy<-rcoal(20)
#splits<-getSimulationSplits(phy)
#char<-convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn=brownianIntrinsic,extrinsicFn=brownianExtrinsic,startingStates=c(3),intrinsicValues=.06,extrinsicValues=0,timeStep=0.001),phy)
fitContinuous(phy,char)
#particledata<-abcprc2(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,.1),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001, toleranceVector=c(100,50,25,10,5),standardDevFactor=0.05, summaryFn= rawValuesSummaryStats,plot=T,numParticles=50)

phy<-rcoal(5)

splits<-getSimulationSplits(phy)
char<-convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn=brownianIntrinsic,extrinsicFn=brownianExtrinsic,startingStates=c(3),intrinsicValues=.06,extrinsicValues=0,timeStep=0.001),phy)

#testApproach(phy=phy,originalData=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,.1),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001,standardDevFactor=0.05, plot=T,nrepSim=3,startingStatesGuess=c(2),intrinsicValuesGuess=c(0.06),extrinsicValuesGuess=c(0))

doRunOutput<-doRun(phy=phy,traits=char,intrinsicFn= brownianIntrinsic,extrinsicFn= brownianExtrinsic,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,.1),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001,standardDevFactor=0.05, plot=T,nrepSim=10,startingStatesGuess=c(2),intrinsicValuesGuess=c(0.06),extrinsicValuesGuess=c(0),epsilonProportion=0.2,epsilonMultiplier=0.5,nStepsPRC=4,numParticles=10)