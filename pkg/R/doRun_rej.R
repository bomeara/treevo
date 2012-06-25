
##TreeYears = 1000 if tree is in thousands of years


doRun_rej<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingValuesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), TreeYears=1e+04, standardDevFactor=0.20, StartSims=NA, jobName=NA, trueStartingState=NA, trueIntrinsicState=NA, abcMethod="rejection", vipthresh=0.8, abcTolerance=0.1, multicore=FALSE, coreLimit=NA) {

	if (!is.binary.tree(phy)) {
		print("Warning: Tree is not fully dichotomous")
	}
	startTime<-proc.time()[[3]]
	timeStep<-1/TreeYears			
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
	
	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]) {
		namesForPriorMatrix<-c(paste("StartingStates", a, sep=""))
	}
	for (b in 1:dim(intrinsicPriorsValues)[2]) {
		namesForPriorMatrix<-append(namesForPriorMatrix, paste("IntrinsicValue", b, sep=""))
	}
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste("ExtrinsicValue", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")
		
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
	if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
		startingValuesGuess<-rep(NA,length(startingPriorsFns))
		for (i in 1:length(startingPriorsFns)) {
			startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
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
	
	if (is.na(StartSims)) {
		StartSims<-1000*numberParametersFree
	}
	
	nrepSim<-StartSims #Used to be multiple tries where nrepSim = StartSims*((2^try)/2).  If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
	cat(paste("Number of initial simulations set to", nrepSim, "\n"))
	
	trueFreeValuesANDSummaryValues<-parallelSimulation(nrepSim, coreLimit, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore)
	
	save(trueFreeValuesANDSummaryValues, file="tFVandSV.Rdata")
	cat("\n\n")
	simTime<-proc.time()[[3]]-startTime
	cat(paste("Initial simulations took", round(simTime, digits=3), "seconds"), "\n")
	
	#separate the simulation results: true values and the summary values
	trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
	
	#boxcox transform summary values
	boxcoxEstimates<-boxcoxEstimation(summaryValuesMatrix)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValuesMatrix<-boxcoxEstimates$boxcoxSummaryValuesMatrix
	boxcoxOriginalSummaryStats<-boxcoxTransformation(summaryStatsLong(phy=phy, data=traits), boxcoxAddition, boxcoxLambda)

	library("mixOmics")
	getVipSingleColumn<-function(trueFreeValuesColumn, boxcoxSummaryValuesMatrix) {
		return(vip(pls(X=boxcoxSummaryValuesMatrix, Y=trueFreeValuesColumn, ncomp=1)))
	}
	
	#note this has vip for summary val 5, for true param 1, in row 5, column 1, and so forth
	allVip<-apply(trueFreeValuesMatrix, 2, getVipSingleColumn, boxcoxSummaryValuesMatrix)
	
	#this will have which summary vals have importance above the threshold
	whichVip<-(allVip>vipthresh)
	
	abcDistancesRaw<-matrix(nrow=dim(trueFreeValuesMatrix)[1], ncol=dim(trueFreeValuesMatrix)[2])  #rep(0,dim(trueFreeValuesMatrix)[1]) #will hold the distances for each particle
	#abcDistancesRawTotal<-vector(length=dim(trueFreeValuesMatrix)[1])
	#now we go true parameter by true parameter, using only the summary values with enough importance for each
	#we get a distance for each particle from the observed particle for each true param
	#then get a euclidean distance for all of these
	#it's like getting delta X, then delta Y, and figuring out distance from the origin using
	#sqrt((X-0)^2 + (Y-0)^2)
	for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
		abcDistancesRaw[,freeParamIndex]<-abc(target=boxcoxOriginalSummaryStats[whichVip[,freeParamIndex]], param=trueFreeValuesMatrix[,freeParamIndex], sumstat= boxcoxSummaryValuesMatrix[,whichVip[,freeParamIndex] ], tol=1, method=abcMethod)$dist^2 #because we're doing euclidean distance, from observed, which has value 0, 0, 0, etc.
		#abcDistancesRawTotal<-abcDistancesRawTotal+abcDistancesRaw
	}
	abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
	abcResults<-vector("list")
	abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
	abcResults$unadj.values<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
	abcResults$dist<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
	print(abcResults)

	particleDataFrame<-data.frame(cbind(rep(1, dim(abcResults$unadj.values)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(abcResults$unadj.values)[1]), rep(0, dim(abcResults$unadj.values)[1]), abcResults$dist, rep(1, dim(abcResults$unadj.values)[1]), abcResults$unadj.values  ))
	colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2])))

	save(abcDistancesRaw, abcDistancesRawTotal, abcDistances, abcResults, particleDataFrame, file="BarbsDistances.Rdata")
	input.data<-rbind(jobName, length(phy[[3]]), timeStep, StartSims, standardDevFactor, abcMethod, abcTolerance)

	rejectionResults<-vector("list")
	#names(rejectionResults)<-c("input.data", "PriorMatrix", "phy", "traits")

	rejectionResults$input.data<-input.data
	rejectionResults$PriorMatrix<-PriorMatrix
	rejectionResults$phy<-phy
	rejectionResults$traits<-traits
	rejectionResults$trueFreeValuesANDSummaryValues<-trueFreeValuesANDSummaryValues
  	rejectionResults$SimTime<-simTime
  	rejectionResults$vipSumStats<-whichVip
  	rejectionResults$abcDistancesRaw<-abcDistancesRaw
	rejectionResults$particleDataFrame<-particleDataFrame
  
	return(rejectionResults)
}





