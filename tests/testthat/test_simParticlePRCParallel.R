test_that("simParticlePRCParallel run correctly", {
    data(simRunExample)
    set.seed(1)
    #
    phy = simPhyExample
    traits = simCharExample
    intrinsicFn=brownianIntrinsic
    extrinsicFn=nullExtrinsic
    startingPriorsFns="normal"
	startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1])))
	intrinsicPriorsFns = c("exponential")
	intrinsicPriorsValues = list(10)
	extrinsicPriorsFns = c("fixed")
	extrinsicPriorsValues = list(0)
    generation.time=1000000
    standardDevFactor=0.2
    nInitialSims=10
    epsilonProportion=0.7
    epsilonMultiplier=0.7
    nStepsPRC=2
    numParticles=3
    jobName="exampleRun"
    stopRule=FALSE
    multicore=FALSE
    coreLimit=1
    TreeYears=max(branching.times(phy)) * 1e6
    taxonDF <- getTaxonDFWithPossibleExtinction(phy)
    timeStep<-generation.time/TreeYears
    edgesRescaled<-phy$edge.length/max(node.depth.edgelength(phy))
    totalGenerations<-sum(sapply(edgesRescaled,function(x) floor(x/timeStep)))
    saveData<-FALSE
    validation="CV" 
    scale=TRUE
    variance.cutoff=95    
    numParticles=3
    dataGenerationStep=1    
    #
    ################
    # CODE TO FUNCTION
    ################
    freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues,
                        intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
                        extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
    numberParametersTotal<-length(freevector)
    numberParametersFree<-sum(freevector)
	namesParFree <- names(freevector)[freevector]
    #
    # get prior list
    priorList<-getPriorList(
        startingPriorsValues=startingPriorsValues,
        intrinsicPriorsValues=intrinsicPriorsValues,
        extrinsicPriorsValues=extrinsicPriorsValues,
        startingPriorsFns=startingPriorsFns,
        intrinsicPriorsFns=intrinsicPriorsFns,
        extrinsicPriorsFns=extrinsicPriorsFns,
        numberParametersTotal=numberParametersTotal)
    #    
    ##    
    #initialize weighted mean sd matrices
    weightedMeanParam <- matrix(nrow = nStepsPRC, 
        ncol = numberParametersFree)
    param.stdev <- matrix(nrow = nStepsPRC, ncol = numberParametersFree)
    #
	colnames(weightedMeanParam) <- colnames(param.stdev) <- namesParFree 
    rownames(weightedMeanParam) <- rownames(param.stdev) <- paste0("Gen ", c(1: nStepsPRC), sep = "") 
    #
    #
    # save input data for use later
    #input.data<-rbind(jobName, Ntips = Ntip(phy), nInitialSims, generation.time, TreeYears, timeStep, totalGenerations,
     #   epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
    #
    # get summary values for observed data
    originalSummaryValues<-summaryStatsLong(
        phy=phy, 
        traits=traits
        #niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200
        )
    #    
    # INITIAL SIMULATIONS
    initialSimsRes<-initialSimsPRC(
        nrepSim=nInitialSims, 
        phy=phy,  
        taxonDF=taxonDF,
        startingPriorsValues=startingPriorsValues, 
        intrinsicPriorsValues=intrinsicPriorsValues, 
        extrinsicPriorsValues=extrinsicPriorsValues,
        startingPriorsFns=startingPriorsFns, 
        intrinsicPriorsFns=intrinsicPriorsFns, 
        extrinsicPriorsFns=extrinsicPriorsFns,
        freevector=freevector, 
        timeStep=timeStep, 
        intrinsicFn=intrinsicFn, 
        extrinsicFn=extrinsicFn, 
        nStepsPRC=nStepsPRC,
        coreLimit=coreLimit, 
        multicore=multicore,
        numberParametersFree=numberParametersFree,
        originalSummaryValues=originalSummaryValues,
        epsilonProportion=epsilonProportion,
        epsilonMultiplier=epsilonMultiplier,
        validation=validation,
        variance.cutoff=variance.cutoff,
        saveData=saveData,
        jobName=jobName
        )        
    #
    # pull the PLS model list out
    pls.model.list<-initialSimsRes$pls.model.list
    #
    toleranceValue<-initialSimsRes$toleranceVector[dataGenerationStep]
    #
    #
    results <- simParticlePRCParallel(
        nSim=3
        ,multicore=FALSE
        ,coreLimit=1
        ,phy=phy
        ,taxonDF=taxonDF
        , timeStep=timeStep
        ,intrinsicFn=intrinsicFn
        , extrinsicFn=extrinsicFn 
        ,startingPriorsValues=startingPriorsValues
        ,intrinsicPriorsValues=intrinsicPriorsValues
        ,extrinsicPriorsValues=extrinsicPriorsValues
        ,startingPriorsFns=startingPriorsFns
        ,intrinsicPriorsFns=intrinsicPriorsFns
        ,extrinsicPriorsFns=extrinsicPriorsFns
        ,originalSummaryValues=originalSummaryValues
        ,pls.model.list=pls.model.list
        ,toleranceValue=toleranceValue
        ,prevGenParticleList=NULL
        ,standardDevFactor=standardDevFactor
        ,numParticles=numParticles
        )
    expect_is(results, "list")
})
    

test_that("simParticlePRC run correctly", {
    data(simRunExample)
    set.seed(1)
    #
    phy = simPhyExample
    traits = simCharExample
    intrinsicFn=brownianIntrinsic
    extrinsicFn=nullExtrinsic
    startingPriorsFns="normal"
	startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1]))) 
	intrinsicPriorsFns = c("exponential")
	intrinsicPriorsValues = list(10)
	extrinsicPriorsFns = c("fixed")
	extrinsicPriorsValues = list(0)
    generation.time=1000000
    standardDevFactor=0.2
    nInitialSims=10
    epsilonProportion=0.7
    epsilonMultiplier=0.7
    nStepsPRC=2
    numParticles=3
    jobName="exampleRun"
    stopRule=FALSE
    multicore=FALSE
    coreLimit=1
    TreeYears=max(branching.times(phy)) * 1e6
    taxonDF <- getTaxonDFWithPossibleExtinction(phy)
    timeStep<-generation.time/TreeYears
    edgesRescaled<-phy$edge.length/max(node.depth.edgelength(phy))
    totalGenerations<-sum(sapply(edgesRescaled,function(x) floor(x/timeStep)))
    saveData<-FALSE
    validation="CV" 
    scale=TRUE
    variance.cutoff=95    
    numParticles=3
    dataGenerationStep=1    
    #
    ################
    # CODE TO FUNCTION
    ################
    freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues,
                        intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
                        extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
    numberParametersTotal<-length(freevector)
    numberParametersFree<-sum(freevector)
	namesParFree <- names(freevector)[freevector]
    #
    # get prior matrix
    priorList<-getPriorList(
        startingPriorsValues=startingPriorsValues,
        intrinsicPriorsValues=intrinsicPriorsValues,
        extrinsicPriorsValues=extrinsicPriorsValues,
        startingPriorsFns=startingPriorsFns,
        intrinsicPriorsFns=intrinsicPriorsFns,
        extrinsicPriorsFns=extrinsicPriorsFns,
        numberParametersTotal=numberParametersTotal)
    #    
    ##    
    #initialize weighted mean sd matrices
    weightedMeanParam <- matrix(nrow = nStepsPRC, 
        ncol = numberParametersFree)
    param.stdev <- matrix(nrow = nStepsPRC, ncol = numberParametersFree)
    #
	colnames(weightedMeanParam) <- colnames(param.stdev) <- namesParFree 
    rownames(weightedMeanParam) <- rownames(param.stdev) <- paste0("Gen ", c(1: nStepsPRC), sep = "") 
    #
    # save input data for use later
    input.data<-rbind(jobName, Ntips = Ntip(phy), nInitialSims, generation.time, TreeYears, timeStep, totalGenerations,
        epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
    #
    # get summary values for observed data
    originalSummaryValues<-summaryStatsLong(
        phy=phy, 
        traits=traits
        #niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200
        )
    #    
    # INITIAL SIMULATIONS
    initialSimsRes<-initialSimsPRC(
        nrepSim=nInitialSims, 
        phy=phy,  
        taxonDF=taxonDF,
        startingPriorsValues=startingPriorsValues, 
        intrinsicPriorsValues=intrinsicPriorsValues, 
        extrinsicPriorsValues=extrinsicPriorsValues,
        startingPriorsFns=startingPriorsFns, 
        intrinsicPriorsFns=intrinsicPriorsFns, 
        extrinsicPriorsFns=extrinsicPriorsFns,
        freevector=freevector, 
        timeStep=timeStep, 
        intrinsicFn=intrinsicFn, 
        extrinsicFn=extrinsicFn, 
        nStepsPRC=nStepsPRC,
        coreLimit=coreLimit, 
        multicore=multicore,
        numberParametersFree=numberParametersFree,
        originalSummaryValues=originalSummaryValues,
        epsilonProportion=epsilonProportion,
        epsilonMultiplier=epsilonMultiplier,
        validation=validation,
        variance.cutoff=variance.cutoff,
        saveData=saveData,
        jobName=jobName
        )        
    #
    # pull the PLS model list out
    pls.model.list<-initialSimsRes$pls.model.list
    #
    toleranceValue<-initialSimsRes$toleranceVector[dataGenerationStep]
    #
    #
    # now test simParticle PRC
    results <- simParticlePRC<-function(
        phy=phy
        ,taxonDF=taxonDF
        ,timeStep=timeStep
        ,intrinsicFn=intrinsicFn
        ,extrinsicFn=extrinsicFn 
        ,startingPriorsValues=startingPriorsValues
        ,intrinsicPriorsValues=intrinsicPriorsValues
        ,extrinsicPriorsValues=extrinsicPriorsValues
        ,startingPriorsFns=startingPriorsFns
        ,intrinsicPriorsFns=intrinsicPriorsFns
        ,extrinsicPriorsFns=extrinsicPriorsFns
        ,originalSummaryValues=originalSummaryValues
        ,pls.model.list=pls.model.list
        ,toleranceValue=toleranceValue
        ,prevGenParticleList=NULL
        ,standardDevFactor=standardDevFactor
        ,numParticles=numParticles
        )
    expect_is(results, "list")
})

