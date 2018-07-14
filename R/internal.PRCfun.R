# internal side functions for use in doRun_prc

# multicore simSumDistancePRC 
simParticlePRCParallel <- function(
    nSim, multicore, coreLimit
    , phy, taxonDF, timeStep 
    , intrinsicFn, extrinsicFn 
    , startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues
    , startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns
    #startingValues, intrinsicValues, extrinsicValues
    , originalSummaryValues, pls.model.list
    , toleranceValue, prevGenParticleList
    , standardDevFactor, numParticles){
    #
    # set up multicore
    cluster <- setupMulticore(multicore, nSim = nSim, coreLimit = coreLimit)
    #        
    repDistFE <- foreach(1:nSim, .combine = c)
    #
    # need a function for parallel doSimulation, and all other associated particle actions
    newParticleList <- (    #makeQuiet(
        repDistFE %dopar% simParticlePRC(
            phy = phy, taxonDF = taxonDF, timeStep = timeStep, 
            intrinsicFn = intrinsicFn, 
            extrinsicFn = extrinsicFn, 
            startingPriorsValues = startingPriorsValues, 
            startingPriorsFns = startingPriorsFns, 
            intrinsicPriorsValues = intrinsicPriorsValues, 
            intrinsicPriorsFns = intrinsicPriorsFns, 
            extrinsicPriorsValues = extrinsicPriorsValues, 
            extrinsicPriorsFns = extrinsicPriorsFns, 
            originalSummaryValues = originalSummaryValues, 
            pls.model.list = pls.model.list, 
            toleranceValue = toleranceValue, 
            prevGenParticleList = prevGenParticleList, 
            standardDevFactor = standardDevFactor, 
            numParticles = numParticles
            )
        #)
        )
    # stop multicore processes
    stopMulticore(cluster)
    #
    return(newParticleList)
    }    

    
# multicore simSumDistancePRC 
simParticlePRCNotParallel <- function(
    nSim, multicore, coreLimit
    , phy, taxonDF, timeStep 
    , intrinsicFn, extrinsicFn 
    , startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues
    , startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns
    #startingValues, intrinsicValues, extrinsicValues
    , originalSummaryValues, pls.model.list
    , toleranceValue, prevGenParticleList
    , standardDevFactor, numParticles){
    #
    newParticleList <- list()
    #
    for(i in 1:nSim){
        newParticleList <- c(newParticleList, 
            simParticlePRC(
                phy = phy, taxonDF = taxonDF, timeStep = timeStep, 
                intrinsicFn = intrinsicFn, 
                extrinsicFn = extrinsicFn, 
                startingPriorsValues = startingPriorsValues, 
                startingPriorsFns = startingPriorsFns, 
                intrinsicPriorsValues = intrinsicPriorsValues, 
                intrinsicPriorsFns = intrinsicPriorsFns, 
                extrinsicPriorsValues = extrinsicPriorsValues, 
                extrinsicPriorsFns = extrinsicPriorsFns, 
                originalSummaryValues = originalSummaryValues, 
                pls.model.list = pls.model.list, 
                toleranceValue = toleranceValue, 
                prevGenParticleList = prevGenParticleList, 
                standardDevFactor = standardDevFactor, 
                numParticles = numParticles
                )
            )
        }
    #
    return(newParticleList)
    }    
    

# internal function for simulating and evaulating particles
    # for doRun_PRC - to be run in parallel
simParticlePRC <- function(
    phy, taxonDF, timeStep, intrinsicFn, extrinsicFn, 
    startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, 
    startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, 
    #startingValues, intrinsicValues, extrinsicValues, 
    originalSummaryValues, pls.model.list
    , toleranceValue, prevGenParticleList, standardDevFactor, numParticles){
    # 
    # get particle parameters
    #
    # can test for dataGenerationStep = 1 if prevGenParticleList is null
    if(is.null(prevGenParticleList)){
        # i.e. if dataGenerationStep  ==  1
            # get new particles from priors
        # for generation = 1, weight all particles the same
        newparticle <- abcparticle(id = NA, generation = 1, weight = 1/numParticles)    #, parentid = NA
        newparticle <- initializeStatesFromMatrices(
            particle = newparticle, 
            startingPriorsValues = startingPriorsValues, 
            startingPriorsFns = startingPriorsFns, 
            intrinsicPriorsValues = intrinsicPriorsValues, 
            intrinsicPriorsFns = intrinsicPriorsFns, 
            extrinsicPriorsValues = extrinsicPriorsValues, 
            extrinsicPriorsFns = extrinsicPriorsFns
            )                            
    }else{
        #print("check 1 - get weighted params for gen >1")
        prevGenParticleWeights <- sapply(prevGenParticleList, function(x) x$weight)
        # use particles from PREVIOUS GENERATION to randomly select a particle
        particleToSelect <- which.max(as.vector(rmultinom(1, size = 1, prob = prevGenParticleWeights)))
        # get that particle's data
        oldparticle <- prevGenParticleList[[particleToSelect]]
        #
        #print("check 1b - mutate")
        #print(list(particle = oldparticle, 
        #    startingPriorsValues = startingPriorsValues, startingPriorsFns = startingPriorsFns, 
        #    intrinsicPriorsValues = intrinsicPriorsValues, intrinsicPriorsFns = intrinsicPriorsFns, 
        #    extrinsicPriorsValues = extrinsicPriorsValues, extrinsicPriorsFns = extrinsicPriorsFns, 
        #    standardDevFactor = standardDevFactor
        #    ))
        #    
        newparticle <- mutateStates(particle = oldparticle, 
            startingPriorsValues = startingPriorsValues, startingPriorsFns = startingPriorsFns, 
            intrinsicPriorsValues = intrinsicPriorsValues, intrinsicPriorsFns = intrinsicPriorsFns, 
            extrinsicPriorsValues = extrinsicPriorsValues, extrinsicPriorsFns = extrinsicPriorsFns, 
            standardDevFactor = standardDevFactor
            )
        newparticle$parentid <- particleToSelect
        #print("check 2a - finished mutate")
        }
    #
    #print("check 2B - startSim")
    #
    # do the simulation
    simTraitsParticle <- doSimulationInternal(
        taxonDF = taxonDF, 
        intrinsicFn = intrinsicFn, 
        extrinsicFn = extrinsicFn, 
        startingValues = newparticle$startingValues, 
        intrinsicValues = newparticle$intrinsicValues, 
        extrinsicValues = newparticle$extrinsicValues, 
        timeStep = timeStep
        )
    #
    #print("check 3 - endSim")
    #
    # get the summary stats    
    simSumMat <- summaryStatsLong(phy = phy, traits = simTraitsParticle)
    # get the distance of the simulation to the original
    simDistance <- abcDistance(summaryValuesMatrix = simSumMat, 
        originalSummaryValues = originalSummaryValues, 
        pls.model.list = pls.model.list)    
    #
    #print("check 4 - got sumstat")
    #
    # get the weights, if it passes the tolerance
    if ((simDistance) < toleranceValue) {
        # record the distance
        newparticle$distance <- simDistance
            #ID doesn't need to be set - that's just the ID of the particle...
        #newparticle$id <- particle
        #particle <- particle+1
        #
        if(!is.null(prevGenParticleList)){
            #for generation>1 - now get weights, using correction in Sisson et al. 2007
            newparticle$weight <- sumLogTranProb(
                 prevGenParticleList = prevGenParticleList
                , newStartingValues = newparticle$startingValues
                , newIntrinsicValues = newparticle$intrinsicValues
                , newExtrinsicValues = newparticle$extrinsicValues
                , startingPriorsFns = startingPriorsFns
                , intrinsicPriorsFns = intrinsicPriorsFns
                , extrinsicPriorsFns = extrinsicPriorsFns
                , startingPriorsValues = startingPriorsValues
                , intrinsicPriorsValues = intrinsicPriorsValues
                , extrinsicPriorsValues = extrinsicPriorsValues
                , standardDevFactor = standardDevFactor
                )
            }
    }else{
        # particle didn't pass, discard it
        newparticle <- NA
        }    
    newparticle <- list(newparticle)
    return(newparticle)
    }
    
boostNsim <- function(nSims, nCores){
    newNsim <- nCores*(nSims%/%nCores + as.logical(nSims%%nCores))
    return(newNsim)
    }

getlnTransitionProb <- function(newvalue, meantouse, Fn, priorValues, stdFactor){
        #newvalue = newparticleList[[1]]$startingValues[j], 
        #meantouse = prevGenParticleList[[i]]$startingValues[j], 
        #Fn = startingPriorsFns[j], 
        #priorValues =  startingPriorsValues[, j], 
        #stdFactor = standardDevFactor
                                        
    if (Fn == "uniform") {
        sdtouse <- stdFactor*((max(priorValues)-min(priorValues))/sqrt(12))
        #message(paste0("Fn is uniform and sdtouse = ", sdtouse))
    }
    else if (Fn == "exponential") {
        sdtouse <- stdFactor*(1/priorValues[1])
        #message(paste0("Fn is exponential and sdtouse = ", sdtouse))
    }
    else {
        sdtouse <- stdFactor*(priorValues[2])
    }
    #
    lnlocalTransitionProb <- dnorm(newvalue, mean = meantouse, sd = sdtouse, log = TRUE
        ) - ((log(1)/pnorm(min(priorValues), mean = meantouse, sd = sdtouse, lower.tail = TRUE, log.p = TRUE))
            * pnorm(max(priorValues), mean = meantouse , sd = sdtouse, lower.tail = FALSE, log.p = TRUE))
    if(length(lnlocalTransitionProb) != 1){
        #message(lnlocalTransitionProb)
        stop("Somehow, multiple lnlocalTransitionProb values produced")
        }
    if (is.nan(lnlocalTransitionProb)) {  #to prevent lnlocalTransitionProb from being NaN (if pnorm = 0)
        lnlocalTransitionProb <- .Machine$double.xmin
    }
    if (min(priorValues) == max(priorValues)) {
        lnlocalTransitionProb = log(1)
    }
    if(!is.finite(lnlocalTransitionProb) || is.na(lnlocalTransitionProb)) {
        message(paste0("issue with lnlocalTransitionProb = ", lnlocalTransitionProb))
        }
    return(lnlocalTransitionProb)
    }    

sumLogTranProb <- function(prevGenParticleList
    , newStartingValues, newIntrinsicValues, newExtrinsicValues
    , startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns
    , startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues
    , standardDevFactor){
    #now get weights, using correction in Sisson et al. 2007
    #
    newWeight = 0
    #
    for (i in 1:length(prevGenParticleList)) {
        #
        LLTPstart <- sapply(length(newStartingValues), 
            function(j) getlnTransitionProb(
                newvalue = newStartingValues[j], 
                meantouse = prevGenParticleList[[i]]$startingValues[j], 
                Fn = startingPriorsFns[j], 
                priorValues =  startingPriorsValues[, j], 
                stdFactor = standardDevFactor
                )
            )
        LLTPintr <- sapply(length(newIntrinsicValues), 
            function(j) getlnTransitionProb(
                newvalue = newIntrinsicValues[j], 
                meantouse = prevGenParticleList[[i]]$intrinsicValues[j], 
                Fn = intrinsicPriorsFns[j], 
                priorValues =  intrinsicPriorsValues[, j], 
                stdFactor = standardDevFactor
                )
            )
        LLTPextr <- sapply(length(newExtrinsicValues), 
            function(j) getlnTransitionProb(
                newvalue = newExtrinsicValues[j], 
                meantouse = prevGenParticleList[[i]]$extrinsicValues[j], 
                Fn = extrinsicPriorsFns[j], 
                priorValues =  extrinsicPriorsValues[, j], 
                stdFactor = standardDevFactor
                )
            )
        #
        lnTransitionProb = log(1)
        #
        lnTransitionProb <- lnTransitionProb+sum(LLTPstart)+sum(LLTPintr)+sum(LLTPextr)
        #
        #if(!is.finite(lnTransitionProb) || is.na(lnTransitionProb)) {
        #    warning(paste0("Issue with lnTransitionProb: ", 
        #        " lnTransitionProb = ", lnTransitionProb))
        #    }
        #
        newWeight <- newWeight+(prevGenParticleList[[i]]$weight)*exp(lnTransitionProb)
        } #for (i in 1:length(prevGenParticleList)) bracket
    #
    if (!is.finite(newWeight) || is.na(newWeight)) {
        warning(paste0("Possible error; newWeight is ", newWeight))
        }
    return(newWeight)
    }
    


