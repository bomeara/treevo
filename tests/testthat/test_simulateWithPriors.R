test_that("simulateWithPriors works", {
  simPhy <- rcoal(30)
  simPhy$edge.length <- simPhy$edge.length * 20
  
  simData <- simulateWithPriors(phy = simPhy, intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, startingPriorsFns = "normal",
    startingPriorsValues = matrix(c(mean(simChar[,
      1]), sd(simChar[, 1]))), intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,
      byrow = FALSE), extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,
      byrow = FALSE), generation.time = 1e+05, freevector = NULL,
    giveUpAttempts = 10, verbose = FALSE
	#,niter.brown = 25,niter.lambda = 25, niter.delta = 25, niter.OU = 25,niter.white = 25
	)
  simData
  simDataParallel <- parallelSimulateWithPriors(
    nrepSim = 2,
    multicore = FALSE, 
	coreLimit = 1, 
	phy = simPhy,
    intrinsicFn = brownianIntrinsic, 
	extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal", 
	startingPriorsValues = matrix(c(mean(simChar[,1]), sd(simChar[, 1]))), 
	intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,byrow = FALSE), 
	extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,byrow = FALSE), 
	generation.time = 1e+05, 
	checkpointFile = NULL,
    checkpointFreq = 24, 
	verbose = FALSE, 
	freevector = NULL,
    taxon.df = NULL
	#,niter.brown = 25, niter.lambda = 25,niter.delta = 25, niter.OU = 25, niter.white = 25
	)
  simDataParallel
  
})
