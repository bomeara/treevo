test_that("methodsPLS works", {

  set.seed(1)
  simPhyExample <- rcoal(5)
  simPhyExample$edge.length <- simPhyExample$edge.length * 20

  nSimulations <- 6

  #expect_warning(
  simDataParallel <- parallelSimulateWithPriors(
    nrepSim = nSimulations,
    multicore = FALSE,
    coreLimit = 1,
    phy = simPhyExample,
    intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal",
	startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1]))), 
	intrinsicPriorsFns = c("exponential"),
	intrinsicPriorsValues = list(10),
	extrinsicPriorsFns = c("fixed"),
	extrinsicPriorsValues = list(0),
    generation.time = 300000,
    checkpointFile = NULL,
    checkpointFreq = 24,
    verbose = FALSE,
    freevector = NULL,
    taxonDF = NULL
    #,niter.brown = 25,
    #niter.lambda = 25,
    #niter.delta = 25,
    #niter.OU = 25,
    #niter.white = 25
    )
    #)
    
  nParFree <- sum(attr(simDataParallel, "freevector"))
  trueFreeValuesMat <- simDataParallel[, 1:nParFree]
  summaryValuesMat <- simDataParallel[, -1:-nParFree]

  expect_warning(
  PLSmodel <- returnPLSModel(
    trueFreeValuesMatrix = trueFreeValuesMat,
    summaryValuesMatrix = summaryValuesMat,
    validation = "CV",
    scale = TRUE,
    variance.cutoff = 95,
    segments = nSimulations)
   )
  expect_silent(
    PLSTransform(
      summaryValuesMatrix = summaryValuesMat,
      pls.model = PLSmodel
      )
  )
    
})
