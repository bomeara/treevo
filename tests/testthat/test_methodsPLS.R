test_that("methodsPLS works", {

  set.seed(1)
  simPhy <- rcoal(5)
  simPhy$edge.length <- simPhy$edge.length * 20

  nSimulations <- 6

  #expect_warning(
  simDataParallel <- parallelSimulateWithPriors(
    nrepSim = nSimulations,
    multicore = FALSE,
    coreLimit = 1,
    phy = simPhy,
    intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal",
    startingPriorsValues = matrix(c(mean(simChar[,
      1]), sd(simChar[, 1]))),
    intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,
      byrow = FALSE),
    extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,
      byrow = FALSE),
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

  PLSTransform(
    summaryValuesMatrix = summaryValuesMat,
    pls.model = PLSmodel
    )
    
})
