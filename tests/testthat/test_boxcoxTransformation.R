test_that("boxcoxTransformation works", {
  set.seed(1)
  simPhy <- rcoal(5)
  simPhy$edge.length <- simPhy$edge.length * 20

  #expect_warning(
  simDataParallel <- parallelSimulateWithPriors(
    nrepSim = 2,
    multicore = FALSE,
    coreLimit = 1,
    phy = simPhy,
    intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal",
    startingPriorsValues = list(c(mean(simChar[, 1]), sd(simChar[, 1]))), 
    intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = list(10),
    extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = list(0),
    generation.time = 100000,
    checkpointFile = NULL,
    checkpointFreq = 24,
    verbose = FALSE,
    freevector = NULL,
    taxonDF = NULL
    #,niter.brown = 25, niter.lambda = 25,
    #niter.delta = 25, niter.OU = 25, niter.white = 25
    )
    #)
    
  nParFree <- sum(attr(simDataParallel, "freevector"))
  summaryValuesMat <- simDataParallel[, -1:-nParFree]
  #
  boxTranMat <- boxcoxTransformationMatrix(summaryValuesMatrix = summaryValuesMat)
  #
  result<-boxcoxTransformation(
    summaryValuesVector = summaryValuesMat[,1],
    boxcoxAddition = boxTranMat$boxcoxAddition,
    boxcoxLambda = boxTranMat$boxcoxLambda)
})
