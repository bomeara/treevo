test_that("boxcoxTransformation works", {
  set.seed(1)
  simPhy <- rcoal(5)
  simPhy$edge.length <- simPhy$edge.length * 20
  
  
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
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2, byrow = FALSE), 
	extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2, byrow = FALSE), 
	timeStep = 1e-04, 
	checkpointFile = NULL,
    checkpointFreq = 24, 
	verbose = FALSE, 
	freevector = NULL,
    taxon.df = NULL 
	#,niter.brown = 25, niter.lambda = 25,
    #niter.delta = 25, niter.OU = 25, niter.white = 25
	)
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
