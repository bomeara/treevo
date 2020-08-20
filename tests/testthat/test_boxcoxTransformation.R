test_that("boxcoxTransformation works", {
	set.seed(1)
	simPhyExample <- rcoal(5)
	simPhyExample$edge.length <- simPhyExample$edge.length * 20

	#expect_warning(
	simDataParallel <- parallelSimulateWithPriors(
		nrepSim = 2,
		multicore = FALSE,
		coreLimit = 1,
		phy = simPhyExample,
		intrinsicFn = brownianIntrinsic,
		extrinsicFn = nullExtrinsic,
		startingPriorsFns = "normal",
		startingPriorsValues = list(c(
			mean(simCharExample[, 1]), 
			sd(simCharExample[, 1])
			)), 
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
		)

	nParFree <- sum(attr(simDataParallel, "freevector"))
	summaryValuesMat <- simDataParallel[, -1:-nParFree]

	boxTranMat <- boxcoxTransformationMatrix(summaryValuesMatrix = summaryValuesMat)

	result<-boxcoxTransformation(
		summaryValuesVector = summaryValuesMat[,1],
		boxcoxAddition = boxTranMat$boxcoxAddition,
		boxcoxLambda = boxTranMat$boxcoxLambda
		)

	expect_is(result, "numeric")
	})
