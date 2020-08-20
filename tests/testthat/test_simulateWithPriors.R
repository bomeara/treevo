test_that("simulateWithPriors works", {
	set.seed(1)
	simPhyExample <- rcoal(5)
	simPhyExample$edge.length <- simPhyExample$edge.length * 20

	#expect_warning(
	simData <- simulateWithPriors(
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
		freevector = NULL,
		giveUpAttempts = 10,
		verbose = FALSE
		)


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

	expect_is(simData, "numeric")
	expect_is(simDataParallel, "matrix")
	
	})
