test_that("doRun_prc runs correctly", {
	#
    data(simRunExample)
    set.seed(1)
    expect_warning(
	#
    results <- doRun_prc(
		phy = simPhy,
		traits = simChar,
		intrinsicFn=brownianIntrinsic,
		extrinsicFn=nullExtrinsic,
		startingPriorsFns="normal",
		startingPriorsValues = list(c(mean(simChar[, 1]), sd(simChar[, 1]))), 
		intrinsicPriorsFns = c("exponential"),
		intrinsicPriorsValues = list(10),
		extrinsicPriorsFns = c("fixed"),
		extrinsicPriorsValues = list(0),
		generation.time=1000000,
		nRuns=1,
		nStepsPRC=2,
		numParticles=10,
		nInitialSimsPerParam=5,
		jobName="exampleRun",
		diagnosticPRCmode=TRUE,
		verboseParticles=TRUE,
		stopRule=FALSE,
		multicore=FALSE,
		coreLimit=1
		)
	#  
    )
    expect_output(print(results))
    expect_is(results, "list")
    expect_warning(expect_false(any(is.na(pairwiseESS(results$particleDataFrame)))))
	#
})



test_that("doRun_rej works", {
  data(simRunExample)
    set.seed(1)
   expect_warning(
   resultsRej <- doRun_rej(
    phy = simPhy,
    traits = simChar,
    intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal",
	startingPriorsValues = list(c(mean(simChar[, 1]), sd(simChar[, 1]))), 
	intrinsicPriorsFns = c("exponential"),
	intrinsicPriorsValues = list(10),
	extrinsicPriorsFns = c("fixed"),
	extrinsicPriorsValues = list(0),
    nInitialSims = 10,
    generation.time=1000000,
    jobName = "examplerun_rej",
    abcTolerance = 0.05,
    multicore = FALSE,
    coreLimit = 1)
    )
    expect_is(resultsRej, "list")
})

