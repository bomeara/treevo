test_that("doRun_prc runs correctly", {
	data(simRunExample)
	set.seed(1)
	expect_warning(
	results <- doRun_prc(
	  phy = simPhy,
	  traits = simChar,
	  intrinsicFn=brownianIntrinsic,
	  extrinsicFn=nullExtrinsic,
	  startingPriorsFns="normal",
	  startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
	  intrinsicPriorsFns=c("exponential"),
	  intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
	  extrinsicPriorsFns=c("fixed"),
	  extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
	  generation.time=10000000,
	  standardDevFactor=0.2,
	  plot=FALSE,
	  StartSims=10,
	  epsilonProportion=0.7,
	  epsilonMultiplier=0.7,
	  nStepsPRC=3,
	  numParticles=3,
	  jobName="exampleRun",
	  stopRule=FALSE,
	  multicore=FALSE,
	  coreLimit=1
	  )
	  )
	expect_is(results, "list")
	expect_false(any(is.na(pairwiseESS(results$particleDataFrame))))
})



test_that("doRun_rej works", {
  data(simRunExample)
	set.seed(1)
   expect_warning(
   resultsRej <- doRun_rej(phy = simPhy, traits = simChar,
    intrinsicFn = brownianIntrinsic, 
	extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal", 
	startingPriorsValues = matrix(c(mean(simChar[,1]), sd(simChar[, 1]))), 
	  intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,
      byrow = FALSE), extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,
      byrow = FALSE), StartSims = 10, jobName = "examplerun_rej",
    abcTolerance = 0.05, multicore = FALSE, coreLimit = 1)
	)
	expect_is(resultsRej, "list")
	expect_false(any(is.na(pairwiseESS(resultsRej$particleDataFrame))))
})

