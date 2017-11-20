test_that("simulation ran", {
	data(simRun)
	charDoSim<-doSimulationWithPossibleExtinction(
		phy=simPhy,
		intrinsicFn=brownianIntrinsic,
		extrinsicFn=nullExtrinsic,
		startingValues=c(10), #root state
		intrinsicValues=c(0.01),
		extrinsicValues=c(0),
		timeStep=0.0001,
		saveHistory=FALSE)
	expect_equal(class(charDoSim[,1]), "numeric")
	expect_equal(dim(charDoSim)[1], 30)
	}
)


test_that("doPRC runs correctly", {
	data(simRun)
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
	  timeStep=0.0001,
	  standardDevFactor=0.2,
	  plot=FALSE,
	  StartSims=10,
	  epsilonProportion=0.7,
	  epsilonMultiplier=0.7,
	  nStepsPRC=3,
	  numParticles=20,
	  jobName="exampleRun",
	  stopRule=FALSE,
	  multicore=FALSE,
	  coreLimit=1
	)
	expect_is(results, "list")
	expect_false(any(is.na(PairwiseESS(results$particleDataFrame))))
})


test_that("plotting works", {
	data(simRun)
	expect_error(plotPosteriors(results$particleDataFrame, results$PriorMatrix), NA)
})
