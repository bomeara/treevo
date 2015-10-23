test_that("doPRC runs correctly", {
	data(simData)
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
	  TreeYears=1000,
	  standardDevFactor=0.2,
	  plot=FALSE,
	  StartSims=10,
	  epsilonProportion=0.7,
	  epsilonMultiplier=0.7,
	  nStepsPRC=5,
	  numParticles=100,
	  jobName="exampleRun",
	  stopRule=FALSE,
	  multicore=FALSE,
	  coreLimit=1
	)
	expect_is(results, "list")	
})
