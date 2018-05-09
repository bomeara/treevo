test_that("doRun_prc runs correctly in parallel", {
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
	  generation.time=1000000,
	  nRuns=2,
	  startSimsPerParam=5,
	  nStepsPRC=2,
	  numParticles=5,
	  jobName="exampleRun",
	  diagnosticPRCmode=TRUE,
	  verboseParticles=TRUE,
	  stopRule=FALSE,
	  multicore=TRUE,
	  coreLimit=2
	  )
	  )
	expect_output(print(results))
	expect_is(results, "list")
	expect_warning(expect_false(any(is.na(pairwiseESS(results[[1]]$particleDataFrame)))))
	expect_false(any(is.na(pairwiseESS(list(results[[1]]$particleDataFrame,results[[2]]$particleDataFrame)))))
	expect_false(any(is.na(pairwiseESS(results))))
})



