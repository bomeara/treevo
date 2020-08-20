test_that("doRun_prc runs correctly in parallel", {
    #skip_on_cran()
    
	data(simRunExample)
    set.seed(1)
	
    expect_warning(
		results <- doRun_prc(
			phy = simPhyExample,
			traits = simCharExample,
			intrinsicFn=brownianIntrinsic,
			extrinsicFn=nullExtrinsic,
			startingPriorsFns="normal",
			startingPriorsValues = list(
				c(mean(simCharExample[, 1]), 
					sd(simCharExample[, 1]))
				), 
			intrinsicPriorsFns = c("exponential"),
			intrinsicPriorsValues = list(10),
			extrinsicPriorsFns = c("fixed"),
			extrinsicPriorsValues = list(0),
			generation.time=1000000,
			nRuns=2,
			nInitialSimsPerParam=5,
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
	
    expect_warning(
		expect_false(any(
			is.na(pairwiseESS(results[[1]]$particleDataFrame))
			))
		)
    
	expect_false(
		any(is.na(
			pairwiseESS(list(
				results[[1]]$particleDataFrame,
				results[[2]]$particleDataFrame
				))
			))
		)
    
	expect_false(
		any(is.na(pairwiseESS(results)))
		)
	
})



