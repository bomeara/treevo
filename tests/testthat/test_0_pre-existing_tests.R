test_that("simulation ran", {
	data(simRunExample)
	set.seed(1)
	charDoSim<-doSimulationWithPossibleExtinction(
		phy=simPhy,
		intrinsicFn=brownianIntrinsic,
		extrinsicFn=nullExtrinsic,
		startingValues=c(10), #root state
		intrinsicValues=c(0.01),
		extrinsicValues=c(0),
		generation.time=100000,
		saveHistory=FALSE)
	expect_equal(class(charDoSim[,1]), "numeric")
	expect_equal(dim(charDoSim)[1], 30)
	}
)




test_that("plotting works", {
	data(simRunExample)
	set.seed(1)
	expect_error(plotPosteriors(results$particleDataFrame, results$PriorMatrix), NA)
})
