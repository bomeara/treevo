




test_that("plotting works", {
	data(simRunExample)
	set.seed(1)
	expect_error(plotPosteriors(results$particleDataFrame, results$PriorMatrix), NA)
})
