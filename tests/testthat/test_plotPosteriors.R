test_that("plotPosteriors works", {
  data(simRunExample)
  plotPosteriors(particleDataFrame = results$particleDataFrame,
    priorsMat = results$PriorMatrix, realParam = TRUE,
    realParamValues = c(ancState, genRate))
})

test_that("plotting correctly gives an error", {
	data(simRunExample)
	set.seed(1)
	expect_error(plotPosteriors(results$particleDataFrame, results$PriorMatrix), NA)
})