test_that("plotPosteriors works", {
  data(simRunExample)
  plotPosteriors(particleDataFrame = results$particleDataFrame,
    priorsMat = results$PriorMatrix, realParam = TRUE,
    realParamValues = c(ancState, genRate))
})
