test_that("plotPriorPost works", {
  data(simRunExample)
  plotPrior(priorFn = "exponential", priorVariables = c(10))
  
  plotPrior(priorFn = "normal", priorVariables = c(1,
    2))
  
  plotPrior(priorFn = "gamma", priorVariables = c(2,
    0.2), plotQuants = FALSE, plotLegend = FALSE)
  priorKernal <- getUnivariatePriorCurve(priorFn = "normal",
    priorVariables = c(28, 2), nPoints = 1e+05, from = NULL,
    to = NULL, prob = 0.95)
  postKernal <- getUnivariatePosteriorCurve(acceptedValues = results$particleDataFrame$StartingStates1,
    from = NULL, to = NULL, prob = 0.95)
  priorKernal
  postKernal
  plotUnivariatePosteriorVsPrior(posteriorCurve = postKernal,
    priorCurve = priorKernal, label = "parameter",
    trueValue = NULL, prob = 0.95)
})
