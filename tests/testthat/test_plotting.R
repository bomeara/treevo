test_that("plotABC_3D works", {
	set.seed(1)
  if (requireNamespace("gpclib", quietly = TRUE) & requireNamespace("rgl", quietly = TRUE)) {
    data(simRunExample)
    plotABC_3D(
		particleDataFrame = results$particleDataFrame,
		parameter = 7,
		show.particles = "none",
		plot.parent = FALSE,
		realParam = FALSE,
		realParamValues = NA)
  }
})

test_that("plotPosteriors works", {
	set.seed(1)
	data(simRunExample)
	plotPosteriors(
		particleDataFrame = results$particleDataFrame,
		priorsMat = results$PriorMatrix,
		realParam = TRUE,
		realParamValues = c(ancState, genRate)
		)
})

test_that("plotPosteriors correctly gives an error", {
	set.seed(1)
	data(simRunExample)
	expect_error(
		plotPosteriors(
			results$particleDataFrame,
			results$PriorMatrix)
		,NA)
})

test_that("plotPriorPost works", {
	set.seed(1)
  data(simRunExample)
  plotPrior(priorFn = "exponential", priorVariables = c(10))

  plotPrior(priorFn = "normal", priorVariables = c(1,
    2))

  plotPrior(priorFn = "gamma", priorVariables = c(2,
    0.2), plotQuants = FALSE, plotLegend = FALSE)
  priorKernal <- getUnivariatePriorCurve(priorFn = "normal",
    priorVariables = c(28, 2), nPoints = 1e+05, from = NULL,
    to = NULL, prob = 0.95)
  postKernal <- getUnivariatePosteriorCurve(acceptedValues = results$particleDataFrame$starting_1,
    from = NULL, to = NULL, prob = 0.95)
  priorKernal
  postKernal
  plotUnivariatePosteriorVsPrior(posteriorCurve = postKernal,
    priorCurve = priorKernal, label = "parameter",
    trueValue = NULL, prob = 0.95)
})

test_that("parentOffspringPlots works", {
	set.seed(1)
  data(simRunExample)
  parentOffspringPlots(results$particleDataFrame)
})
