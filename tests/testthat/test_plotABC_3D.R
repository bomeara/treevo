test_that("plotABC_3D works", {
  if (requireNamespace("gpclib", quietly = TRUE) & requireNamespace("rgl",
    quietly = TRUE)) {
    data(simRunExample)
    plotABC_3D(particleDataFrame = results$particleDataFrame,
      parameter = 7, show.particles = "none", plot.parent = FALSE,
      realParam = FALSE, realParamValues = NA)
  }
})
