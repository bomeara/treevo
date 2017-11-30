test_that("intrinsicModels works", {
  tree <- rcoal(30)
  simPhy$edge.length <- simPhy$edge.length * 20
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01), extrinsicValues = c(0),
    timeStep = 1e-04, plot = TRUE, saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = boundaryIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01, 0, 15), extrinsicValues = c(0),
    timeStep = 1e-04, plot = TRUE, saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = minBoundaryAutoregressiveIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01, 3, 0.1, 0), extrinsicValues = c(0),
    timeStep = 1e-04, plot = TRUE, saveHistory = FALSE)
})
