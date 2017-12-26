test_that("extrinsicModels works", {
  tree <- rcoal(30)
  simPhy$edge.length <- simPhy$edge.length * 20
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = nearestNeighborDisplacementExtrinsic,
    startingValues = c(10), intrinsicValues = c(0),
    extrinsicValues = c(0.1, 0.1, 0.1), generation.time = 1e+05,
    plot = TRUE, saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = everyoneDisplacementExtrinsic, startingValues = c(10),
    intrinsicValues = c(0), extrinsicValues = c(0.1,
      0.1, 0.1), generation.time = 1e+05, plot = TRUE,
    saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), intrinsicValues = c(0),
    extrinsicValues = c(0.1, 0.1, 2), generation.time = 1e+05,
    plot = TRUE, saveHistory = FALSE)
})
