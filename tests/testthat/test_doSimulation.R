test_that("doSimulation works", {
  tree <- rcoal(30)
  simPhy$edge.length <- simPhy$edge.length * 20
  char <- doSimulation(phy = tree, intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01), extrinsicValues = c(0),
    saveHistory = FALSE)
  char <- doSimulation(phy = tree, intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), intrinsicValues = c(0.05,
      10, 0.01), extrinsicValues = c(0, 0.1, 0.25),
    saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01), extrinsicValues = c(0),
    plot = FALSE, saveHistory = FALSE)
  char <- doSimulationForPlotting(phy = tree, intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), intrinsicValues = c(0.05,
      10, 0.01), extrinsicValues = c(0, 0.1, 0.25),
    plot = TRUE, saveHistory = FALSE)
  char <- doSimulationWithPossibleExtinction(phy = tree,
    intrinsicFn = brownianIntrinsic, extrinsicFn = nullExtrinsic,
    startingValues = c(10), intrinsicValues = c(0.01),
    extrinsicValues = c(0), saveHistory = FALSE)
})
