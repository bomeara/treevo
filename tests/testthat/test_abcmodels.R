test_that("ABC intrinsicModels works", {
  set.seed(1)
  tree <- rcoal(5)
  simPhy$edge.length <- simPhy$edge.length * 20
  #
  expect_warning(
  char <- doSimulation(phy = tree, intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01), extrinsicValues = c(0),
    generation.time = 100000)
    )
    expect_warning(
  char <- doSimulation(phy = tree, intrinsicFn = boundaryIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01, 0, 15), extrinsicValues = c(0),
    generation.time = 100000)
    )
    expect_warning(
  char <- doSimulation(phy = tree, intrinsicFn = minBoundaryAutoregressiveIntrinsic,
    extrinsicFn = nullExtrinsic, startingValues = c(10),
    intrinsicValues = c(0.01, 3, 0.1, 0), extrinsicValues = c(0),
    generation.time = 100000)
    )
})

test_that("ABC extrinsicModels works", {
  set.seed(1)
  tree <- rcoal(5)
  simPhy$edge.length <- simPhy$edge.length * 20
  #
  
    expect_warning(
  char <- doSimulation(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = nearestNeighborDisplacementExtrinsic,
    startingValues = c(10), intrinsicValues = c(0),
    extrinsicValues = c(0.1, 0.1, 0.1), generation.time = 100000)
    )
    expect_warning(
  char <- doSimulation(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = everyoneDisplacementExtrinsic, startingValues = c(10),
    intrinsicValues = c(0), extrinsicValues = c(0.1,
      0.1, 0.1), generation.time = 100000)
    )
    expect_warning(  
  char <- doSimulation(phy = tree, intrinsicFn = nullIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), intrinsicValues = c(0),
    extrinsicValues = c(0.1, 0.1, 2), generation.time = 100000)
    )
})
