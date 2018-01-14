test_that("doSimulation works", {
  set.seed(1)
  tree <- rcoal(5)
  tree$edge.length <- tree$edge.length * 20
  
  char <- doSimulation(
	phy = tree, 
	intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, 
	startingValues = c(10),
    intrinsicValues = c(0.01), 
	extrinsicValues = c(0),
	generation.time=1000000,
	saveHistory = FALSE
	)
	expect_equal(class(char[,1]), "integer")
	expect_equal(dim(char)[1], Ntip(tree))
	
  char <- doSimulation(
	phy = tree, 
	intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.05,10, 0.01), 
	extrinsicValues = c(0, 0.1, 0.25),
	generation.time=1000000,
	saveHistory = FALSE
	)
	expect_equal(class(char[,1]), "integer")
	expect_equal(dim(char)[1], Ntip(tree))
	
})

test_that("doSimulationForPlotting works", {
  set.seed(1)
  tree <- rcoal(5)
  tree$edge.length <- tree$edge.length * 20
	
  char <- doSimulationForPlotting(
	phy = tree, 
	intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, 
	startingValues = c(10),
    intrinsicValues = c(0.01), 
	extrinsicValues = c(0),
	generation.time=1000000,
    plot = FALSE, 
	saveHistory = FALSE
	)
	expect_equal(class(char[,1]), "integer")
	expect_equal(dim(char)[1], Ntip(tree))
	
  char <- doSimulationForPlotting(
	phy = tree, 
	intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.05,10, 0.01), 
	extrinsicValues = c(0, 0.1, 0.25),
	generation.time=1000000,
    plot = TRUE, 
	saveHistory = FALSE
	)
	expect_equal(class(char[,1]), "integer")
	expect_equal(dim(char)[1], Ntip(tree))
})

	
test_that("doSimulationWithPossibleExtinction works", {
  set.seed(1)
  tree <- rcoal(5)
  tree$edge.length <- tree$edge.length * 20
	
	
  expect_warning(
  charDoSim <- doSimulationWithPossibleExtinction(
	phy = tree,
    intrinsicFn = brownianIntrinsic, 
	extrinsicFn = nullExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.01),
	generation.time=1000000,
    extrinsicValues = c(0), 
	saveHistory = FALSE
	)
	)
	
	expect_equal(class(charDoSim[,1]), "integer")
	expect_equal(dim(charDoSim)[1], Ntip(tree))
})

