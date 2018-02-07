test_that("doSimulation works", {
  set.seed(1)
  tree <- rcoal(5)
  tree$edge.length <- tree$edge.length * 20

  charDoSim <- doSimulation(
	phy = tree,
	intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic,
	startingValues = c(10),
    intrinsicValues = c(0.01),
	extrinsicValues = c(0),
	generation.time=1000000)
	
	expect_equal(class(charDoSim[,1]), "numeric")
	expect_equal(dim(charDoSim)[1], Ntip(tree))
	
  charDoSim <- doSimulation(
	phy = tree,
	intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10),
	intrinsicValues = c(0.05,10, 0.01),
	extrinsicValues = c(0, 0.1, 0.25),
	generation.time=1000000)
	
	expect_equal(class(charDoSim[,1]), "numeric")
	expect_equal(dim(charDoSim)[1], Ntip(tree))
	
})

	
