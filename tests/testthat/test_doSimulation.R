test_that("doSimulation works", {
  tree <- rcoal(30)
  tree$edge.length <- tree$edge.length * 20
  char <- doSimulation(
  phy = tree, 
	intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, 
	startingValues = c(10),
    intrinsicValues = c(0.01), 
	extrinsicValues = c(0),
	generation.time=100000,
	saveHistory = FALSE)
  char <- doSimulation(
  phy = tree, 
	intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.05,
      10, 0.01), 
	extrinsicValues = c(0, 0.1, 0.25),
	generation.time=100000,
	saveHistory = FALSE)
  char <- doSimulationForPlotting(
	phy = tree, 
	intrinsicFn = brownianIntrinsic,
    extrinsicFn = nullExtrinsic, 
	startingValues = c(10),
    intrinsicValues = c(0.01), 
	extrinsicValues = c(0),
	generation.time=100000,
    plot = FALSE, 
	saveHistory = FALSE)
  char <- doSimulationForPlotting(
	phy = tree, 
	intrinsicFn = boundaryMinIntrinsic,
    extrinsicFn = ExponentiallyDecayingPushExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.05,10, 0.01), 
	extrinsicValues = c(0, 0.1, 0.25),
	generation.time=100000,
    plot = TRUE, 
	saveHistory = FALSE)
  char <- doSimulationWithPossibleExtinction(
	phy = tree,
    intrinsicFn = brownianIntrinsic, 
	extrinsicFn = nullExtrinsic,
    startingValues = c(10), 
	intrinsicValues = c(0.01),
	generation.time=100000,
    extrinsicValues = c(0), 
	saveHistory = FALSE)
})

test_that("simulation ran", {
	data(simRunExample)
	set.seed(1)
	charDoSim<-doSimulationWithPossibleExtinction(
		phy=simPhy,
		intrinsicFn=brownianIntrinsic,
		extrinsicFn=nullExtrinsic,
		startingValues=c(10), #root state
		intrinsicValues=c(0.01),
		extrinsicValues=c(0),
		generation.time=100000,
		saveHistory=FALSE)
	expect_equal(class(charDoSim[,1]), "numeric")
	expect_equal(dim(charDoSim)[1], 30)
	}
)