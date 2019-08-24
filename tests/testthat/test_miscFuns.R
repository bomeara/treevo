#test_that("bayesCoverageProb works", {
#  set.seed(1)
#  data(simRunExample)
#  genPar <- c(ancStateExample, genRateExample)
#  HPDs <- list(
#	resultsBMExample[[1]]$HPD, 
#	resultsBoundExample[[1]]$HPD
#	)
#  result<-bayesCoverageProb(RealParam = genParExample, HPD = HPDs, verbose = FALSE)
#})

test_that("compareListIPD works", {
  set.seed(1)
  data(simRunExample)
  pdfList <- list(
	resultsBMExample[[1]]$particleDataFrame, 
	resultsBoundExample[[1]]$particleDataFrame
	)
  result<-compareListIPD(
	particleDataFrameList = pdfList, 
	verbose = FALSE
	)
})

#test_that("credibleInt works", {
#  set.seed(1)
#  data(simRunExample)
#  result<-credibleInt(resultsBMExample[[1]]$particleDataFrame)
#})

test_that("getBMRatePrior works", {
  set.seed(1)
  data(simRunExample)
  result<-getBMRatePrior(phy = simPhyExample, 
	traits = simCharExample, timeStep = 1
	)
})

#test_that("highestPostDensity works", {
#  set.seed(1)
#  data(simRunExample)
#  result<-highestPostDensity(resultsBMExample[[1]]$particleDataFrame, percent = 0.95,
#    returnData = FALSE)
#})
