test_that("bayesCoverageProb works", {
  set.seed(1)
  data(simRunExample)
  genPar <- c(ancState, genRate)
  HPDs <- list(results$HPD, resultsBound$HPD)
  result<-bayesCoverageProb(RealParam = genPar, HPD = HPDs, verbose = FALSE)
})

test_that("compareListIPD works", {
  set.seed(1)
  data(simRunExample)
  pdfList <- list(results$particleDataFrame, resultsBound$particleDataFrame)
  result<-compareListIPD(particleDataFrameList = pdfList, verbose = FALSE)
})

test_that("credibleInt works", {
  set.seed(1)
  data(simRunExample)
  result<-credibleInt(results$particleDataFrame)
})

test_that("getBMRatePrior works", {
  set.seed(1)
  data(simRunExample)
  result<-getBMRatePrior(phy = simPhy, traits = simChar, timeStep = 1)
})

test_that("highestPostDens works", {
  set.seed(1)
  data(simRunExample)
  result<-highestPostDens(results$particleDataFrame, percent = 0.95,
    returnData = FALSE)
})
