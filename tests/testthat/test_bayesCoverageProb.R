test_that("bayesCoverageProb works", {
  data(simRunExample)
  genPar <- c(ancState, genRate)
  HPDs <- list(results$HPD, resultsBound$HPD)
  result<-bayesCoverageProb(RealParam = genPar, HPD = HPDs, verbose = FALSE)
})
