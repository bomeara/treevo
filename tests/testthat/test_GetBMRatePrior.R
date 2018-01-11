test_that("getBMRatePrior works", {
  data(simRunExample)
  result<-getBMRatePrior(phy = simPhy, traits = simChar, timeStep = 1)
})
