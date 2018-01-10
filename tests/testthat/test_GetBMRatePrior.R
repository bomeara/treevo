test_that("getBMRatePrior works", {
  data(simRunExample)
  xxx<-getBMRatePrior(phy = simPhy, traits = simChar, timeStep = 1)
})
