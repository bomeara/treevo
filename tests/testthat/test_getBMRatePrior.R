test_that("getBMRatePrior works", {
  data(simRunExample)
  getBMRatePrior(phy = simPhy, traits = simChar, timeStep = 1)
})
