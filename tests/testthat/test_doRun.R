test_that("doRun works", {
  data(simRunExample)
  resultsPRC <- doRun_prc(phy = simPhy, traits = simChar,
    intrinsicFn = brownianIntrinsic, extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal", startingPriorsValues = matrix(c(mean(simChar[,
      1]), sd(simChar[, 1]))), intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,
      byrow = FALSE), extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,
      byrow = FALSE), generation.time = 1e+05, standardDevFactor = 0.2,
    plot = FALSE, StartSims = 10, epsilonProportion = 0.7,
    epsilonMultiplier = 0.7, nStepsPRC = 3, numParticles = 20,
    jobName = "examplerun_prc", stopRule = FALSE, multicore = FALSE,
    coreLimit = 1)
  resultsPRC
  resultsRej <- doRun_rej(phy = simPhy, traits = simChar,
    intrinsicFn = brownianIntrinsic, extrinsicFn = nullExtrinsic,
    startingPriorsFns = "normal", startingPriorsValues = matrix(c(mean(simChar[,
      1]), sd(simChar[, 1]))), intrinsicPriorsFns = c("exponential"),
    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2,
      byrow = FALSE), extrinsicPriorsFns = c("fixed"),
    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2,
      byrow = FALSE), StartSims = 10, jobName = "examplerun_rej",
    abcTolerance = 0.05, multicore = FALSE, coreLimit = 1)
  resultsRej
})
