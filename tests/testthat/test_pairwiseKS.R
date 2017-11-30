test_that("pairwiseKS works", {
  data(simRunExample)
  pdfList <- list(Brownian = results$particleDataFrame,
    Bounded = resultsBound$particleDataFrame)
  pairwiseKS(particleDataFrameList = pdfList)
})
