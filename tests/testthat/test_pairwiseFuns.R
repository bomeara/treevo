test_that("pairwiseKS works", {
  data(simRunExample)
  pdfList <- list(Brownian = results$particleDataFrame,
    Bounded = resultsBound$particleDataFrame)
  pairwiseKS(particleDataFrameList = pdfList)
})

test_that("pairwiseESS works", {
  data(simRunExample)
  expect_warning(
  res<-pairwiseESS(results$particleDataFrame)
  )
})
