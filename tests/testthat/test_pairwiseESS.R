test_that("pairwiseESS works", {
  data(simRunExample)
  pairwiseESS(results$particleDataFrame)
})
