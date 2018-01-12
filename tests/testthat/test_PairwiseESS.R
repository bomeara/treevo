test_that("pairwiseESS works", {
  data(simRunExample)
  expect_warning(
  res<-pairwiseESS(results$particleDataFrame)
  )
})
