test_that("highestPostDens works", {
  data(simRunExample)
  highestPostDens(results$particleDataFrame, percent = 0.95,
    returnData = FALSE)
})
