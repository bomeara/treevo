test_that("compareListIPD works", {
  data(simRunExample)
  pdfList <- list(results$particleDataFrame, resultsBound$particleDataFrame)
  compareListIPD(particleDataFrameList = pdfList, verbose = TRUE)
})
