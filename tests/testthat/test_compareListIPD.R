test_that("compareListIPD works", {
  data(simRunExample)
  pdfList <- list(results$particleDataFrame, resultsBound$particleDataFrame)
  xxx<-compareListIPD(particleDataFrameList = pdfList, verbose = FALSE)
})
