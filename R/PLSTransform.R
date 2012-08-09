PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,ncol=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  print(paste("dim of summaryValuesMatrix is ",dim(summaryValuesMatrix)))
  predict(pls.model,summaryValuesMatrix,type="scores")
}
