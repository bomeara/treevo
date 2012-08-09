PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,nrow=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  predict(pls.model,summaryValuesMatrix,type="scores")
}
