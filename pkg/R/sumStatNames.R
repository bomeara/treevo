sumStatNames<-function(ntips=Ntip(phy)) {
  return(c("BML", "BMbeta", "BMaic", "lambdaL", "lambdabeta", "lambdalambda", "lambdaaic", "deltaL", "deltabeta", "deltadelta", "deltaaic", "ouL", "oubeta", "oualpha", "ouaic", "whiteL", "whiteaic", "rawmean", "rawmax", "rawmin", "rawvar", "rawmedian", paste("tipVal", 1:ntips, sep=""), paste("pic", 1:ntips-1, sep=""), paste("ancState", 1:ntips-1, sep=""), paste("ancCIrange", 1:ntips-1, sep="")))
}
