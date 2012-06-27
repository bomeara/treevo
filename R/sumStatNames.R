sumStatNames<-function(ntips=Ntip(phy)) {
  return(c("BML", "BMbeta", "BMaic", "lambdaL", "lambdabeta", "lambdalambda", "lambdaaic", "deltaL", "deltabeta", "deltadelta", "deltaaic", "ouL", "oubeta", "oualpha", "ouaic", "whiteL", "whiteaic", "rawmean", "rawmax", "rawmin", "rawvar", "rawmedian", paste("tipVal", sequence(ntips), sep=""), paste("pic", sequence(ntips-1), sep=""), paste("ancState", sequence(ntips-1), sep=""), paste("ancCIrange", sequence(ntips-1), sep="")))
}
