ntax<-NTAXPLUGIN
intrinsicFn<-INTRINSICPLUGIN
extrinsicFn<-EXTRINSICPLUGIN


phy<-rcoal(ntax)
phy$edge.length<-phy$edge.length/max(branching.times(phy))
splits<-getSimulationSplits(phy)
char<-convertTaxonFrameToGeigerData (doSimulation(splits=splits,intrinsicFn= intrinsicFn,extrinsicFn= extrinsicFn,startingStates=c(3),intrinsicValues=.06,extrinsicValues=0,timeStep=0.001),phy)


a<-doRun(phy=phy,traits=char,intrinsicFn= intrinsicFn,extrinsicFn= extrinsicFn,startingMatrix=matrix(data=c(0,15),nrow=2),intrinsicMatrix=matrix(data=c(0.0001,.1),nrow=2),extrinsicMatrix=matrix(data=c(0,0),nrow=2),timeStep=0.001,standardDevFactor=0.05, plot=FALSE,nrepSim=10,startingStatesGuess=c(2),intrinsicValuesGuess=c(0.06),extrinsicValuesGuess=c(0),epsilonProportion=0.2,epsilonMultiplier=0.5,nStepsPRC=3,numParticles=10)
print(a)
presentABCOutput(a)