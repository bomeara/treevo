#################
# 
#####################################################
# cd C:\Users\notDave\Desktop\remoteWorkplace && R CMD BATCH --no-save --no-restore treevo_ABC_vignette.R 
################

## ------------------------------------------------------------------------
# Control Box
devtools::install_github("bomeara/treevo",ref="vignette_third_try")
setwd("C:/Users/notDave/Desktop/remoteWorkplace")

multicore<-TRUE
coreLimit<-6
generation.time<-10000
nStepsPRC=3
numParticles=100


## ----echo=FALSE----------------------------------------------------------
library(paleotree)
data(RaiaCopesRule)
set.seed(444)


## ------------------------------------------------------------------------
plot(ladderize(ammoniteTreeRaia));axisPhylo()

## ------------------------------------------------------------------------
tree<-multi2di(ammoniteTreeRaia)
tree<-addTermBranchLength(tree,0.01)
#plot(ladderize(tree));axisPhylo()

## ---- echo=FALSE, fig.show='hold'----------------------------------------
#plotTraitgram(tree=tree, trait=sutureComplexity,
#	 conf.int=FALSE, main="Ammonite Suture Complexity")

#plotTraitgram(tree=tree, trait=shellSize,
#	 conf.int=FALSE, main="Ammonite Shell Diameter")

## ------------------------------------------------------------------------

library(TreEvo)

# character data for doRun must be in matrix form
  # with rows labeled with taxon names

charData<-matrix(sutureComplexity,ncol=1)
rownames(charData)<-names(sutureComplexity)

 save.image(file="treevo_vignette_workspace.Rdata")



## ------------------------------------------------------------------------

# typical Brownian Motion
resultsBM<-doRun_prc(
  phy = tree,
  traits = charData,
  intrinsicFn=brownianIntrinsic,
  extrinsicFn=nullExtrinsic,
  startingPriorsFns="normal",
  startingPriorsValues=matrix(c(mean(charData[,1]), sd(charData[,1]))),
  intrinsicPriorsFns=c("exponential"),
  intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
  extrinsicPriorsFns=c("fixed"),
  extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
  generation.time=generation.time,
  nRuns=3,
  nStepsPRC=nStepsPRC,
  numParticles=numParticles,
  jobName="typicalBMrun",
  stopRule=FALSE,
  multicore=multicore,
  coreLimit=coreLimit,
  verboseParticles=TRUE
  )

 save.image(file="treevo_vignette_workspace.Rdata")



## ------------------------------------------------------------------------

resultsBound<-doRun_prc(
  phy = tree,
  traits = charData,
  intrinsicFn=boundaryMinIntrinsic,
  extrinsicFn=nullExtrinsic,
  startingPriorsFns="normal",
  startingPriorsValues=matrix(c(mean(charData[,1]), sd(charData[,1]))),
  intrinsicPriorsFns=c("exponential","normal"),
  intrinsicPriorsValues=matrix(c(10, 10, -10, 1), nrow=2, byrow=FALSE),
  extrinsicPriorsFns=c("fixed"),
  extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
  generation.time=generation.time,
  nRuns=3,
  nStepsPRC=nStepsPRC,
  numParticles=numParticles,
  jobName="BMwithBoundRun",
  stopRule=FALSE,
  multicore=multicore,
  coreLimit=coreLimit,
  verboseParticles=FALSE
  )


 save.image(file="treevo_vignette_workspace.Rdata")
