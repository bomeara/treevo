#' Example Analysis Output of a Simulated Dataset
#' 
#' Simulated 30-taxon coalescent tree and simulated character data from a Brownian Motion
#' intrinsic model (\code{\link{brownianIntrinsic}}), with saved generating parameters.
#' Character data was generated under this model using \code{doSimulation}.
#' Also includes results from an example analysis.
#' 
#' 
#' @name simRunExample

#' @aliases simRunExample simPhy simChar results resultsBound genRate ancState 
#simCharOut

#' @docType data

#' @format
#' Loading the \code{simRunExample} example dataset adds seven new objects to the namespace:
#' 
#' \describe{

#' \item{\code{simPhy}}{A simulated 30-tip coalescent phylogeny in typical \code{phylo} format.}

#' \item{\code{ancState}}{The starting ancestral value,
#' used for generating the simulated continuous trait data.}

#' \item{\code{genRate}}{The true rate of trait change under Brownian Motion,
#' used for generating the simulated continuous trait data.}

# \item{\code{simCharOut}}{The raw output of \code{\link{doSimulation}} on \code{simPhy},
# under the model \code{\link{brownianIntrinsic}}.}

#' \item{\code{simChar}}{The output of \code{\link{doSimulation}} on \code{simPhy}, 
#'  under the model \code{\link{brownianIntrinsic}}. composed of just the simulated
#' trait values as a one-column matrix with row names indicating tip labels, as desired by \code{doRun} functions.}

#' \item{\code{results}}{The results of \code{\link{doRun_prc}}, under the generating model of \code{\link{brownianIntrinsic}}}

#' \item{\code{resultsBound}}{The results of \code{\link{doRun_prc}}, under the incorrect model of \code{\link{boundaryMinIntrinsic}}}

#' }
#' 
#' The objects \code{results} and \code{resultsBound} are lists composed of a number
#' of elements (see the documentation for the \code{\link{doRun_prc}} function for more detail). These elements are
#' respectively \code{input.data}, \code{PriorMatrix}, \code{particleDataFrame}, \code{toleranceVector}, \code{phy},
#' \code{traits}, \code{simTime}, \code{time.per.gen}, \code{credibleInt}, and \code{HPD}.
#' 


#' @examples
#' 
#' data(simRunExample)
#' 
#' # ...things to do with this data?
#' 
#' ###################
#' 
#' # This data was generated using this process:
#' 
#' \dontrun{
#' 
#' library(TreEvo)
#' 
#' set.seed(1)
#' simPhy<-rcoal(20)
#' # get realistic edge lengths
#' simPhy$edge.length<-simPhy$edge.length*20
#' 
#' # plot with time axis (root is about ~15 Ma)
#' plot(simPhy)
#' axisPhylo()
#' 
#' genRate<-c(0.001)
#' ancState<-c(10)
#' 
#' #Simple Brownian motion
#' simChar<-doSimulation(
#' 	phy=simPhy,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=ancState, #root state
#' 	intrinsicValues=genRate,
#' 	extrinsicValues=c(0),
#' 	generation.time=10000
#' 	)

# 
# # clean for use with doRun
# simChar<-simCharOut[,"statesmatrix",drop=FALSE]
# rownames(simChar)<-simPhy$tip.label[simCharOut$taxonid]

#' 
#' results<-doRun_prc(
#' 	phy = simPhy,
#' 	traits = simChar,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingPriorsFns="normal",
#' 	startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#' 	intrinsicPriorsFns=c("exponential"),
#' 	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#' 	extrinsicPriorsFns=c("fixed"),
#' 	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' 	generation.time=10000,
#' 	nRuns=2,
#' 	nStepsPRC=3,
#' 	numParticles=20,
#' 	jobName="examplerun_prc",
#' 	stopRule=FALSE,
#' 	multicore=FALSE,
#' 	verboseParticles=TRUE,
#' 	coreLimit=1
#' 	)
#' 
#' resultsBound<-doRun_prc(
#' 	phy = simPhy,
#' 	traits = simChar,
#' 	intrinsicFn=boundaryMinIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingPriorsFns="normal",
#' 	startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#' 	intrinsicPriorsFns=c("exponential","normal"),
#' 	intrinsicPriorsValues=matrix(c(10, 10, -10, 1), nrow=2, byrow=FALSE),
#' 	extrinsicPriorsFns=c("fixed"),
#' 	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' 	generation.time=10000,
#' 	nRuns=2
#' 	nStepsPRC=3,
#' 	numParticles=20,
#' 	jobName="examplerun_prc_bound",
#' 	stopRule=FALSE,
#' 	multicore=FALSE,
#' 	verboseParticles=TRUE,
#' 	coreLimit=1
#' 	)
#' 
#' rm(.Random.seed)
#' save.image(file="simRunExample.rdata")
#' 
#' 
#' }
#' 
#' 

# @keywords datasets
NULL

