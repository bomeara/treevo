#' Example Analysis Output of a Simulated Dataset
#' 
#' Simulated 30-taxon coalescent tree and simulated character data from a Brownian Motion
#' intrinsic model (\code{\link{brownianIntrinsic}}), with saved generating parameters.
#' Character data was generated under this model using \code{doSimulation}.
#' Also includes results from an example analysis.
#' 
#' 
#' @name simRunExample

#' @aliases 
#' simRunExample 
#' simPhyExample simCharExample
#' resultsBMExample resultsBoundExample
#' genRateExample ancStateExample

#simCharOut

#' @docType data

#' @format
#' Loading the \code{simRunExample} example dataset adds seven new 
#' objects to the namespace:
#' 
#' \describe{

#' \item{\code{simPhyExample}}{A simulated 30-tip coalescent
#' phylogeny in typical \code{phylo} format.}

#' \item{\code{ancStateExample}}{The starting ancestral value, 
#' used for generating the simulated continuous trait data.}

#' \item{\code{genRateExample}}{The true rate of trait change under Brownian Motion, 
#' used for generating the simulated continuous trait data.}

# \item{\code{simCharOutExample}}{The raw output 
# of \code{\link{doSimulation}} on \code{simPhyExample}, 
# under the model \code{\link{brownianIntrinsic}}.}

#' \item{\code{simCharExample}}{The output of \code{\link{doSimulation}} on \code{simPhyExample}, 
#'  under the model \code{\link{brownianIntrinsic}}. composed of just the simulated
#' trait values as a one-column matrix with row names indicating tip 
#' labels, as desired by \code{doRun} functions.}

#' \item{\code{resultsBMExample}}{The results of \code{\link{doRun_prc}},
#'  under the generating model of \code{\link{brownianIntrinsic}}}

#' \item{\code{resultsBoundExample}}{The results of 
#' \code{\link{doRun_prc}}, under the incorrect model
#'  of \code{\link{boundaryMinIntrinsic}}}

#' }
#' 
#' The objects \code{resultsBMExample} and \code{resultsBoundExample} 
#' are lists composed of a number of elements (see the documentation 
#' for the \code{\link{doRun_prc}} function for more detail).


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
#' simPhyExample <- rcoal(20)
#' # get realistic edge lengths
#' simPhyExample$edge.length <- simPhyExample$edge.length*20
#' 
#' # plot with time axis (root is about ~15 Ma)
#' plot(simPhyExample)
#' axisPhylo()
#' 
#' genRateExample <- c(0.001)
#' ancStateExample <- c(10)
#' 
#' #Simple Brownian motion
#' simCharExample <- doSimulation(
#'     phy = simPhyExample, 
#'     intrinsicFn = brownianIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingValues = ancStateExample, #root state
#'     intrinsicValues = genRateExample, 
#'     extrinsicValues = c(0), 
#'     generation.time = 10000
#'     )
#' 
#' resultsBMExample <- doRun_prc(
#'     phy = simPhyExample, 
#'     traits = simCharExample, 
#'     intrinsicFn = brownianIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingPriorsFns = "normal", 
#'     startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1]))), 
#'     intrinsicPriorsFns = c("exponential"), 
#'     intrinsicPriorsValues = list(10), 
#'     extrinsicPriorsFns = c("fixed"), 
#'     extrinsicPriorsValues = list(0), 
#'     generation.time = 10000, 
#'     nRuns = 2, 
#'     nStepsPRC = 3, 
#'     numParticles = 20, 
#'     nInitialSimsPerParam = 10, 
#'     jobName = "examplerun_prc", 
#'     stopRule = FALSE, 
#'     multicore = FALSE, 
#'     verboseParticles = TRUE, 
#'     coreLimit = 1
#'     )
#' 
#' resultsBoundExample <- doRun_prc(
#'     phy = simPhyExample, 
#'     traits = simCharExample, 
#'     intrinsicFn = boundaryMinIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingPriorsFns = "normal", 
#'     startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1]))), 
#'     intrinsicPriorsFns = c("exponential", "normal"), 
#'     intrinsicPriorsValues = list(10,c(-10,1)), 
#'     extrinsicPriorsFns = c("fixed"), 
#'     extrinsicPriorsValues = list(0), 
#'     generation.time = 10000, 
#'     nRuns = 2, 
#'     nStepsPRC = 3, 
#'     numParticles = 20, 
#'     nInitialSimsPerParam = 10, 
#'     jobName = "examplerun_prc_bound", 
#'     stopRule = FALSE, 
#'     multicore = FALSE, 
#'     verboseParticles = TRUE, 
#'     coreLimit = 1
#'     )
#' 
#' rm(.Random.seed)
#' save.image(file = "simRunExample.rdata")
#' if(interactive()){savehistory("simRunExample.Rhistory")}
#' 
#' }
#' 
#' 

# @keywords datasets
NULL

