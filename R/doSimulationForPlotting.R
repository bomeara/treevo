#' Discrete Time Character Simulation
#' 
#' This function evolves continuous characters in a discrete time process
#' 
#' When saveHistory is TRUE, processor time will increase quite a bit.
#' SaveRealParams is useful for tracking the "real" true values if simulating
#' data for abc runs.  It is not useful for empirical abc runs.
#' 

# @param splits Output from the function getSimulationSplits; is a data frame
# of branching times, ancestor and descendant vectors

#' @param phy A phylogenetic tree, in package \code{ape}'s \code{phylo} format.

#' @param intrinsicFn Name of intrinsic function characters should be simulated
#' under

#' @param extrinsicFn Name of extrinsic function characters should be simulated
#' under

#' @param startingValues State at the root

#' @param intrinsicValues Vector of values corresponding to the params of the
#' intrinsic model

#' @param extrinsicValues Vector of values corresponding to the params of the
#' extrinsic model

#' @param timeStep This value corresponds to the number of discrete time steps
#' on the shortest branch

#' @param saveHistory Saves the character history throughout the simulation

#' @param saveRealParams Saves intrinsicValues and extrinsicValues as both a
#' matrix and a vector

#' @param jobName Optional name for the job

#' @return A data frame of species character (tip) values in the tree.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords doSimulation doSimulationForPlotting


#' @examples
#' 
#' 
#' tree<-rcoal(30)
#' 


#' @name doSimulationForPlotting
#' @rdname doSimulationForPlotting
#' @export
