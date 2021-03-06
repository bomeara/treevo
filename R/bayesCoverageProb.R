#* Bayesian Coverage Probability
#* 
#* This function calculates coverage probability
#* for a list of highest posterior densities
#* (HPDs), calculated as Highest Density Regions (HDRs)
#* via the function \code{highestDensityInterval}, for a set of parameters.
#* 
#* Only for use with simulated data to test models.
#* 

#* @param RealParam Real parameter values.

#* @param HPD list of highest posterior density from \code{\link{doRun}} results.

#* @param verbose If \code{TRUE}, commented screen output is produced.

#* @return Returns a value for each free parameter that describes the
#* percentage that the real value falls within the HPD.

#* @seealso 
#* This function simply wraps \code{\link{highestDensityInterval}} for a multivariate dataset.

#* @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished


#* @examples
#* 
#* data(simRunExample)
#* 
#* # real (generating) parameters
#* genPar <- c(ancStateExample, genRateExample)
#* 
#* HPDs <- list(resultsBMExample[[1]]$HPD, 
#*    resultsBoundExample[[1]]$HPD)
#* 
#* bayesCoverageProb(
#*    RealParam = genPar, HPD = HPDs, verbose = TRUE
#*    )
#* 


#used to be BCP = Bayesian Coverage Probability, now bayesCoverageProb (08-29-17)
#RealParam can be "RealParams$vector" from doSimulation c(x1, x2, ...) or a list
#HPD should be a list of HPD output from different runs
#Calculates what percent of the time the real parameter falls into the HPD

# 07-03-18
# This functions seems problematic -
	# it basically asks about coverage for each parameter as if they were indep of each other
# And reports the percentage of parameters that fall within the already-calculated HPD
# that is fundamentally troublesome (see testMultivarOutlierHDR)
# It'd be better to ask if a particular particle was
	# within the multivariate data cloud or not 
	# again testMultivarOutlierHDR
# but that's just TRUE or FALSE, no percentage

#####

# We can use `bayesCoverageProb` to compare true / generating parameter values to the posteriors of analyses done on that data. This function calculates what percent of the time the real parameter falls into the HPD. However, this function might be problematic, as it basically asks about coverage for each parameter as if they were indep of each other, and reports the percentage of parameters that fall within the already-calculated HPD. This is fundamentally troublesome, as explained in `testMultivarOutlierHDR`. It'd be better to ask if a particular particle was within the multivariate data cloud or not. See `testMultivarOutlierHDR` for more discussion.

#####

#* @name bayesCoverageProb
#* @rdname bayesCoverageProb
#* @export
bayesCoverageProb <- function(RealParam, HPD, verbose = FALSE){
    if(inherits(RealParam, "numeric")){
        rps <- vector("list", length = length(HPD))
        for (i in 1: length(HPD)){
            rps[[i]] <- RealParam
            }
        }else{
            rps <- RealParam
            }
    # should we allow for only one HPD to be evaluated?
    #If(is.data.frame(HPD)){
    #    HPD <- list(HPD)
    #    }
    #
    #if(length(RealParam)  !=  dim(HPD[[1]])[1]){ warning("RealParams and HPD do not match")}
        #need something like this, but it will have to be after changing the RealParam to take a list
    Covered <- matrix(nrow = length(HPD), ncol = length(rps[[i]]))
    colnames(Covered) <- rownames(HPD[[1]])
    for(i in 1:length(HPD)){
        for(j in 1:length(rps[[i]])){
            if(!is.na(HPD[[i]][j, 4])) { #keep only parameters that vary (ie, not fixed)
                if(HPD[[i]][j, 4]  >=  rps[[i]][j] && rps[[i]][j]  >=  HPD[[i]][j, 3]) {
                    Covered[i, j] <- 1
                }
                else{
                    Covered[i, j] <- 0
                }
            }
        }
    }
    CoverProb <- apply(Covered, 2, mean)
    if(verbose){
        CoverProb <- vector("list")
        CoverProb$byRun <- Covered
        CoverProb$BCP <- apply(Covered, 2, mean)
    }
    return(CoverProb)
}
