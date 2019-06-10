#  Calculate frequency of finding best solution in Geiger
#
#  This function is taken from an internal Geiger function
# probably really REALLY not needed now that geiger is no longer imported
#
#  @param fitContResult An object returned by geiger::fitContinuous()

#  @param tol Tolerance for equality of solutions

#  @return The frequency with which the best solution was found

# @name solnfreq
# @rdname solnfreq
# @export
solnfreq <- function(fitContResult, tol = .Machine$double.eps^0.5){
    ll = logLik(fitContResult)
    aa = abs(fitContResult$res[, "lnL"]-ll) <= tol
    result <- max(1, sum(aa[!is.na(aa)]))/length(aa)
    return(result)
    }




#  Calculate summary statistics
#
#  This function creates a vector of summary statistics for TreEvo analysis
#
#  Calculates 17 summary statistics from the geiger::fitContinuous.hacked function
#  (brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, 
#  lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, 
#  ou.alpha, ou.aic, white.lnl, white.aic), plus raw.mean, raw.max, raw.min, 
#  raw.var, raw.median, and all tip character values, phylogenetic independent
#  contrasts, ancestral state reconstruction values, and the range of ancestral
#  state reconstruction confidence interval.
#

# actually the above isn't true - isn't clear it uses geiger::fitContinuous.hacked at all!

# @inheritParams doSimulation
# @inheritParams doRun

#  @param niter.brown Number of random starts for BM model (min of 2)

#  @param niter.lambda Number of random starts for lambda model (min of 2)

#  @param niter.delta Number of random starts for delta model (min of 2)

#  @param niter.OU Number of random starts for OU model (min of 2)

#  @param niter.white Number of random starts for white model (min of 2)

#  @param do.CI Use confidence interval? By default, only for ultrametric trees

#  @return Returns a vector of summary statistics

#  @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

# @keywords summaryStatsLong

#  @examples
#
#  #summaryStatsLong(phy, char)
#

# @name summaryStatsLong
# @rdname summaryStatsLong
# @export
summaryStatsLong <- function(phy, traits
        #, niter.brown = 25, niter.lambda = 25, niter.delta = 25, niter.OU = 25, niter.white = 25, 
        #, do.CI = is.ultrametric(phy)
        ) {
    #
    if (any(phy$edge.length == 0)){
        if(!any(phy$edge[which(phy$edge.length == 0), 2] %in% phy$edge[, 1])){
        #if(any(phy$edge.length == 0)){
            stop(paste0(
				"There are zero-length terminal branches in your trees\n",
				"   Summary statistics cannot be calculated on such trees, please rescale edges."
				))
        }
    }
	#
	
	
	matchTreeWithTrait <- function(tree, traitsOld){
		# first, convert trait from a matrix or data.frame, if it is such
			# refine down to a single trait
		if(inherits(traitsOld, "matrix", "data.frame")) {
			my.names <- rownames(traitsOld)
			trait <- as.numeric(traitsOld[, 1])
			names(trait) <- my.names
			}
		# okay now it should be a vector
		# 
		# is it named
		if(is.null(names(trait))){
			stop("input trait data is not labeled with taxon unit names")
			}
		# check length
			# test that tree and trait has same number of tips / same number of taxa
		if(Ntip(tree) != length(trait)){
			warning(
				"Number of trait values input is not equal to the number of tip taxa used"
				)
			
			
			}
		
		
		return(trait)
		}

	
	


    #    if(is.null(names(traits)))
    #        names(traits) <- rownames(traits)
    #    traits <- as.data.frame(traits)
	#
    #it actually runs faster without checking for cores. And we parallelize elsewhere
	#
    brown <- getBM(phy = phy, dat = traits)    #, niterN = niter.brown
    brown.lnl <- brown$lnl
    brown.beta  <- brown$beta
    brown.aic  <- brown$aic
    #
    lambda <- getLambda(phy = phy, dat = traits)    #, niterN = niter.lambda
    lambda.lnl <- lambda$lnl
    lambda.beta  <- lambda$beta
    lambda.aic  <- lambda$aic
    lambda.lambda  <- lambda$lambda
    #    
    delta <- getDelta(phy = phy, dat = traits)    #, niterN = niter.delta
    delta.lnl <- delta$lnl
    delta.beta  <- delta$beta
    delta.aic  <- delta$aic
    delta.delta  <- delta$delta
    #
    ou <- getOU(phy = phy, dat = traits)    #, niterN = niter.OU
    ou.lnl <- ou$lnl
    ou.beta  <- ou$beta
    ou.aic  <- ou$aic
    ou.alpha  <- ou$alpha
    #
    white <- getWhite(phy = phy, dat = traits)    #, niterN = niter.white
    white.lnl <- white$lnl
    white.aic  <- white$aic
	#
    raw.mean <- as.numeric(mean(traits))
    raw.max <- as.numeric(max(traits))
    raw.min <- as.numeric(min(traits))
    raw.var <- as.numeric(var(traits))
    raw.median <- as.numeric(median(traits))    #message("summaryStatsLong")
	#
    pic <- makeQuiet(as.vector(pic.ortho(as.matrix(traits), phy)))  #independent contrasts
	#
    #aceResults <- makeQuiet(ace(traits, phy))
    #anc.states <- as.vector(aceResults$ace) #ancestral states
    #if(do.CI) {
    #    anc.CIrange <- as.vector(aceResults$CI95[, 2]-aceResults$CI95[, 1]) #range between upper and lower 95% CI
    #    summarystats <- c(summarystats, anc.CIrange)
    #}
	#
    #
    #fastAnc is much faster than ape's ace for our purposes
    ancResults <- phytools::fastAnc(tree = phy, x = traits, CI = TRUE)
    anc.states <- ancResults$ace
    anc.CIrange <- ancResults$CI
    #
	#
    #combined summary stats
    summarystats <- c(
		brown.lnl, brown.beta, brown.aic, 
		lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, 
        delta.lnl, delta.beta, delta.delta, delta.aic, 
		ou.lnl, ou.beta, ou.alpha, ou.aic, 
		white.lnl, white.aic, 
        raw.mean, raw.max, raw.min, raw.var, raw.median, 
		traits[[1]], pic, anc.states, anc.CIrange
		)
	#
	#
    summarystats[which(is.finite(summarystats) == FALSE)] <- NA
    #
    #while(sink.number()>0) {sink()}
    summarystats
    }
