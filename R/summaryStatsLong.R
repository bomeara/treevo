#' Calculate frequency of finding best solution in Geiger
#'
#' This function is taken from an internal Geiger function
#'
#' @param x A returned object from fitContinuous()
#' @param tol Tolerance for equality of solutions
#' @return The frequency with which the best solution was found
solnfreq <- function(x, tol = .Machine$double.eps^0.5){
			ll=logLik(x)
			aa=abs(x$res[,"lnL"]-ll)<=tol
			max(1,sum(aa[!is.na(aa)]))/length(aa)
}




#' Calculate summary statistics
#'
#' This function creates a vector of summary statistics for TreEvo analysis
#'
#' Calculates 17 summary statistics from the fitContinuous.hacked function
#' (brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda,
#' lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta,
#' ou.alpha, ou.aic, white.lnl, white.aic), plus raw.mean, raw.max, raw.min,
#' raw.var, raw.median, and all tip character values, phylogenetic independent
#' contrasts, ancestral state reconstruction values, and the range of ancestral
#' state reconstruction confidence interval.
#'
#' @param phy Tree (Phylogenetic tree in phylo format)
#' @param traits data matrix with rownames equal to phy
#' @param niter.brown Number of random starts for BM model (min of 2)
#' @param niter.lambda Number of random starts for lambda model (min of 2)
#' @param niter.delta Number of random starts for delta model (min of 2)
#' @param niter.OU Number of random starts for OU model (min of 2)
#' @param niter.white Number of random starts for white model (min of 2)
#' @param do.CI Use confidence interval? By default, only for ultrametric trees
#' @return Returns a vector of summary statistics
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
#' @keywords summaryStatsLong
#' @examples
#'
#' #summaryStatsLong(phy, char)
#'
summaryStatsLong<-function(phy, traits, niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25, do.CI=is.ultrametric(phy)) {
	if (any(phy$edge.length==0)){
		if(!any(phy$edge[which(phy$edge.length==0),2] %in% phy$edge[,1])){
		#if(any(phy$edge.length==0)){
			return("There are zero branch lengths at the tips of your trees--will not run properly")
		}
	}
	if(class(traits)=="matrix" | class(traits)=="data.frame") {
		my.names <- rownames(traits)
		traits <- as.numeric(traits[,1])
		names(traits) <- my.names
	}

#	if(is.null(names(traits)))
#		names(traits) <- rownames(traits)
#	traits<-as.data.frame(traits)
	brown<-fitContinuous(phy=phy, dat=traits, model="BM", ncores=1, control=list(niter=niter.brown)) #it actually runs faster without checking for cores. And we parallelize elsewhere
	brown.lnl<-as.numeric(brown$opt$lnL)
	brown.beta <-as.numeric(brown$opt$sigsq)
	brown.aic <-as.numeric(brown$opt$aic)

	lambda<-fitContinuous(phy=phy, dat=traits, model="lambda", ncores=1, control=list(niter=niter.lambda))
	lambda.lnl <-as.numeric(lambda$opt$lnL)
	lambda.beta <-as.numeric(lambda$opt$sigsq)
	lambda.lambda <-as.numeric(lambda$opt$lambda)
	lambda.aic <-as.numeric(lambda$opt$aic)

	delta<-fitContinuous(phy=phy, dat=traits, model="delta", ncores=1, control=list(niter=niter.delta))
	delta.lnl <-as.numeric(delta$opt$lnL)
	delta.beta <-as.numeric(delta$opt$sigsq)
	delta.delta <-as.numeric(delta$opt$delta)
	delta.aic <-as.numeric(delta$opt$aic)

	ou<-fitContinuous(phy=phy, dat=traits, model="OU", ncores=1, control=list(niter=niter.OU))
	ou.lnl <-as.numeric(ou$opt$lnL)
	ou.beta <-as.numeric(ou$opt$sigsq)
	ou.alpha <-as.numeric(ou$opt$alpha)
	ou.aic <-as.numeric(ou$opt$aic)

	white<-fitContinuous(phy=phy, dat=traits, model="white", ncores=1, control=list(niter=niter.white))
	white.lnl<-as.numeric(white$opt$lnL)
	white.aic<-as.numeric(white$opt$aic)


	raw.mean<-as.numeric(mean(traits))
	raw.max<-as.numeric(max(traits))
	raw.min<-as.numeric(min(traits))
	raw.var<-as.numeric(var(traits))
	raw.median<-as.numeric(median(traits))	#cat("summaryStatsLong")

	pic<-as.vector(pic.ortho(as.matrix(traits), phy))  #independent contrasts
	aceResults<-ace(traits, phy)
	anc.states<-as.vector(aceResults$ace) #ancestral states

	#combined summary stats
	summarystats<-c(brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic, raw.mean, raw.max, raw.min, raw.var, raw.median, traits[[1]], pic, anc.states)

	if(do.CI) {
		anc.CIrange<-as.vector(aceResults$CI95[,2]-aceResults$CI95[,1]) #range between upper and lower 95% CI
		summarystats <- c(summarystats, anc.CIrange)
	}


	summarystats[which(is.finite(summarystats)==FALSE)]<-NA

	while(sink.number()>0) {sink()}
	summarystats
}
