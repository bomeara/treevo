#' Calculate summary statistics
#' 
#' This function creates a vector of summary statistics for TreEvo analysis
#' 
#' Caluclates 17 summary statistics from the fitContinuous.hacked function
#' (brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda,
#' lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta,
#' ou.alpha, ou.aic, white.lnl, white.aic), plus raw.mean, raw.max, raw.min,
#' raw.var, raw.median, and all tip character values, phylogenetic independent
#' contrasts, ancestral state reconstruction values, and the range of ancestral
#' state reconstruction confidence interval.
#' 
#' @param phy Tree (Phylogenetic tree in phylo format)
#' @param traits data matrix with rownames equal to phy
#' @return Returns a vector of summary statistics
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
#' @keywords summaryStatsLong
#' @examples
#' 
#' #summaryStatsLong(phy, char)
#' 
summaryStatsLong<-function(phy, traits) {
	if (any(phy$edge.length==0)){
		if(!any(phy$edge[which(phy$edge.length==0),2] %in% phy$edge[,1])){
		#if(any(phy$edge.length==0)){
			return("There are zero branch lengths at the tips of your trees--will not run properly")
		}
	}	
	if(is.null(names(traits)))
		names(traits) <- rownames(traits)
	tratis<-as.data.frame(traits)
	brown<-fitContinuous(phy=phy, dat=traits, model="BM") 
	brown.lnl<-as.numeric(brown$opt$lnL) 
	brown.beta <-as.numeric(brown$opt$sigsq)
	brown.aic <-as.numeric(brown$opt$aic)

	lambda<-fitContinuous(phy=phy, dat=traits, model="lambda")
	lambda.lnl <-as.numeric(lambda$opt$lnL)
	lambda.beta <-as.numeric(lambda$opt$sigsq)
	lambda.lambda <-as.numeric(lambda$opt$lambda)
	lambda.aic <-as.numeric(lambda$opt$aic)

	delta<-fitContinuous(phy=phy, dat=traits, model="delta")
	delta.lnl <-as.numeric(delta$opt$lnL)
	delta.beta <-as.numeric(delta$opt$sigsq)
	delta.delta <-as.numeric(delta$opt$delta)
	delta.aic <-as.numeric(delta$opt$aic)
	
	ou<-fitContinuous(phy=phy, dat=traits, model="OU")
	ou.lnl <-as.numeric(ou$opt$lnL)
	ou.beta <-as.numeric(ou$opt$sigsq)
	ou.alpha <-as.numeric(ou$opt$alpha)
	ou.aic <-as.numeric(ou$opt$aic)
	
	white<-fitContinuous(phy=phy, dat=traits, model="white")
	white.lnl<-as.numeric(white$opt$lnL)
	white.aic<-as.numeric(white$opt$aic)
	
	
	raw.mean<-as.numeric(mean(traits))
	raw.max<-as.numeric(max(traits))
	raw.min<-as.numeric(min(traits))
	raw.var<-as.numeric(var(traits))
	raw.median<-as.numeric(median(traits[,]))	#cat("summaryStatsLong")
	
	pic<-as.vector(pic.ortho(as.matrix(traits), phy))  #independent contrasts
	aceResults<-ace(as.matrix(traits), phy)
	anc.states<-as.vector(aceResults$ace) #ancestral states 
	anc.CIrange<-as.vector(aceResults$CI95[,2]-aceResults$CI95[,1]) #range between upper and lower 95% CI
	
	#combined summary stats	
	summarystats<-c(brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic, raw.mean, raw.max, raw.min, raw.var, raw.median, traits[[1]], pic, anc.states, anc.CIrange)


	summarystats[which(is.finite(summarystats)==FALSE)]<-NA
	
	while(sink.number()>0) {sink()}
	summarystats
}
