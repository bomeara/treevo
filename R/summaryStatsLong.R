summaryStatsLong<-function(phy, data) {
	if (any(phy$edge.length==0)){
		if(!any(phy$edge[which(phy$edge.length==0),2] %in% phy$edge[,1])){
		#if(any(phy$edge.length==0)){
			return("There are zero branch lengths at the tips of your trees--will not run properly")
		}
	}	
	#sink(file="/dev/null")
	
	brown<-fitContinuous.hacked(phy=phy, data=data, model="BM")[[1]] #will only run if want to do brownian summary stats
	brown.lnl<-as.numeric(brown$lnl) #divide by zero so we get Inf if we don't want that summary stat
	brown.beta <-as.numeric(brown$beta)
	brown.aic <-as.numeric(brown$aic)

	lambda<-fitContinuous.hacked(phy=phy, data=data, model="lambda")[[1]]
	lambda.lnl <-as.numeric(lambda$lnl)
	lambda.beta <-as.numeric(lambda$beta)
	lambda.lambda <-as.numeric(lambda$lambda)
	lambda.aic <-as.numeric(lambda$aic)

	delta<-fitContinuous.hacked(phy=phy, data=data, model="delta")[[1]]
	delta.lnl <-as.numeric(delta$lnl)
	delta.beta <-as.numeric(delta$beta)
	delta.delta <-as.numeric(delta$delta)
	delta.aic <-as.numeric(delta$aic)
	
	ou<-fitContinuous.hacked(phy=phy, data=data, model="OU")[[1]]
	ou.lnl <-as.numeric(ou$lnl)
	ou.beta <-as.numeric(ou$beta)
	ou.alpha <-as.numeric(ou$alpha)
	ou.aic <-as.numeric(ou$aic)
	
	white<-fitContinuous.hacked(phy=phy, data=data, model="white")[[1]]
	white.lnl<-as.numeric(white$lnl)
	white.aic<-as.numeric(white$aic)
	
	
	raw.mean<-as.numeric(mean(data))
	raw.max<-as.numeric(max(data))
	raw.min<-as.numeric(min(data))
	raw.var<-as.numeric(var(data))
	raw.median<-as.numeric(median(data[,]))	#cat("summaryStatsLong")
	summarystats<-c(brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic, raw.mean, raw.max, raw.min, raw.var, raw.median, data[[1]] )

#add contrasts (2n-1)
#add ancestral states (2n-1)
#add variance on ancestral states (2n-1)


	summarystats[which(is.finite(summarystats)==FALSE)]<-NA
	
	while(sink.number()>0) {sink()}
	summarystats
}
