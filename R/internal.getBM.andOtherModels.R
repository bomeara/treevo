# internal functions for fitting models of trait evolution using ML functions in other packages for use as summary statistics
# 
# 01-11-18: now centralizing all calls to geiger::fitContinuous so that we can replace with usage from phylolm::phylolm

	
# models available in phylolm	
# "BM", "OUrandomRoot", "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend"
#
# "OUfixedRoot" gives same answer as geiger::fitContinuous("OU") for both ultrametric and non-ultrametric trees
#
# will have to build a star tree and fit BM to do White Noise with phylolm
#

#White noise with phylolm, compared to fitCotinuous:
#
#tree<-rtree(100)
#trait<-rTraitCont(tree)
#
#fitCont_WN<-geiger::fitContinuous(phy=tree,dat=trait,model="white")
#fitCont_WN
#
#star<-stree(Ntip(tree))
#star$edge.length<-rep(1,nrow(star$edge))
#phylolm::phylolm(trait~1,phy=star,model="BM")
#


########################################################
		
getBM<-function(phy,dat){			#,niterN,niter.goal=NA
	res<-list()
	#fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat, 
	#	model="BM", ncores=1, control=list(niter=niterN)))
	#res$lnl<-as.numeric(fit$logLik)
	#res$beta <-as.numeric(fit$opt$sigsq)
	#res$aic <-as.numeric(fit$aic)

	
	fit<-phylolm::phylolm(formula=dat~1,phy=phy,model="BM")
	res$lnl<-as.numeric(fit$logLik)
	res$beta <-as.numeric(fit$sigma2)
	res$aic <-as.numeric(fit$aic)

	
	#if(!is.na(niter.goal)){
	#	res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
	#	}
	return(res)
	}
	
	
getLambda<-function(phy,dat){			#,niterN,niter.goal=NA
	res<-list()
	#fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat, 
	#	model="lambda", ncores=1, control=list(niter=niterN)))
	fit<-phylolm::phylolm(formula=dat~1,phy=phy,model="lambda")
	res$lnl<-as.numeric(fit$logLik)
	res$beta <-as.numeric(fit$sigma2)
	res$lambda<-as.numeric(fit$optpar)
	res$aic <-as.numeric(fit$aic)
	#if(!is.na(niter.goal)){
	#	res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
	#	}
	return(res)
	}		

getDelta<-function(phy,dat){			#,niterN,niter.goal=NA
	res<-list()
	#fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
	#	model="delta", ncores=1, control=list(niter=niterN)))
	fit<-phylolm::phylolm(formula=dat~1,phy=phy,model="delta")
	res$lnl<-as.numeric(fit$logLik)
	res$beta <-as.numeric(fit$sigma2)
	res$delta<-as.numeric(fit$optpar)
	res$aic <-as.numeric(fit$aic)
	#if(!is.na(niter.goal)){
	#	res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
	#	}
	return(res)
	}			

getOU<-function(phy,dat){			#,niterN,niter.goal=NA
	res<-list()
	#fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
	#	model="OU", ncores=1, control=list(niter=niterN)))
	fit<-phylolm::phylolm(formula=dat~1,phy=phy,model="OUfixedRoot")
	res$lnl<-as.numeric(fit$logLik)
	res$beta <-as.numeric(fit$sigma2)
	res$alpha<-as.numeric(fit$optpar)
	res$aic <-as.numeric(fit$aic)
	#if(!is.na(niter.goal)){
	#	res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
	#	}
	return(res)
	}		

getWhite<-function(phy,dat){		#,niterN,niter.goal=NA
	res<-list()
	#fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
	#	model="white", ncores=1, control=list(niter=niterN)))
	#
	# will have to build a star tree and fit BM to do White Noise with phylolm
	#
	star<-stree(Ntip(tree))
	star$edge.length<-rep(1,nrow(star$edge))
	fit<-phylolm::phylolm(formula=trait~1,phy=star,model="BM")	
	#
	res$lnl<-as.numeric(fit$logLik)
	res$beta <-as.numeric(fit$sigma2)
	res$aic <-as.numeric(fit$aic)
	#if(!is.na(niter.goal)){
	#	res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
	#	}
	return(res)
	}			
	