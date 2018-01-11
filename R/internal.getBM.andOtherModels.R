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


########################################################
		
getBM<-function(phy,dat,niterN,niter.goal=NA){
	fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat, 
		model="BM", ncores=1, control=list(niter=niterN)))
	res$lnl<-as.numeric(fit$opt$lnL)
	res$beta <-as.numeric(fit$opt$sigsq)
	res$aic <-as.numeric(fit$opt$aic)
	if(!is.na(niter.goal)){
		res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
		}
	return(res)
	}
	
	
getLambda<-function(phy,dat,niterN,niter.goal=NA){
	fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat, 
		model="lambda", ncores=1, control=list(niter=niterN)))
	res$lnl<-as.numeric(fit$opt$lnL)
	res$beta <-as.numeric(fit$opt$sigsq)
	res$lambda<-as.numeric(fit$opt$lambda)
	res$aic <-as.numeric(fit$opt$aic)
	if(!is.na(niter.goal)){
		res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
		}
	return(res)
	}		

getDelta<-function(phy,dat,niterN,niter.goal=NA){
	fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
		model="delta", ncores=1, control=list(niter=niterN)))
	res$lnl<-as.numeric(fit$opt$lnL)
	res$beta <-as.numeric(fit$opt$sigsq)
	res$delta<-as.numeric(fit$opt$delta)
	res$aic <-as.numeric(fit$opt$aic)
	if(!is.na(niter.goal)){
		res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
		}
	return(res)
	}			

getOU<-function(phy,dat,niterN,niter.goal=NA){
	fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
		model="OU", ncores=1, control=list(niter=niterN)))
	res$lnl<-as.numeric(fit$opt$lnL)
	res$beta <-as.numeric(fit$opt$sigsq)
	res$alpha<-as.numeric(fit$opt$alpha)
	res$aic <-as.numeric(fit$opt$aic)
	if(!is.na(niter.goal)){
		res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
		}
	return(res)
	}		

getWhite<-function(phy,dat,niterN,niter.goal=NA){
	fit<-makeQuiet(geiger::fitContinuous(phy=phy, dat=dat,
		model="white", ncores=1, control=list(niter=niterN)))
	res$lnl<-as.numeric(fit$opt$lnL)
	res$beta <-as.numeric(fit$opt$sigsq)
	res$alpha<-as.numeric(fit$opt$alpha)
	res$aic <-as.numeric(fit$opt$aic)
	if(!is.na(niter.goal)){
		res$niter.g <- round(max(10, min(niter.goal/solnfreq(fit),100)))
		}
	return(res)
	}			
	