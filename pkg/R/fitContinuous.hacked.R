library(corpcor)


`fitContinuous.hacked` <-
function(phy, data, data.names=NULL, model=c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend"), bounds=NULL,  meserr=NULL, userstart=NULL, badLnL=10000)
{
	#print(badLnL)
	# sort is T because sub-functions assume data are in
	# this particular order
	
	model<-match.arg(model)
	
	td<-treedata(phy, data, data.names, sort=T)

	ntax=length(td$phy$tip.label)
	
	
	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me	
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o<-match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr<-as.matrix(meserr[o,])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me<-td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o<-match(rownames(td$data), rownames(meserr))
			meserr=meserr[o,]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
				stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")	
		}
	}

	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	ds			<- list()
   		ds$tree 		<- td$phy          # TIP data 
   		ds$badLnL <- badLnL
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    cat("Fitting ", model, "model:\n")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00000001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 2.999999, 0.0000001, 50, -3, 0, 0.0000000001, 100, -100, 100), nrow=8, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "a", "nv", "mu");
    colnames(bounds.default) <- c("min", "max")

 	#---- USER DEFINED PARAMETER BOUNDS
 	if (is.null(bounds)) {
 		bounds <- bounds.default       # USE DEFAULTS
 	}else{
 		if (class(bounds)!="list"){
 			stop("Please specify user defined parameter bounds as a list()")
 		}else{
 			specified   <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
 							  is.null(bounds$kappa), is.null(bounds$delta),  is.null(bounds$alpha), is.null(bounds$a),
 							  is.null(bounds$nv), is.null(bounds$mu)
 							  )
 			bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$a, bounds$nv, bounds$mu), 
 								  nrow=sum(specified), ncol=2, byrow=TRUE
 								  )
 			rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "a",  "nv", "mu")[specified]
   	 		colnames(bounds.user) <- c("min", "max")
  
   	 		#----  SET FINAL SEARCH BOUNDS
 			bounds <- bounds.default
 			bounds[specified,] <- bounds.user     # Final Bounds
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))

  	ds$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    
    ds$userstart<-userstart
    
    result<-list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[,i]
    	ds$meserr=meserr[,i]
  		result[[i]]<-fitContinuousModel(ds, print=print)
  		if(!is.null(colnames(td$data))) names(result)[i]<-colnames(td$data)[i] else names(result)[i]<-paste("Trait", i, sep="")

  	}
  	result
}


`fitContinuousModel` <-
function(ds, print=TRUE)
{
	bounds 	<- ds$bounds
	model 	<- ds$model
	badLnL	<- ds$badLnL
	n 		<- length(ds$data)
	userstart <- ds$userstart

	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	
	beta.start<-var(ds$data)/max(branching.times(ds$tree))


	out         <- NULL
	
	y			<- ds$data			# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	n			<- length(y)
	

	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (model=="BM") {
		k<-2
	
		vcv<-vcv.phylo(tree)

		start=log(beta.start)
		#lower=log(bounds[1,"beta"])
		#upper=log(bounds[2,"beta"])
		
		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#print ("in tryFoo")
			#badLnL=100000
			result<-try(foo(x), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				#print("in badLnL")
				return(badLnL)
			}
		options(warn=0) 
		}
		
		foo<-function(x) {
			vv<-exp(x)*vcv
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		
		o<-nlm(tryFoo, p=start)
		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")

		results<-list(lnl=-o$minimum, beta=exp(o$estimate[1]))		
		#results<-list(lnl=-o$value, beta= exp(o$par))

	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model=="lambda"){
		k<-3
		
		start=log(c(beta.start, 0.5))
		#lower=log(bounds[1,c("beta","lambda")])
		#upper=log(bounds[2,c("beta","lambda")])
		
		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#print ("in tryFoo")
			#badLnL=100000
			result<-try(foo(x), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				#print(paste("in is.finite, result =", result))
				return(result)
			}
			else {
				#print("in badLnL")
				return(badLnL)
			}
		options(warn=0) 
		}
			
		foo<-function(x) {

			#print("in foo")
			vcv<-vcv.phylo(tree) #change back to tree
			index			<-	matrix(TRUE, n,n)
			diag(index)		<- FALSE
			#print("vcv is ")
			#print(vcv)
			#print(paste("x[1] = ",x[1],"x[2]=",x[2]))
			vcv[index] 	<- vcv[index]*exp(x[2])
			
			
			vv<-exp(x[1])*vcv
			#print("vv before meserr")
			#print(vv)
			
			diag(vv)<-diag(vv)+meserr^2
			#print("vv after meserr")
			#print(vv)

			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			#print("about to call dmvnormPseudoinverse")
			-dmvnormPseudoinverse(y, mu, vv, log=T)
		}
		#print("starting nlm")
		o<-nlm(tryFoo, p=start)
		#print("ending nlm")
		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$minimum, beta=exp(o$estimate[1]), lambda=exp(o$estimate[2]))	
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))

	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model=="kappa"){
		k<-3
		start=log(c(beta.start, 0.5))
		#lower=log(bounds[1,c("beta","kappa")])
		#upper=log(bounds[2,c("beta","kappa")])
				
		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#badLnL=100000
			result<-try(foo(x), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				return(badLnL)
			}
		options(warn=0) 
		}
		
		foo<-function(x) {

			t<-kappaTree(tree, kappa=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnormPseudoinverse(y, mu, vv, log=T)
		}

		o<-nlm(tryFoo, p=start)
		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$minimum, beta= exp(o$estimate[1]), kappa=exp(o$estimate[2]))	#used to say lambda?		
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))


	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model=="delta"){
		
		k<-3
				start=log(c(beta.start, 0.5))

	if(!is.null(userstart)) {
		start<-log(userstart)
	
	}
		#lower=log(bounds[1,c("beta","delta")])
		#upper=log(bounds[2,c("beta","delta")])
		
		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#badLnL=100000
			result<-try(foo(x), silent=T)
			#print(paste("badLnL in delta is set to",badLnL))
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				#print(paste("in is.finite, result =", result))
				return(result)
			}
			else {
				return(badLnL)
				print("in badLnL")
			}
		options(warn=0) 
		}
		
		foo<-function(x) {
			t<-deltaTree(tree, delta=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			#determinantVCV<-det(vv)
			badLnL=100000
			#if (determinantVCV==0){ #old delta had bounds, so couldn't get very low values and so didn't get singular matrices. Now that can happen, so have to guard against it
			#	warning("Possibly singular variance-covariance matrix, so giving this particular parameter combination a very bad likelihood score (rather than crashing)")
			#	return(badLnL)
			#}
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		#o<-optim(foo, p=start, lower=lower, upper=upper, method="SANN")
		o<-nlm(tryFoo, p=start)
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), delta=exp(o$par[2]))	
		results<-list(lnl=-o$minimum, beta= exp(o$estimate[1]), delta=exp(o$estimate[2]))	
	#----------------------------------
	#-----        WHITE NOISE     -----
	#----------------------------------	
	} else if (model=="white"){
		
		k<-2
		start=c(mean(y), log(var(y)))
		#lower=c(-Inf, log(bounds[1,"nv"]))
		#upper=c(Inf, log(bounds[2, "nv"]))
		
		tryLnl.noise<-function (p, x, se) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#badLnL=100000
			result<-try(Lnl.noise(p, x, se), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				return(badLnL)
			}
		}
		
		Lnl.noise<- function (p, x, se)
		# p is the vector of parameters, tree is not needed
		# x and se are trait means and std errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree (also optimum)
  			vs<- exp(p[2])		# white noise variance 
  			n<- length(x)
  			VV<- diag(vs, nrow=n)	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
    
  			## logl
  			M<- rep(root,n)
  			-dmvnorm(x, M, VV, log=TRUE)
		}

		o<- nlm(tryLnl.noise, p=start, x=y, se=meserr)		
		#o<- optim(start, fn=lnl.noise, x=y, se=meserr, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$minimum, beta=exp(o$estimate[1]), nv=exp(o$estimate[2]))	
		#results<-list(lnl=-o$value, mean= o$par[1], nv=exp(o$par[2]))	
	#----------------------------------
	#-----        TREND           -----
	#----------------------------------	
	} else if (model=="trend"){
		
		k<-3
		vcv<- vcv.phylo(tree)
  		ww<- lm(y ~ diag(vcv))
  		p0<- c(phylogMean(vcv, y), var(y)/max(branching.times(tree)), coef(ww)[2])
		if(is.na(p0[3])) {
			p0[3]<-0
			if(is.ultrametric(tree))
				cat("WARNING: Cannot estimate a trend with an ultrametric tree; lnl will be the same as the BM model")
		}
		#lower=c(-Inf, log(bounds[1,"beta"]), bounds[1,"mu"])
		#upper=c(Inf, log(bounds[2,"beta"]), bounds[2,"mu"])
		
		tryLnl.BMtrend <-function(p, vcv, x, se) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#badLnL=100000
			result<-try(Lnl.BMtrend(p, vcv, x, se), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				return(badLnL)
			}
		}
		
		Lnl.BMtrend<- function(p, vcv, x, se)
		# p is vector of parameters, tr is tree
		# x and se are vectors of trait means and standard errors
		{
  			## prep parameters
  			root<- p[1]	# trait value of root of tree
  			vs<- exp(p[2])		# BM variance 
  			ms<- p[3]		# BM trend
  			VV<- vs*vcv	
  			diag(VV)<- diag(VV) + se^2	# add sampling variance to diagonal
  			n<- length(x)
  
 			## logl
  			M<- root+ ms*diag(vcv)
  			- dmvnorm(x, M, VV, log=TRUE)
  		}

		o<- nlm(tryLnl.BMtrend, p=p0, vcv=vcv, x=y, se=meserr)
		#o<- optim(p0, fn=lnl.BMtrend, vcv=vcv, x=y, se=meserr, lower=lower, upper=upper, method="L")
		
		names(o$par)<-NULL
		
		results<-list(lnl=-o$minimum, mean=exp(o$estimate[1]), beta=exp(o$estimate[2]), mu=exp(o$estimate[3]))			
		#results<-list(lnl=-o$value, mean= o$par[1], beta=exp(o$par[2]), mu=o$par[3])		
	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model=="OU"){
	## modified 12 dec 07 to call ouMatrix(x) instead of vcv.phylo(ouTree(x))

		k<-3
		nlm.print.level<-0
		failureCountSwitch<-50
		start=log(c(beta.start, 0.5))
		#lower=log(bounds[1,c("beta","alpha")])
		#upper=log(bounds[2,c("beta","alpha")])
	
		vcvOrig<-vcv.phylo(tree)

		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#print ("in tryFoo")
			#badLnL=100000
			result<-try(foo(x), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				#print("in badLnL")
				return(badLnL)
			}
		options(warn=0) 
		}
		
		foo<-function(x) {
			vcv <- ouMatrix(vcvOrig, exp(x[2]))
			## t<-ouTree(tree, exp(x[2]))
			##vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnormPseudoinverse(y, mu, vv, log=T)
		}
		
		outTries<-list()
		
		# First one: try near BM solution
		start=c(log(beta.start), -50)
		outTries[[1]]<-nlm(foo, p=start)
		#outTries[[1]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		# Second one: try one with very strong constraints
		tv<-var(y)
		start=log(c(tv*2000, 1000))
		outTries[[2]]<-nlm(tryFoo, p=start,print.level=nlm.print.level)
		#outTries[[2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L")
	

	
		# Try ten random ones
		for(i in 1:10){
			failureCount<-0
			while(1) {
				lower=c(runif(2, min=-20, max=-1))
				upper=lower+runif(2, min=0, max=10)
				start=c(runif(1, min=lower[1], max=upper[1]), runif(1, min=lower[2], max=upper[2]))
				te<-try(outTries[[i+2]]<-nlm(foo, p=start,print.level=nlm.print.level,interlim=25), silent=T)
				#print(te)
				#te<-try(outTries[[i+2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				failureCount<-failureCount+1
				if (failureCount>failureCountSwitch) { #nlm isn't working, let's try optim
					#toptim<-try(outTries[[i+2]]<-nlm(tryFoo, p=start), silent=T)
					toptim<-try(outTries[[i+2]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
					if (class(toptim)!="try-error") {
						outTries[[i+2]]$estimate<-outTries[[i+2]]$par
						outTries[[i+2]]$minimum<-outTries[[i+2]]$value
						break
					}
				}
				}
				
		}
		
		# Try range of alphas
		atry<- -5:4
		stry<- log(tv*2*exp(atry))
		for(i in 1:10){
			failureCount<-0
			while(1) {
				lower=c(-20, -20)
				upper=c(10, 10)
				start=c(stry[i], atry[i])
				te<-try(outTries[[i+12]]<-nlm(foo, p=start,print.level=nlm.print.level,interlim=25), silent=T)
				#print(te)
				#te<-try(outTries[[i+12]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
				if(class(te)!="try-error") break
				failureCount<-failureCount+1
				if (failureCount>failureCountSwitch) { #nlm isn't working, let's try optim
					#toptim<-try(outTries[[i+12]]<-nlm(tryFoo, p=start), silent=T)
					toptim<-try(outTries[[i+12]]<-optim(foo, p=start, lower=lower, upper=upper, method="L"), silent=T)
					if (class(toptim)!="try-error") {
						outTries[[i+12]]$estimate<-outTries[[i+12]]$par
						outTries[[i+12]]$minimum<-outTries[[i+12]]$value
						break
					}
				}

				}
				
		}
		
		
		
		ntries<-22
		ltry<-numeric(ntries)
		lsol<-matrix(nrow= ntries, ncol=2)
		for(j in 1:ntries) {
				ltry[j]<-outTries[[j]]$minimum
				lsol[j,]<-exp(outTries[[j]]$estimate)
			}

		ltd<-ltry-min(ltry)
		b<-min(which(ltry==min(ltry)))

		gc<-which(ltd<0.01)
		us<-lsol[gc,1]
		usc<-sum((us-min(us))>0.01)			
		out<-outTries[[b[1]]]	
		print(out)
	
		if(usc>1) {out$message="Warning: likelihood surface is flat."}
			

		results<-list(lnl=-out$minimum, beta=exp(out$estimate[1]), alpha=exp(out$estimate[2]))	
		#results<-list(lnl=-out$value, beta= exp(out$par[1]), alpha=exp(out$par[2]), convergence=out$convergence, message=out$message, k=k)


	#----------------------------------
	#-----        EB ONLY      -----
	#----------------------------------	
	} else if(model=="EB"){

		k<-3
		start=c(log(beta.start), 0.01)
		#lower=c(log(bounds[1,"beta"]),bounds[1,"a"])
		#upper=c(log(bounds[2,"beta"]),bounds[2,"a"])
		
		tryFoo<-function(x) {
			options(warn=-1) #do this so that warnings and error get surpressed. Gets generated from exp(param)=Inf in phylogmean
			#print ("in tryFoo")
			#badLnL=100000
			result<-try(foo(x), silent=T)
			if (is.finite(result) && result < log(.Machine$double.xmax)) {
				return(result)
			}
			else {
				#print("in badLnL")
				return(badLnL)
			}
		options(warn=0) 
		}
		
		foo<-function(x) {
			t<-exponentialchangeTree(tree, a=x[2])

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}

		o<-nlm(tryFoo, p=start)
		#o<-optim(foo, p=start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$minimum, beta=exp(o$estimate[1]), a=exp(o$estimate[2]))				
		#results<-list(lnl=-o$value, beta= exp(o$par[1]), a=o$par[2])	
	}
	
	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k<-k
	return(results) 

}


dmvnormPseudoinverse<-function (x, mean, sigma, log = FALSE) {
	#print("in dmvnormPseudoinverse")
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    #print("sigma")
    #print(sigma)
    #print(" starting mahalanobis") 
    distval <- mahalanobis(x, center = mean, cov = pseudoinverse(sigma),inverted=TRUE)
    #print(" finishing mahalanobis") 

    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    #print("leaving dmvnormPseudoinverse")
    if (log) 
        return(logretval)
    exp(logretval)
}

phylogMean<-function(phyvcv, data) {
	o<-rep(1, length(data))
	#print("phyvcv in phylogMean")
	#print(phyvcv)
	ci<-pseudoinverse(phyvcv)

	#print("ci in phylogMean")
	#print(ci)
	
	#print("t(o) %*% ci %*% o")
	#print(t(o) %*% ci %*% o)
	m1<-pseudoinverse(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	return(m1 %*% m2)
	
}
	
	
ouMatrix <- function(vcvMatrix, alpha) {
## follows Hansen 1997; does not assume ultrametricity (AH 12 dec 07)
## vectorized by LJH
  vcvDiag<-diag(vcvMatrix)
  diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
  diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)

  Tij = diagi + diagj - (2 * vcvMatrix)
    
  vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
  return(vcvRescaled) 
}
    	