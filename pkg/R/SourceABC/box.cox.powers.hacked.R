
box.cox.powers.hacked<-function(X, start=NULL, hypotheses=NULL, ...){
	#modified from the car package to change the optim fn. The rest of the code is from J. Fox, author of car.
    modified.power<-function(x, lambda, gm){
        if (lambda == 0) log(x)*gm
        else (gm^(1-lambda))*((x^lambda)-1)/lambda
        }
    neg.kernel.profile.logL<-function(X, lambda, gm){
        for (j in 1:ncol(X)){
            X[, j]<-modified.power(X[, j], lambda[j], gm[j])
            }
        (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
        }
    univ.neg.kernel.logL <- function(x, lambda, gm){
        x <- modified.power(x, lambda, gm)
        (length(x)/2)*log(((length(x)-1)/length(x))*var(x))
        }
    X<-as.matrix(X)
    nc <- ncol(X)
    if(any(X<=0)) stop("All values must be > 0")
    gm<-apply(X, 2, function(x) exp(mean(log(x))))
    if (is.null(start)) {
        start <- rep(1, nc)
        for (j in 1:nc){
            res<- optimize(
                f = function(lambda) univ.neg.kernel.logL(x=X[, j], lambda=lambda, gm=gm[j]), 
                lower=-50, upper=+50)
            start[j] <- res$minimum
            }
        }
    res<-optim(start, neg.kernel.profile.logL, hessian=TRUE, method="BFGS", X=X, gm=gm, ...)
    result<-list()
    result$start<-start
    result$criterion<-res$value
    result$names<-colnames(X)
    result$lambda<-res$par
    result$stderr<-sqrt(diag(inv(res$hessian)))
    result$LR0<-2*(neg.kernel.profile.logL(X, rep(0, nc), gm)-res$value)
    result$LR1<-2*(neg.kernel.profile.logL(X, rep(1, nc), gm)-res$value)
    if (!is.null(hypotheses)) {
        for (i in 1:length(hypotheses)){
            if (length(hypotheses[[i]]) != nc) 
                stop(paste("hypothesis", i, "that powers =", hypotheses[[i]], "does not have", nc, "values"))
            hypotheses[[i]] <- list(test=2*(neg.kernel.profile.logL(X, hypotheses[[i]], gm)-res$value), 
                hypothesis=hypotheses[[i]])
            }
        result$hypotheses <- hypotheses
        }
    result$return.code<-res$convergence
    if(result$return.code != 0) 
        warning(paste("Convergence failure: return code =", 
            result$return.code))
    class(result)<-"box.cox.powers"
    result
    }
