library(geiger)
library(stats)
summaryStatsLong<-function(phy, data, todo=c(), jobName="") {
	sink(file="/dev/null")
	if (length(todo)==0) {
		todo=rep(1, 22+dim(data)[1]) #by default, include everything -- the 22 summary stats and the raw tip data
	}
	
	#thought here: include brown<-try(fitContinouous()), store brown, then do try(brown$lnl), etc. Faster than calling each fn
	brown<-try(fitContinuous.hacked(phy=phy, data=data[max(todo[1:3])], model="BM")[[1]]) #will only run if want to do brownian summary stats
	#brown.lnl<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[1]], model="BM")[[1]]$lnl)) #if todo[i]==0, will cause an error right away, saving on computation time
	brown.lnl<-as.numeric(try(brown$lnl/todo[1])) #divide by zero so we get Inf if we don't want that summary stat
	#brown.beta<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[2]], model="BM")[[1]]$beta))
	brown.beta <-as.numeric(try(brown$beta/todo[2]))
	#brown.aic<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[3]], model="BM")[[1]]$aic))
	brown.aic <-as.numeric(try(brown$aic/todo[3]))

	lambda<-try(fitContinuous.hacked(phy=phy, data=data[max(todo[4:7])], model="lambda")[[1]])
	#lambda.lnl<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[4]], model="lambda")[[1]]$lnl))
	lambda.lnl <-as.numeric(try(lambda$lnl/todo[4]))
	#lambda.beta<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[5]], model="lambda")[[1]]$beta))
	lambda.beta <-as.numeric(try(lambda$beta/todo[5]))
	#lambda.lambda<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[6]], model="lambda")[[1]]$lambda))
	lambda.lambda <-as.numeric(try(lambda$lambda/todo[6]))
	#lambda.aic<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[7]], model="lambda")[[1]]$aic))
	lambda.aic <-as.numeric(try(lambda$aic/todo[7]))

	delta<-try(fitContinuous.hacked(phy=phy, data=data[max(todo[8:11])], model="delta")[[1]])
	#delta.lnl<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[8]], model="delta")[[1]]$lnl))
	delta.lnl <-as.numeric(try(delta$lnl/todo[8]))
	#delta.beta<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[9]], model="delta")[[1]]$beta))
	delta.beta <-as.numeric(try(delta$beta/todo[9]))
	#delta.delta<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[10]], model="delta")[[1]]$delta))
	delta.delta <-as.numeric(try(delta$delta/todo[10]))
	#delta.aic<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[11]], model="delta")[[1]]$aic))
	delta.aic <-as.numeric(try(delta$aic/todo[11]))
	
	ou<-try(fitContinuous.hacked(phy=phy, data=data[max(todo[12:15])], model="OU")[[1]])
	#ou.lnl<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[12]], model="OU")[[1]]$lnl))
	ou.lnl <-as.numeric(try(ou$lnl/todo[12]))
	#ou.beta<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[13]], model="OU")[[1]]$beta))
	ou.beta <-as.numeric(try(ou$beta/todo[13]))
	#ou.alpha<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[14]], model="OU")[[1]]$alpha))
	ou.alpha <-as.numeric(try(ou$alpha/todo[14]))
	#ou.aic<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[15]], model="OU")[[1]]$aic))
	ou.aic <-as.numeric(try(ou$aic/todo[15]))
	
	white<-try(fitContinuous.hacked(phy=phy, data=data[max(todo[16:17])], model="white")[[1]])
	#white.lnl<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[16]], model="white")[[1]]$lnl))
	white.lnl<-as.numeric(try(white$lnl/todo[16]))
	#white.aic<-as.numeric(try(fitContinuous.hacked(phy=phy, data=data[todo[17]], model="white")[[1]]$aic))
	white.aic<-as.numeric(try(white$aic/todo[17]))
	
	
	raw.mean<-as.numeric(try(mean(data)/todo[18]))
	raw.max<-as.numeric(try(max(data)/todo[19]))
	raw.min<-as.numeric(try(min(data)/todo[20]))
	raw.var<-as.numeric(try(var(data)/todo[21]))
	raw.median<-as.numeric(try(median(data[,])/todo[22]))	#cat("summaryStatsLong")
	summarystats<-c(brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic, raw.mean, raw.max, raw.min, raw.var, raw.median, data[[1]] )
	#cat("\n summaryStatsLong summarystats1\n")
	#print(summarystats)
	summarystats[which(todo==0)]<-NA
	#cat("\n summaryStatsLong summarystats2\n")
	#print(summarystats)
	summarystats[which(is.finite(summarystats)==FALSE)]<-NA
	
	somethingfailed=FALSE
	failurevector<-rep(NA,length(summarystats))
	for (i in 1:length(summarystats)) {
		if (todo[i]==1) {
			if (is.na(summarystats[i])) {
				somethingfailed=TRUE
				failurevector[i]=2 #failure!
				#todo=0
			}
			else {
				failurevector[i]=1 #success!
			}
		}
		else {
			failurevector[i]=0 #means was not tried
		}
	}
	print (paste("failurevector = ",failurevector,collapse="",sep=""))
	if (somethingfailed) {
		failed.summarystats<-c(dput(brown), dput(alpha), dput(white), dput(white))
		GeigerFailure<-vector("list", 2)
		GeigerFailure[[1]]<-failurevector
		GeigerFailure[[2]]<-failed.summarystats
		save(GeigerFailure, file=paste("geigerFailure", jobName, sep=""))
	}
	
	#cat("\n summaryStatsLong summarystats3\n")
	#print(sstummarystats)
	while(sink.number()>0) {sink()}
	summarystats
}
