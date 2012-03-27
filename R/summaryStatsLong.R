library(geiger)
library(stats)
summaryStatsLong<-function(phy, data, todo=c(), jobName="") {
	if(any(phy$edge.length==0)){
		return("There are zero branch lengths in your tree--will not run properly")
	}
	sink(file="/dev/null")
	if (length(todo)==0) {
		todo=rep(1, 22+dim(data)[1]) #by default, include everything -- the 22 summary stats and the raw tip data
	}
	
	brown<-fitContinuous.hacked(phy=phy, data=data, model="BM")[[1]] #will only run if want to do brownian summary stats
	brown.lnl<-as.numeric(brown$lnl/todo[1]) #divide by zero so we get Inf if we don't want that summary stat
	brown.beta <-as.numeric(brown$beta/todo[2])
	brown.aic <-as.numeric(brown$aic/todo[3])

	lambda<-fitContinuous.hacked(phy=phy, data=data, model="lambda")[[1]]
	lambda.lnl <-as.numeric(lambda$lnl/todo[4])
	lambda.beta <-as.numeric(lambda$beta/todo[5])
	lambda.lambda <-as.numeric(lambda$lambda/todo[6])
	lambda.aic <-as.numeric(lambda$aic/todo[7])

	delta<-fitContinuous.hacked(phy=phy, data=data, model="delta")[[1]]
	delta.lnl <-as.numeric(delta$lnl/todo[8])
	delta.beta <-as.numeric(delta$beta/todo[9])
	delta.delta <-as.numeric(delta$delta/todo[10])
	delta.aic <-as.numeric(delta$aic/todo[11])
	
	ou<-fitContinuous.hacked(phy=phy, data=data, model="OU")[[1]]
	ou.lnl <-as.numeric(ou$lnl/todo[12])
	ou.beta <-as.numeric(ou$beta/todo[13])
	ou.alpha <-as.numeric(ou$alpha/todo[14])
	ou.aic <-as.numeric(ou$aic/todo[15])
	
	white<-fitContinuous.hacked(phy=phy, data=data, model="white")[[1]]
	white.lnl<-as.numeric(white$lnl/todo[16])
	white.aic<-as.numeric(white$aic/todo[17])
	
	
	raw.mean<-as.numeric(mean(data)/todo[18])
	raw.max<-as.numeric(max(data)/todo[19])
	raw.min<-as.numeric(min(data)/todo[20])
	raw.var<-as.numeric(var(data)/todo[21])
	raw.median<-as.numeric(median(data[,])/todo[22])	#cat("summaryStatsLong")
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
		failed.summarystats<-c(dput(brown), dput(lambda), dput(delta), dput(ou), dput(white))
		GeigerFailure<-vector("list", 4)
		GeigerFailure[[1]]<-failurevector
		GeigerFailure[[2]]<-failed.summarystats
		GeigerFailure[[3]]<-phy
		GeigerFailure[[4]]<-data
		save(GeigerFailure, file=paste("geigerFailure", jobName, ".Rdata", sep=""))
	}
	
	#cat("\n summaryStatsLong summarystats3\n")
	#print(sstummarystats)
	while(sink.number()>0) {sink()}
	summarystats
}
