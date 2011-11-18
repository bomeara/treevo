plotPosteriors<-function(particleDataFrame, PriorMatrix, realParam=FALSE, realParamValues=NA) {
# particleDataFrame can be single or a list of particleDataFrames (1:n)
# priors can also be single matrix or a list of matrices (Note that priors have to be the same to make comparison across runs, therefore if a list of priors is given, this function will extract only the first matrix)
#priors should be $PriorMatrix from TreEvo output
#realParamValues should include a real value for every prior (fixed or not).

x<-particleDataFrame
priors<-PriorMatrix
if(class(x)=="data.frame"){
	data1<-subset(x[which(x[,6]>0),], generation==max(x[,1])) #make generation and other names by column so it works for partial and complete 
	run<-rep(1, dim(data1)[1])
	all<-cbind(run, data1)
}
	
if(class(x)=="list"){
	all<-data.frame()
	for (list in 1:length(x)) {
		data1<-subset(x[[list]][which(x[[list]][,6]>0),], generation==max(x[[list]]$generation)) #make generation and other names by column so it works for partial and complete 
		run<-rep(list, dim(data1)[1])
		x[[list]]<-cbind(run, data1)
		all<-rbind(all, x[[list]])		
	}
}
if (class(priors)=="list"){
	priors<-priors[[1]]	
}


freeParams<-dim(subset(priors[,which(priors[1,] != "fixed")])[,])[2] 
dev.new(width=2.5*freeParams, height=3)
nf<-layout(matrix(1:freeParams, nrow=1, byrow=TRUE), respect=TRUE)
layout.show(nf)

#alternatively, we can plotPrior above the posterior density to check if it is moving
#nf<-layout(matrix(1:(2*paramsToPlot), nrow=2, byrow=FALSE), respect=TRUE)
#layout.show(nf)

v<-vector("list", max(all$run))
nParticles<-dim(subset(all[which(all$weight>0),], run==max(all$run)))[1]

for (param in 1:dim(priors)[2]) {
	#print(param)
	#v<-vector("list", max(all$run))
	which.param<-param+7
	r<-c()
	q<-c()

	if (priors[1,param] != "fixed") {

		if (priors[1, param] == "uniform" || priors[1, param] == "exponential") {
			for (dens in 1:max(all$run)) {
				v[[dens]]<-density(subset(all[which(all$run==dens),],)[,which.param], weights=nParticles*subset(all[which(all$weight>0),], run==dens)[,7]/sum(nParticles*subset(all[which(all$weight>0),], run==dens)[,7]), from=min(subset(all[which(all$run==dens),],)[,which.param]), to=max(subset(all[which(all$run==dens),],)[,which.param]))
				v[[dens]]$x<-c(min(v[[dens]]$x), v[[dens]]$x, max(v[[dens]]$x))
				v[[dens]]$y<-c(0, v[[dens]]$y, 0)
				q<-c(q, v[[dens]]$x)
				r<-c(r, v[[dens]]$y)
			}
		}

		else if (priors[1, param] == "normal" || priors[1, param] == "lognormal" || priors[1, param] == "gamma"){
			for (dens in 1:max(all$run)) {
				v[[dens]]<-density(subset(all[which(all$run==dens),],)[,which.param], weights=nParticles*subset(all[which(all$weight>0),], run==dens)[,7]/sum(nParticles*subset(all[which(all$weight>0),], run==dens)[,7]))
				q<-c(q, v[[dens]]$x)
				r<-c(r, v[[dens]]$y)
			}
		}

		if (priors[1, param]=="uniform") {
			#print("made it to uniform priorFn")
			min=as.numeric(min(priors[2, param], priors[3, param]))
			max=as.numeric(max(priors[2, param], priors[3, param]))
			x<-runif(1000, min, max)
			poly<-curve(dunif(x), xlim=c(min(min-(.3*(max-min)), min(q)), max(max+(.3*(max-min)), max(q))), ylim=c(0, max(1.5, max(r))), type="n", xlab=names(all)[which.param], ylab="Density", bty="n")
			rect(min, 0, max, 1, border=rgb(1, 0, 0), lwd=1.5)
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[1, param]=="normal") {
			#print("made it to normal priorFn")
			mean=as.numeric(priors[2, param])
			stdev=as.numeric(priors[3, param])
			x<-rnorm(1000, mean, stdev)
			w<-density(x)
			poly<-curve(dnorm(x, mean, stdev), from=min(x), to=max(x), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n")
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(max(w$y), max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[1, param]=="lognormal") {  
			mean=as.numeric(priors[2, param])
			stdev=as.numeric(priors[3, param])
			x<-rlnorm(1000, mean, stdev)
			w<-density(x)
			poly<-curve(dlnorm(x, mean, stdev), from=0, to=qlnorm(0.99, mean, stdev), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n")
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(max(w$y), max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[1, param]=="gamma") {
			shape=as.numeric(priors[2, param])
			scale=as.numeric(priors[3, param])
			x<-rgamma(1000, shape, scale)
			w<-density(x)
			poly<-curve(dgamma(x, shape, scale), from=0, to=qgamma(0.99, shape, scale), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n")
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(max(w$y), max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[1, param]=="exponential") {
			#print("made it to exponential priorFn")
			rate=as.numeric(priors[2, param])
			x<-rexp(1000, rate)
			w<-density(x)
			poly<-curve(dexp(x, rate), from=0, to=qexp(0.99, rate), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n") #, sub=names(data)[which.param]
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(max(w$y), max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		
		
	## plotting	
		
		#plot(q, r, type="n", ylab="density", xlim=c(min(q), max(q)), ylim=c(min(r), max(r)), sub=names(all)[which.param])			
		for (lines in 1:length(v)) {
			polygon(v[[lines]]$x, v[[lines]]$y, border=NA, col=rgb(0, 0, 0, (.7/max(all$run))))
		}
	} #if (priors[1,param] != "fixed") 
} #for

	
}