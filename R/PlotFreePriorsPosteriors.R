plotParamPosteriors<-function(data, priors, realParam=FALSE, realParamValues=NA) {
# data can be a single particleDataFrame or a list of particleDataFrames (1:n)
# priors can also be single matrix or a list of matrices (Note that priors have to be the same to make comparison across runs, therefore if a list of priors is given, this function will extract only the first matrix)


if(class(data)=="data.frame"){
	data1<-subset(data[which(data[,6]>0),], generation==max(data[,1])) #make generation and other names by column so it works for partial and complete 
	run<-rep(1, dim(data1)[1])
	all<-cbind(run, data1)
}
	
if(class(data)=="list"){
	all<-data.frame()
	for (list in 1:length(data)) {
		data1<-subset(data[[list]][which(data[[list]][,6]>0),], generation==max(data[[list]][,1])) #make generation and other names by column so it works for partial and complete 
		run<-rep(list, dim(data1)[1])
		data[[list]]<-cbind(run, data1)
		all<-rbind(all, data[[list]])		
	}
}
if (class(priors)=="list"){
	priors<-priors[[1]]	
}


dev.new(width=2.5*(dim(priors)[1]), height=3)
nf<-layout(matrix(1:dim(priors)[1], nrow=1, byrow=TRUE), respect=TRUE)
layout.show(nf)

#alternatively, we can plotPrior above the posterior density to check if it is moving
#nf<-layout(matrix(1:(2*paramsToPlot), nrow=2, byrow=FALSE), respect=TRUE)
#layout.show(nf)

v<-vector("list", max(all$run))
nParticles<-dim(subset(all[which(all$weight>0),], run==max(all$run)))[1]

for (param in 1:dim(priors)[1]) {
	#print(param)
	v<-vector("list", max(all$run))
	which.param<-param+7
	print(which.param)
	r<-c()
	q<-c()

		if (priors[param, 1] == "uniform" || priors[param, 1] == "exponential") {
			
			for (dens in 1:max(all$run)) {
				v[[dens]]<-density(subset(all[which(all$run==dens),],)[,which.param], weights=nParticles*subset(all[which(all$weight>0),], run==dens)[,7]/sum(nParticles*subset(all[which(all$weight>0),], run==dens)[,7]), from=min(subset(all[which(all$run==dens),],)[,which.param]), to=max(subset(all[which(all$run==dens),],)[,which.param]))
				v[[dens]]$x<-c(min(v[[dens]]$x), v[[dens]]$x, max(v[[dens]]$x))
				v[[dens]]$y<-c(0, v[[dens]]$y, 0)
				q<-c(q, v[[dens]]$x)
				r<-c(r, v[[dens]]$y)
			}
		}
		else {
			for (dens in 1:max(all$run)) {
				v[[dens]]<-density(subset(all[which(all$run==dens),],)[,which.param], weights=nParticles*subset(all[which(all$weight>0),], run==dens)[,7]/sum(nParticles*subset(all[which(all$weight>0),], run==dens)[,7]))
				q<-c(q, v[[dens]]$x)
				r<-c(r, v[[dens]]$y)
			}
		}


		if (priors[param,1]=="uniform") {
			#print("made it to uniform priorFn")
			min=as.numeric(priors[param, 2])
			max=as.numeric(priors[param,3])
			x<-runif(1000, min, max)
			poly<-curve(dunif(x), xlim=c(min(min-(.3*(max-min)), min(q)), max(max+(.3*(max-min)), max(q))), ylim=c(0, max(1.5, max(r))), type="n", xlab=names(all)[which.param], ylab="Density", bty="n", sub=names(data)[which.param])
			rect(min, 0, max, 1, border=rgb(1, 0, 0), lwd=1.5)
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[param,1]=="normal") {
			#print("made it to normal priorFn")
			mean=as.numeric(priors[param, 2])
			stdev=as.numeric(priors[param,3])
			x<-rnorm(1000, mean, stdev)
			w<-density(x)
			poly<-curve(dnorm(x, mean, stdev), from=min(x), to=max(x), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n", sub=names(data)[which.param])
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[param,1]=="lognormal") {  
			mean=as.numeric(priors[param, 2])
			stdev=as.numeric(priors[param, 3])
			x<-rlnorm(1000, mean, stdev)
			w<-density(x)
			poly<-curve(dlnorm(x, mean, stdev), from=0, to=qlnorm(0.99, mean, stdev), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n", sub=names(data)[which.param])
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[param,1]=="gamma") {
			shape=as.numeric(priors[param, 2])
			scale=as.numeric(priors[param, 3])
			x<-rgamma(1000, shape, scale)
			w<-density(x)
			poly<-curve(dgamma(x, shape, scale), from=0, to=qgamma(0.99, shape, scale), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n", sub=names(data)[which.param])
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		else if (priors[param,1]=="exponential") {
			#print("made it to exponential priorFn")
			rate=as.numeric(priors[param, 2])
			x<-rexp(1000, rate)
			w<-density(x)
			poly<-curve(dexp(x, rate), from=0, to=qexp(0.99, rate), xlim=c(min(min(w$x), min(q)), max(max(w$x), max(q))), ylim=c(0, max(max(w$y), max(r))), xlab=names(all)[which.param], ylab="Density", col=rgb(1, 0, 0), lwd=1.5, bty="n", sub=names(data)[which.param])
			if (realParam) {
				segments(realParamValues[param], 0, realParamValues[param], max(1.5, max(r)), col=rgb(0, 0, 1), lwd=1.5)
			}
		}
		
		
	## plotting	
		
		#plot(q, r, type="n", ylab="density", xlim=c(min(q), max(q)), ylim=c(min(r), max(r)), sub=names(all)[which.param])			
		for (lines in 1:length(v)) {
			polygon(v[[lines]]$x, v[[lines]]$y, border=NA, col=rgb(0, 0, 0, (.7/max(all$run))))
		}
} #for

	
}