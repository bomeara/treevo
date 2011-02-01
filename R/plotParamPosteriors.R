plotParamPosteriors<-function(data) {
# data can be a single particleDataFrame or a list of particleDataFrames (1:n)

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

paramsToPlot<-dim(all)[2]-7
dev.new(width=15, length=.1)
nf<-layout(matrix(1:paramsToPlot, nrow=1, byrow=FALSE), respect=TRUE)
layout.show(nf)

#alternatively, we can plotPrior above the posterior density to check if it is moving
#nf<-layout(matrix(1:(2*paramsToPlot), nrow=2, byrow=FALSE), respect=TRUE)
#layout.show(nf)

v<-vector("list", max(all$run))
nParticles<-dim(subset(all[which(all$weight>0),], run==max(all$run)))[1]

	for (param in 1: paramsToPlot) {
		v<-vector("list", max(all$run))
		which.param<-param+7
		r<-c()
		q<-c()
		
		if (range(all[,which.param])!=0) {
			for (dens in 1:max(all$run)){
				v[[dens]]<-density(subset(all[which(all$run==dens),],)[,which.param], weights=nParticles*subset(all[which(all$weight>0),], run==dens)[,7]/sum(nParticles*subset(all[which(all$weight>0),], run==dens)[,7]))
			
				q<-c(q, v[[dens]]$x)
				r<-c(r, v[[dens]]$y)
			}
		
			plot(q, r, type="n", ylab="density", xlim=c(min(q), max(q)), ylim=c(min(r), max(r)), sub=names(all)[which.param])			
			for (lines in 1:length(v)) {
			polygon(v[[lines]]$x, v[[lines]]$y, border=NA, col=rgb(0, 0, 0, (.7/max(all$run))))
			}
		}
	
		else {
			plot(q, r, type="n", ylab="density", xlim=c(-1, 1), ylim=c(0, 1), sub=names(all)[which.param])
		}
	}
}