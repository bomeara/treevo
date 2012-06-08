#calculates the distance between each possible pairs of points, where one point must come from each vector
#Brians function
interparticleDistance<-function(x,y,abs=TRUE) {
   distances<-matrix(nrow=length(x),ncol=length(y))
   for (i in sequence(length(x))) {
      for (j in sequence(length(y))) {
        distances[i,j]<-x[i]-y[j] 
      }
   }
   if (abs==TRUE) {
      distances<-abs(distances) 
   }
   return(distances)
}



#Calculate IPD either within each generation or between parent-offspring gens
IPD<-function(particleDataFrame, plot=match.arg(arg=plot,choices=c("within", "between"),several.ok=FALSE)){
	if(plot=="within"){
		runIPD<-matrix(nrow=max(particleDataFrame$generation), ncol=dim(particleDataFrame)[2]-6)
		for (gen in 1:dim(runIPD)[1]){ #each generation
			subpDF<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==gen)[7:dim(particleDataFrame)[2]])
			rownames(runIPD)<-paste("gen", c(1:max((particleDataFrame$generation))), sep="")
			colnames(runIPD)<-names(subpDF)
			for (param in 1:dim(runIPD)[2]){ #each param
				IPD<-median(interparticleDistance(subpDF[,param], subpDF[,param]))
				runIPD[gen, param]<-IPD
			}
		}
	}
	else if(plot=="between"){
		runIPD<-matrix(nrow=max(particleDataFrame$generation)-1, ncol=dim(particleDataFrame)[2]-6)
		for (gen in 1:dim(runIPD)[1]){ #each parent-offspring generation
			subpDF_parent<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==gen)[7:dim(particleDataFrame)[2]])
			subpDF_offspring<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==gen+1)[7:dim(particleDataFrame)[2]])
			rownames(runIPD)<-paste("gen", c(1:max(particleDataFrame$generation-1)), "-", c(2:max(particleDataFrame$generation)), sep="")
			colnames(runIPD)<-names(subpDF_parent)
			for (param in 1:dim(runIPD)[2]){ #each param
				IPD<-median(interparticleDistance(subpDF_parent[,param], subpDF_offspring[,param]))
				print(IPD)
				runIPD[gen, param]<-IPD
	#			#paramIPD<-append(runIPD, IPD)
			}
		}
	}
	
	
	nf<-layout(matrix(1:ncol(runIPD), nrow=1, byrow=TRUE), respect=TRUE)
	layout.show(nf)
	for(i in 1:param){
		plot(c(1:gen), runIPD[,i])
	}
	return(runIPD)
}



#Calculate IPD as a ratio of between/within(parent gen)
IPD_ratio<-function(particleDataFrame, plotDensity=T){  #include list of L?
	runIPD<-matrix(nrow=max(particleDataFrame$generation)-1, ncol=dim(particleDataFrame)[2]-6)
	for (gen in 1:dim(runIPD)[1]){ #each parent-offspring generation
		subpDF_parent<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==gen)[7:dim(particleDataFrame)[2]])
		subpDF_offspring<-as.data.frame(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==gen+1)[7:dim(particleDataFrame)[2]])
		rownames(runIPD)<-paste("gen", c(1:max(particleDataFrame$generation-1)), "-", c(2:max(particleDataFrame$generation)), "/", c(1:max(particleDataFrame$generation-1)), sep="")
		colnames(runIPD)<-names(subpDF_parent)
		for (param in 1:dim(runIPD)[2]){ #each param
			IPD<-median(interparticleDistance(subpDF_parent[,param], subpDF_offspring[,param]))/median(interparticleDistance(subpDF_parent[,param], subpDF_parent[,param]))
			if (is.na(IPD)){ #divide by 0, get NaN
				IPD<-0
			}
			print(IPD)
			runIPD[gen, param]<-IPD
	#		#paramIPD<-append(runIPD, IPD)
		}
	}
	if (plotDensity){
		nf<-layout(matrix(1:(2*ncol(runIPD)), nrow=2, byrow=F), respect=TRUE)
		layout.show(nf)
		for(i in 1:param){
			plot(c(1:gen), runIPD[,i])
			plot(density(subpDF_parent[,i]))
		}	
	}
	
	else{
		nf<-layout(matrix(1:ncol(runIPD), nrow=1, byrow=TRUE), respect=TRUE)
		layout.show(nf)
		for(i in 1:param){
			plot(c(1:gen), runIPD[,i])
		}
	}
	return(runIPD)
}





#move towards calculating IPD between gens of different runs to see if they are converging or not. 
#Calculate IPD as a ratio of medIPDs from two runs
compareIPD<-function(particleDataFrame1, particleDataFrame2){
	if (any(names(particleDataFrame1) != names(particleDataFrame2))){
		warning("particleDataFrame names do not match, likely trying to compare different models")
	}
	particleDataFrame1<-convertabcResultsToTreEvoResults(particleDataFrame1)
	particleDataFrame2<-convertabcResultsToTreEvoResults(particleDataFrame2)
	runIPD<-matrix(nrow=max(particleDataFrame1$generation), ncol=dim(particleDataFrame1)[2]-6)
	for (gen in 1:dim(runIPD)[1]){ #each parent-offspring generation
		subpDF_1<-as.data.frame(subset(particleDataFrame1[which(particleDataFrame1$weight>0),], generation==gen)[7:dim(particleDataFrame1)[2]])
		subpDF_2<-as.data.frame(subset(particleDataFrame2[which(particleDataFrame2$weight>0),], generation==gen)[7:dim(particleDataFrame2)[2]])
		rownames(runIPD)<-paste("gen", c(1:max(particleDataFrame1$generation)), sep="")
		colnames(runIPD)<-names(subpDF_1)
		for (param in 1:dim(runIPD)[2]){ #each param
			IPD<-(median(interparticleDistance(subpDF_1[,param], subpDF_2[,param]))^2)  / median(interparticleDistance(subpDF_1[,param], subpDF_2[,param]))*2
			print(IPD)
			runIPD[gen, param]<-IPD
#			#paramIPD<-append(runIPD, IPD)
		}
	}
	
	nf<-layout(matrix(1:ncol(runIPD), nrow=1, byrow=TRUE), respect=TRUE)
	layout.show(nf)
	for(i in 1:param){
		plot(c(1:gen), runIPD[,i])
	}
	return(runIPD)
}

#will compare IPD from a list of particleDataFrames.  Calculates the ratio of IPD
compareListIPD<-function(particleDataFrame, verbose=F){  #list of particleDataFrames
	params<-dim(all.a[[1]][7:dim(all.a[[1]])[2]])[2]
plot.new()
nf<-layout(matrix(1:params, nrow=1, byrow=TRUE), respect=TRUE)
	layout.show(nf)	
	data1<-vector("list")
	maxgen<-c()
	for (list in 1:length(particleDataFrame)) {
		data1[[list]]<-subset(particleDataFrame[[list]][which(particleDataFrame[[list]][,6]>0),], ) 
		maxgen<-append(maxgen, max(data1[[list]]$generation))
	}	
	for (param in 1:params){
		genvector<-vector()
		IPDvector<-vector()
		for (gen in 1:max(maxgen)){
			IPDmatrix<-matrix(nrow=length(data1), ncol=length(data1))
			for (row in 1: dim(IPDmatrix)[1]){
				for (col in 1: dim(IPDmatrix)[2]){
					IPDmatrix[row, col]<-(median(interparticleDistance(data1[[row]][which(data1[[row]]$generation==gen),6+param], data1[[col]][which(data1[[col]]$generation==gen),6+param]))^2)  / median(interparticleDistance(data1[[row]][which(data1[[row]]$generation==gen),6+param], data1[[col]][which(data1[[col]]$generation==gen),6+param]))*2
					if (is.na(IPDmatrix[row, col])){  #protect against NAs, since we are dividing above (if sd(A and or B[param] = 0 then it will NA))
						IPDmatrix[row, col]<-0
					}
				}
			}
			#print(paste("param = ", param))
			#print(paste("generation = ", gen))
			genvector<-append(genvector, rep(gen, length(unique(as.vector(IPDmatrix)))))
			IPDvector<-append(IPDvector, unique(as.vector(IPDmatrix)))
			#print(genvector)
			#print(IPDvector)
			if(verbose){
				print(paste("param=", names(data1[[1]][6+param])))
				print(cbind(genvector, IPDvector))
			}
		}
		#pdf(paste("compareListIPD", jobName, ".pdf", sep=""))
		plot(genvector, IPDvector, xlab="generation", ylab="IPD", sub=names(data1[[1]][6+param]))
	}	
}

