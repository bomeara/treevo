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
#Calculate IPD either within each generation or between parent-offspring gens
compareIPD<-function(particleDataFrame1, particleDataFrame2){
	if (any(names(particleDataFrame1) != names(particleDataFrame2))){
		warning("particleDataFrame names do not match, likely trying to compare different models")
	}
	runIPD<-matrix(nrow=max(particleDataFrame1$generation), ncol=dim(particleDataFrame1)[2]-6)
	for (gen in 1:dim(runIPD)[1]){ #each parent-offspring generation
		subpDF_1<-as.data.frame(subset(particleDataFrame1[which(particleDataFrame1$weight>0),], generation==gen)[7:dim(particleDataFrame1)[2]])
		subpDF_2<-as.data.frame(subset(particleDataFrame2[which(particleDataFrame2$weight>0),], generation==gen)[7:dim(particleDataFrame2)[2]])
		rownames(runIPD)<-paste("gen", c(1:max(particleDataFrame1$generation)), sep="")
		colnames(runIPD)<-names(subpDF)
		for (param in 1:dim(runIPD)[2]){ #each param
			IPD<-median(interparticleDistance(subpDF_1[,param], subpDF_2[,param]))
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
