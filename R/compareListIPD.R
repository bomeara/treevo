#' Plot Interparticle Distance between Runs
#' 
#' This function plots interparticle distance (IPD) between runs for each free
#' parameter.
#' 
#' 

#' @inheritParams credibleInt

#' @param verbose Commented screen output.

#' @return Produces a plot with IPD between runs per generation.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished


#' @examples
#'
#' data(simRun)
#' 
#' compareListIPD(results$particleDataFrame, verbose=FALSE)
#' 


#will compare IPD from a list of particleDataFrames.  Calculates the ratio of IPD


#' @name compareListIPD
#' @rdname compareListIPD
#' @export
compareListIPD<-function(particleDataFrame, verbose=FALSE){  #list of particleDataFrames
	params<-dim(particleDataFrame[[1]][7:dim(particleDataFrame[[1]])[2]])[2]
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
					IPDmatrix[row, col]<-(median(interparticleDistance(data1[[row]][which(data1[[row]]$generation==gen),6+param],
						data1[[col]][which(data1[[col]]$generation==gen),6+param]))^2)  / median(
						interparticleDistance(data1[[row]][which(data1[[row]]$generation==gen),6+param],
						data1[[col]][which(data1[[col]]$generation==gen),6+param]))*2
					#protect against NAs, since we are dividing above (if sd(A and or B[param] = 0 then it will NA))
					if (is.na(IPDmatrix[row, col])){  
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

