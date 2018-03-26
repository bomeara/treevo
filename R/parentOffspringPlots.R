#' Plotting Parent-Offspring Particles
#' 
#' This function uses the \code{particleDataFrame} output by TreEvo ABC analyses and plots
#' parent-offspring particles from generation to generation.
#' 
#' Each parameter is plotted twice for parent-offspring relationships through
#' the generations.  In the top row, particles are plotted as a measure of
#' distance to the observed data; the farther away the particle the bigger the
#' circle.  In the bottom row,
#' particles are plotted as a measure of their weights; larger circles are
#' closer to observed data and therefore carry more weight in the analysis.
#' 
#' As of version 0.6.0, rejected particles are not saved for outputting by the parallelized algorithm,
#' and thus they are no longer displayed by this function.

# Grayed out stars represent rejected particles.
#Gray circles indicate rejected particles.  

#' @param particleDataFrame \code{particleDataFrame} output by TreEvo ABC analyses.

#' @return Creates a set of parent-offspring plots.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords parentOffspringPlots

#' @examples
#' 
#' data(simRunExample)
#' parentOffspringPlots(results$particleDataFrame)
#' 

#' @name parentOffspringPlots
#' @rdname parentOffspringPlots
#' @export
parentOffspringPlots<-function(particleDataFrame){

	
	nparams<-dim(particleDataFrame)[2]-6
	nb<-2*nparams
	nf<-layout(matrix(1:nb, nrow=2, byrow=TRUE), respect=TRUE)
	#layout.show(nf)

	for (param in 1:nparams) {
		param.position<-param+6
				
		plot(particleDataFrame[,param.position], particleDataFrame$generation, 
			xlab=colnames(particleDataFrame)[param.position], ylab="Generation", type="n")
		title("size as measure of distance")
		#kept<-subset(particleDataFrame[which(particleDataFrame$id>0),])[,]	
		#reject<-subset(particleDataFrame[which(particleDataFrame$id<0),])[,]
		#
		kept<-particleDataFrame[particleDataFrame$id>0,]
		short.kept<-kept[kept$generation>1,]

		#for (i in 1:(dim(reject)[1])) {
		#	circle.size<-(reject[i, 5]/max(reject[,5]))*(0.05*(max(particleDataFrame[,param.position])-min(particleDataFrame[,param.position])))
		#	symbols(x=reject[i, param.position], y=reject[i, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="gray")
		#}	
		
		for (j in 1:(dim(kept)[1])) {
			#circle.size<-(kept[j, 5]/max(kept[,5]))*mean(particleDataFrame[,param.position])
			circle.size<-(kept[j, 5]/max(kept[,5]))*(0.05*(max(kept[,param.position])-min(kept[,param.position])))
			symbols(x=kept[j, param.position], y=kept[j, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="black")
		}		
		
		for (k in 1:(dim(short.kept)[1])) {
			prev.gen<-kept[kept$generation==short.kept[k, 1]-1,]  #works to retreive prev gen
			lines(c(short.kept[k, param.position], 
				prev.gen[short.kept[k,]$parentid, param.position]), 
				c(short.kept[k, 1], short.kept[k, 1]-1))
		}	
	}		
		
		
	##-----For particle weights

	for (param in 1:nparams) {
		param.position<-param+6

		plot(particleDataFrame[,param.position], particleDataFrame$generation, 
			xlab=colnames(particleDataFrame)[param.position], ylab="Generation", type="n")
		title("size as measure of particle weights")
		
		#kept<-subset(particleDataFrame[which(particleDataFrame$id>0),])[,]	
		#reject<-subset(particleDataFrame[which(particleDataFrame$id<0),])[,]
		#short.kept<-subset(kept[which(kept$generation>1),])[,]

		#for (i in 1:(dim(reject)[1])) {
		#	points(x=reject[i, param.position], y=reject[i, 1], col="gray", pch=8)
		#}	
		
		for (j in 1:(dim(kept)[1])) {
			circle.size<-(kept[j, 6]/max(kept[,6]))*(0.05*(max(kept[,param.position])-min(kept[,param.position])))
			symbols(x=kept[j, param.position], y=kept[j, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="black")
		}	
		
		for (k in 1:(dim(short.kept)[1])) {
			prev.gen<-kept[kept$generation==short.kept[k, 1]-1,]  #works to retreive prev gen
			lines(c(short.kept[k, param.position], 
				prev.gen[short.kept[k,]$parentid, param.position]), 
				
				c(short.kept[k, 1], short.kept[k, 1]-1))
		}	
	}	
	layout(1)
}		
