#' 3D ABCplots
#' 
#' Plot posterior density distribution for each generation in 3d plot window
#' 
#' This opens a new interactive 3d plotting window and plots the posterior
#' density distribution of accepted particles from each generation.  Several
#' options are available to add to the plot: plotting particles by weight or
#' distance, plotting particle parantage, and plotting the real parameter
#' values (if known).
#' 
#' @param particleDataFrame particleDataFrame output from doRun
#' @param parameter column number of parameter of interest from
#' particleDataFrame
#' @param show.particles option to show particles on 3d plot as "none" or as a
#' function of "weights" or "distance"
#' @param plot.parent option to plot lines on the floor of the 3d plot to show
#' particle parantage
#' @param realParam option to display real parameter value as a solid line,
#' also must give actual value for this (realParamValues).  Note: this should
#' only be done with simulated data where real param values are recorded
#' @param realParamValues Value for (realParam)
#' @author Barb Banbury
#' @references O'Meara and Banbury, unpublished
#' @keywords ThreeD.ABCplots
#' @examples
#' 
#' #data(res)
#' #ThreeD.ABCplots(particleDataFrame=res$particleDataFrame, parameter=7, show.particles="none", plot.parent=FALSE, realParam=FALSE, realParamValues=NA) 
#' 
ThreeD.ABCplots<-function(particleDataFrame, parameter, show.particles="none", plot.parent=FALSE, realParam=FALSE, realParamValues=NA) {
generation<-NULL #to appease R CMD CHECK
#library(gpclib)
#library(rgl)
#library(geometry)
	
x<-particleDataFrame	
param.position<-parameter
	nParticles<-dim(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame$generation)))[1]
	nparams<-dim(x)[2]-6
	
	q<-vector() #vector of x vals
	r<-vector() #vector of y vals
	s<-vector() #generation each x-y coord is found

	if (max(subset(particleDataFrame[which(particleDataFrame$weight>0),])[,param.position])-min(subset(particleDataFrame[which(particleDataFrame$weight>0),])[,param.position])!=0) {
	v<-vector("list", max(particleDataFrame$generation))
		for (i in 1:max(particleDataFrame$generation)){
			which.gen<-(i+1)-1
			v[[i]]<-density(subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==i)[,param.position], weights=nParticles*subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==i)[,6]/sum(nParticles*subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==i)[,6]))
			q<-c(q, v[[i]]$x)
			r<-c(r, v[[i]]$y)
			#s<-c(s, v[[i]]$x)	# return a$generation which v[[i]]
				for (i in 1:length(v[[i]]$x)){
					s<-c(s,which.gen)
				}
			T<-as.matrix(cbind(q, r, s))
		}
				
		x<-T[,1]
		y<-T[,2]
		z<-T[,3]
		open3d()  #make bigger window
		#bg3d("color)  #gives background color for plot

		plot3d(x, y, z, col="black", box=FALSE, type="n", xlab="", ylab="", zlab="", zlim=c(0, max(particleDataFrame$generation)), ylim=c(0, max(y)))
print(paste("HERE"))
		rgl.viewpoint(35, 1, 90)  #sets viewpoint for initial plot
		title3d(colnames(x)[param.position], col='red', pos=c(NA, -2, max(z))) 
		#text3d(x=min(x), y=mean(y), z=max(z), text="Density" col='blue') 
		#title3d("Starting States", col='purple', pos=c(NA, 0, max(which.gen))) 
		for (i in 1:max(s)){
			ngen<-(i+1)-1
			triangles<-triangulate(as(cbind(x[which(z==i)],y[which(z==i)]), "gpc.poly"))
			zfit<-predict(lm(z[which(z==i)] ~ x[which(z==i)] + y[which(z==i)]), newdata=data.frame(x=triangles[,1], y=triangles[,2]))
			opacity<-0.8*(ngen/length(v))
			rgl.material(color="black", alpha=opacity, lit=FALSE)
			#return(dim(triangles))
			#return(dim(zfit))
			triangles3d(cbind(triangles, zfit), col="red")
			#readline(prompt="hit enter ")
			if (realParam) {
				rgl.material(color="blue", lwd=2)
				lines3d(x=c(realParamValues[1], realParamValues[1]), y=c(0, 0), z=c(min(s), max(s)))	
			}
		}
	}		
		
	else {
		warning(paste("You are attempting to plot",colnames(x)[param.position],", which is a fixed value")) #return fixed val
		}
		
		
show.particles<-match.arg(arg=show.particles, choices=c("none", "weights", "distance"),several.ok=FALSE)
kept<-subset(particleDataFrame[which(particleDataFrame$id>0),])[,]	
reject<-subset(particleDataFrame[which(particleDataFrame$id<0),])[,]
short.kept<-subset(kept[which(kept$generation>1),])[,]

	if (show.particles=="none"){
		cat("currently not plotting particles.  To plot particles modify the show.particles= argument.")
		}
	else if (show.particles=="weights") {
		for (j in 1:nrow(kept)) {
			circle.size<-(kept[j, 6]/max(kept[,6]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))*10
			spheres3d(x=kept[j, param.position], y=-1*(.05*(max(y)-min(y))), z=kept[j, 1], radius=circle.size, type="s", col="black", add=TRUE, aspect=TRUE)
		}	
	}		
	else if (show.particles=="distance") {
		for (j in 1:nrow(kept)) {
			circle.size<-(kept[j, 5]/max(kept[,5]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))*15
			spheres3d(x=kept[j, param.position], y=-1*(.05*(max(y)-min(y))), z=kept[j, 1], radius=circle.size, type="s", col="black", add=TRUE, aspect=TRUE)
		}	
	}	
	if (plot.parent){
		for (k in 1:nrow(short.kept)) {
		prev.gen<-subset(kept[which(kept$generation==short.kept[k, 1]-1),])[,]  #works to retreive prev gen
		rgl.material(color="black", alpha=0.5)
		lines3d(x=c(short.kept[k, param.position], prev.gen[short.kept[k,]$parentid, param.position]), y=c(0, 0), z=c(short.kept[k, 1], short.kept[k, 1]-1))		}	
	}	
}
