ThreeD.ABCplots<-function(x, parameter, show.particles=none, plot.parent=FALSE) {
	
library(gpclib)
library(rgl)
library(geometry)
	
param.position<-parameter
	nParticles<-dim(subset(a[which(a$weight>0),], generation==max(a$generation)))[1]
	nparams<-dim(a)[2]-6
	#nf<-layout(matrix(1:nparams, nrow=1), respect=TRUE)  #can't do layout with 3d?
	#layout.show(nf)
	
	q<-vector() #vector of x vals
	r<-vector() #vector of y vals
	s<-vector() #generation each x-y coord is found

	if (range(subset(a[which(a$weight>0),])[,param.position])!=0) {
	v<-vector("list", max(a$generation))
		for (i in 1:max(a$generation)){
			which.gen<-(i+1)-1
			v[[i]]<-density(subset(a[which(a$weight>0),], generation==i)[,param.position], weights=nParticles*subset(a[which(a$weight>0),], generation==i)[,6]/sum(nParticles*subset(a[which(a$weight>0),], generation==i)[,6]))
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

		plot3d(x, y, z, col="black", box=FALSE, type="n", xlab="", ylab="", zlab="", zlim=c(0, max(a$generation)), ylim=c(0, max(y)))
		rgl.viewpoint(35, 1, 90)  #sets viewpoint for initial plot
		title3d(colnames(a)[param.position], col='red', pos=c(NA, 0, max(z))) 
		#text3d(x=min(x), y=mean(y), z=max(z), text="Density" col='blue') 
		#title3d("Starting States", col='purple', pos=c(NA, 0, max(which.gen))) 

		for (i in 1:max(s)){
			ngen<-(i+1)-1
			triangles<-triangulate(as(cbind(x[which(z==i)],y[which(z==i)]), "gpc.poly"))
			zfit<-predict(lm(z[which(z==i)] ~ x[which(z==i)] + y[which(z==i)]), newdata=data.frame(x=triangles[,1], y=triangles[,2]))
			opacity<-0.8*(ngen/length(v))
			rgl.material(color="black", alpha=opacity, lit=FALSE)
			triangles3d(cbind(triangles, zfit), col="red")
		}
	}		
		
	else {
		warning(paste("You are attempting to plot",colnames(a)[param.position],", which is a fixed value")) #return fixed val
		}
		
		
show.particles<-match.arg(arg=show.particles, choices=c("none", "weights", "distance"),several.ok=FALSE)
kept<-subset(a[which(a$id>0),])[,]	
reject<-subset(a[which(a$id<0),])[,]
short.kept<-subset(kept[which(kept$generation>1),])[,]

	if (show.particles=="none"){
		cat("currently not plotting particles.  To plot particles modify the show.particles= argument.")
		}
	else if (show.particles=="weights") {
				#for (i in 1:nrow(reject)) {	
		#	particles3d(x=reject[i, param.position], y=-1*(.05*(max(y)-min(y))), z=reject[i, 1], radius=.3, type="s", col="red", add=TRUE, aspect=TRUE)
		#}
		for (j in 1:nrow(kept)) {
			circle.size<-(kept[j, 6]/max(kept[,6]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))
			particles3d(x=kept[j, param.position], y=-1*(.05*(max(y)-min(y))), z=kept[j, 1], radius=circle.size, type="s", col="black", add=TRUE, aspect=TRUE)
		}	
	}		
	else if (show.particles=="distance") {
				#for (i in 1:nrow(reject)) {	
		#	particles3d(x=reject[i, param.position], y=-1*(.05*(max(y)-min(y))), z=reject[i, 1], radius=.3, type="s", col="red", add=TRUE, aspect=TRUE)
		#}
		for (j in 1:nrow(kept)) {
			circle.size<-(kept[j, 5]/max(kept[,5]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))
			particles3d(x=kept[j, param.position], y=-1*(.05*(max(y)-min(y))), z=kept[j, 1], radius=circle.size, type="s", col="black", add=TRUE, aspect=TRUE)
		}	
	}	
	if (plot.parent){
		for (k in 1:nrow(short.kept)) {
		prev.gen<-subset(kept[which(kept$generation==short.kept[k, 1]-1),])[,]  #works to retreive prev gen
		rgl.material(color="black", alpha=0.5)
		lines3d(x=c(short.kept[k, param.position], prev.gen[short.kept[k,]$parentid, param.position]), y=c(0, 0), z=c(short.kept[k, 1], short.kept[k, 1]-1))		}	
	}	
}














