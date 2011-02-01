ABCplots<-function(a) {
gen1<-a[which(a$generation==1),]
nParticles<-dim(gen1[which(gen1$weight>0),])[1]
nparams<-dim(a)[2]-6
nf<-layout(matrix(1:nparams, nrow=1), respect=TRUE)
#layout.show(nf)


	for (param in 1:nparams){
		param.position<-param+6
		v<-vector("list", max(a$generation))
				for (i in 1:max(a$generation)){
					v[[i]]<-density(subset(a[which(a$weight>0),], generation==i)[,param.position], weights=nParticles*subset(a[which(a$weight>0),], generation==i)[,6]/sum(nParticles*subset(a[which(a$weight>0),], generation==i)[,6]))
				}	
				#call v should return the Call: from all density functions for nStepsPRC
				#call v[[1]]$x should call the x var
				
				maxy=0
				for (i in 1:max(a$generation)){
					maxy=max(maxy, max(v[[i]]$y))
				}
			
				plot(v[[1]]$x, v[[1]]$y, xlab=(names(a[param.position])), ylim=c(0, maxy), ylab="density", type="n")
				subA<-subset(a[which(a$weight>0),])
			 
			
				if (range(subA[,param.position])!=0) {
					for (i in 1:max(a$generation)){
						ngen<-(i+1)-1
						#shade<-rgb(0, 0, 0, 0.5*(length(v)-i)/length(v))
						polygon(v[[i]], col=rgb(0, 0, 0, 0.8*(ngen/length(v))))  #make darker every time by 0.5*ngen i/max for proportion of gens
					}
				}
				else {
					lines(x=rep((subA[,param.position])[1],2),y=c(0, maxy*.9))
				}
		}	
}
