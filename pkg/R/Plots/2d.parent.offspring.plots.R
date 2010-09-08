parent.offspring.plots<-function(x){

nparams<-dim(a)[2]-6
nb<-2*nparams
nf<-layout(matrix(1:nb, nrow=2, byrow=TRUE), respect=TRUE)
layout.show(nf)


for (param in 1:nparams) {
	param.position<-param+6
			
	plot(a[,param.position], a$generation, xlab=colnames(a)[param.position], ylab="Generation", type="n")
	title("size as measure of distance")
	kept<-subset(a[which(a$id>0),])[,]	
	reject<-subset(a[which(a$id<0),])[,]
	short.kept<-subset(kept[which(kept$generation>1),])[,]

	#symbols(x, y, circles=value)

	for (i in 1:nrow(reject)) {
		#circle.size<-(kept[j, 5]/max(kept[,5]))*mean(a[,param.position])
		circle.size<-(reject[i, 5]/max(reject[,5]))*(0.05*(max(a[,param.position])-min(a[,param.position])))
		symbols(x=reject[i, param.position], y=reject[i, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="gray")
	}	
	
	for (j in 1:nrow(kept)) {
		#circle.size<-(kept[j, 5]/max(kept[,5]))*mean(a[,param.position])
		circle.size<-(kept[j, 5]/max(kept[,5]))*(0.05*(max(kept[,param.position])-min(kept[,param.position])))
		symbols(x=kept[j, param.position], y=kept[j, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="black")
	}		
	for (k in 1:nrow(short.kept)) {
		#lines(c(x1, x2), c(y1, y2))
		#lines(c(param.offspring, param.parent), c(generation.offspring, generation.parent))
		prev.gen<-subset(kept[which(kept$generation==short.kept[k, 1]-1),])[,]  #works to retreive prev gen
		#parent.param<-subset(prev.gen[which(kept$generation==short.kept[k, 1]-1),])[,]
		lines(c(short.kept[k, param.position], prev.gen[short.kept[k,]$parentid, param.position]), c(short.kept[k, 1], short.kept[k, 1]-1))
	}	
}		
	
	
##-----For particle weights

for (param in 1:nparams) {
	param.position<-param+6

	plot(a[,param.position], a$generation, xlab=colnames(a)[param.position], ylab="Generation", type="n")
	title("size as measure of particle weights")
	kept<-subset(a[which(a$id>0),])[,]	
	reject<-subset(a[which(a$id<0),])[,]
	short.kept<-subset(kept[which(kept$generation>1),])[,]

	#symbols(x, y, circles=value)

	for (i in 1:nrow(reject)) {
		points(x=reject[i, param.position], y=reject[i, 1], col="gray", pch=8)

	}	
	
	for (j in 1:nrow(kept)) {
		#circle.size<-(kept[j, 6]/max(kept[,6]))/mean(a[,param.position])
		circle.size<-(kept[j, 6]/max(kept[,6]))*(0.05*(max(kept[,param.position])-min(kept[,param.position])))
		symbols(x=kept[j, param.position], y=kept[j, 1], circles=circle.size, inches=FALSE, add=TRUE, fg="black")
	}		
	for (k in 1:nrow(short.kept)) {
		#lines(c(x1, x2), c(y1, y2))
		#lines(c(param.offspring, param.parent), c(generation.offspring, generation.parent))
		prev.gen<-subset(kept[which(kept$generation==short.kept[k, 1]-1),])[,]  #works to retreive prev gen
		#parent.param<-subset(prev.gen[which(kept$generation==short.kept[k, 1]-1),])[,]
		lines(c(short.kept[k, param.position], prev.gen[short.kept[k,]$parentid, param.position]), c(short.kept[k, 1], short.kept[k, 1]-1))
	}	
}	

}		
