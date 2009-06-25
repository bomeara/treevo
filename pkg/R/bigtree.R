bigcontinuous<-function(firstlevel=1000,secondlevel=500,nchar=2) {
	library(geiger)
	library(multicore)
	tipsexamined<-0;
	print("starting")
	basalportion<-birthdeath.tree(b=1,d=0,taxa.stop=firstlevel)
	class(basalportion)<-"phylo"
	
	basalportion<-chronogram(basalportion,scale=1,minEdgeLength=1.0/(firstlevel*10))
	#basalportion<-evolve.phylo(basalportion,value=c(0,0),var=1)
		#print(dput(basalportion))	#class(basalportion)<-"phylo"
	derivedportion<-mclapply(1:firstlevel,maketreeandchars,ntax=secondlevel, parenttree= basalportion,mc.cores=3)
	#use runif for seed rather than geiger's default (time-based) just in case
	basalportion$tip.label<-paste("b", basalportion $tip.label,sep="")
	write.tree(basalportion,file="basal.tre")
	for (i in 1:firstlevel) {
		toadd<-derivedportion[[i]]
		class(toadd)<-"phylo"
		toadd$tip.label<-paste("taxon",as.character(tipsexamined+as.numeric(toadd$tip.label)),sep="") #birthdeath.tree results in trees with labels 1, 2, 3...
		write.tree(toadd,file=paste("derived",as.character(i),".tre"))
		tipsexamined<-secondlevel+tipsexamined
		basalportionTEMP<-bind.tree(x=basalportion,y=toadd,where=match(x=paste("b",i,sep=""),table=basalportion$tip.label))
		basalportion<-basalportionTEMP
		plot(basalportion)
		print(c(i,tipsexamined))
		}
		basalportion
}

maketreeandchars<-function(parenttaxon,ntax, parenttree) {
	print(parenttaxon)
	#ancestralstates<-as.numeric(parenttree$tip.character[parenttaxon,])
	#print(ancestralstates)
	newtree<-birthdeath.tree(b=1,d=0,taxa.stop=ntax,seed=round(runif(n=1,min=1,max=10000000)))
	class(newtree)<-"phylo"
	newtree<-chronogram(newtree,scale=1,minEdgeLength=1.0/(ntax*10))
	#a<-evolve.phylo(newtree,value= ancestralstates,var=1)
		#class(a)<-"phylo"

	#print(mean(as.numeric(a$tip.character)))
	newtree
	}
	
bigcontinuoussaving<-function(firstlevel=1000,secondlevel=500,nchar=2) {
	library(geiger)
	#library(multicore)
	tipsexamined<-0;
	print("starting")
	basalportion<-birthdeath.tree(b=1,d=0,taxa.stop=firstlevel)
	class(basalportion)<-"phylo"
	
	basalportion<-chronogram(basalportion,scale=1,minEdgeLength=1.0/(firstlevel*10))
	#basalportion<-evolve.phylo(basalportion,value=c(0,0),var=1)
		#print(dput(basalportion))	#class(basalportion)<-"phylo"
			basalportion$tip.label<-paste("b", basalportion $tip.label,sep="")
	write.tree(basalportion,file="basal.tre")

	#derivedportion<-mclapply(1:firstlevel, maketreesaving,ntax=secondlevel,mc.cores=3)
		derivedportion<-lapply(1:firstlevel, maketreesaving,ntax=secondlevel)

	#use runif for seed rather than geiger's default (time-based) just in case
}

maketreesaving<-function(parenttaxon,ntax) {
	print(parenttaxon)
	#ancestralstates<-as.numeric(parenttree$tip.character[parenttaxon,])
	#print(ancestralstates)
	newtree<-birthdeath.tree(b=1,d=0,taxa.stop=ntax,seed=round(runif(n=1,min=1,max=10000000)))
	class(newtree)<-"phylo"
	newtree<-chronogram(newtree,scale=1,minEdgeLength=1.0/(ntax*10))
	#a<-evolve.phylo(newtree,value= ancestralstates,var=1)
		#class(a)<-"phylo"
		newtree$tip.label<-paste("taxon",as.character((parenttaxon-1)*ntax +as.numeric(newtree$tip.label)),sep="") #birthdeath.tree results in trees with labels 1, 2, 3...

write.tree(newtree,file=paste("derived",as.character(parenttaxon),".tre",sep=""))
	#print(mean(as.numeric(a$tip.character)))
	newtree
	}