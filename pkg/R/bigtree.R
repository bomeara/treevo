bigcontinuous<-function(firstlevel=1000,secondlevel=500,nchar=2) {
	basalportion<-evolve.phylo(birthdeath.tree(b=1,d=0,taxa.stop=firstlevel),value=c(0,0),var=1)
	derivedportion<-lapply(1:firstlevel,maketreeandchars,ntax=secondlevel, parenttree= basalportion)
	#use runif for seed rather than geiger's default (time-based) just in case
	
}

maketreeandchars<-function(parenttaxon,ntax, parenttree) {
	print(parenttaxon)
	ancestralstates<-as.numeric(parenttree$tip.character[parenttaxon,])
	print(ancestralstates)
	a<-evolve.phylo(birthdeath.tree(b=1,d=0,taxa.stop=ntax,seed=round(runif(n=1,min=1,max=10000000))),value= ancestralstates,var=1)
	print(mean(as.numeric(a$tip.character)))
	a
	}