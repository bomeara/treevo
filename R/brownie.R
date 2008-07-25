browniecensor<-function(phyd,chosenchar,cladeA,reps=0) {
	#cladeA is a vector of taxa
	
	objectA<-subset(phyd,mrca=cladeA)
	objectB<-subset(phyd,tips.exclude=descendants(phyd,MRCA(phyd,getnodes(phyd,cladeA)),which=c("tips")))
	
	dataA<-tdata(objectA,which=c("tip"))[,chosenchar]
	dataB<-tdata(objectB,which=c("tip"))[,chosenchar]
	dataCombined<-c(dataA,dataB)
	
	vcvA<-brownievcv(objectA)
	vcvB<-brownievcv(objectB)
	vcvCombined<-rbind(cbind(vcvA,array(0,dim=c(nTips(objectA),nTips(objectB)))),cbind(array(0,dim=c(nTips(objectB),nTips(objectA))),vcvB)) 
	#vcvCombined=[  vcvA  zeros ]
	#			 [ zeros  vcvB  ]
	
	ancstatevectorA<-rep(brownieanc(dataA,vcvA),times=nTips(objectA))
	ancstatevectorB<-rep(brownieanc(dataB,vcvB),times=nTips(objectB))
	ancstatevectorCombined<-c(ancstatevectorA,ancstatevectorB)
	
	sigmasquaredA<-browniesigmasquared(dataA, ancstatevectorA,vcvA)
	sigmasquaredB<-browniesigmasquared(dataB, ancstatevectorB,vcvB)
	sigmasquaredCombined<-browniesigmasquared(dataCombined, ancstatevectorCombined,vcvCombined)
	
	neglnLA<-brownieneglnL(dataA, ancstatevectorA,vcvA, sigmasquaredA)
	neglnLB<-brownieneglnL(dataB, ancstatevectorB,vcvB, sigmasquaredB)
	neglnL_double<-neglnLA+neglnLB
	neglnL_single<-brownieneglnL(dataCombined, ancstatevectorCombined,vcvCombined, sigmasquaredCombined)
	
	aic_both<-c(brownieaic(neglnL= neglnL_single,K=3),brownieaic(neglnL=neglnL_double,K=4))
	aic_both<-aic_both-min(aic_both)
	
	aicc_both<-c(brownieaicc(neglnL=neglnL_single,K=3,N=nTips(phyd)),brownieaicc(neglnL=neglnL_double,K=4,N=nTips(phyd)))
	aicc_both<-aicc_both-min(aicc_both)
	
	result<-c(ancstatevectorA[1],sigmasquaredA,neglnLA,ancstatevectorB[1],sigmasquaredB,neglnLB,sigmasquaredCombined,neglnL_single,aic_both,aicc_both)
	names(result)<-c("anc A","sigma-squared A", "negLnL A","anc B","sigma-squared B", "negLnL B","sigma-squared 1-rate","negLnL 1-rate","AIC 1-rate", "AIC 2-rate","AICc 1-rate","AICc 2-rate")
	result
}

brownieaic<-function(neglnL,K) {
	aic<-2*(neglnL+K)
	aic
}
	
brownieaicc<-function(neglnL,K,N) {
	aicc<-2*(neglnL+K)+2*K*(K+1)/(N-K-1)
	aicc
}

browniesigmasquared<-function(data,ancstatevector,vcv) {
	Ntax<-dim(vcv)[1]
	invvcv<-chol2inv(chol(vcv))
	sigmasquared=((data-ancstatevector)%*%invvcv%*%(data-ancstatevector))/Ntax
	sigmasquared
}

brownieanc<-function(data,vcv) {
	Ntax<-dim(vcv)[1]
	invvcv<-chol2inv(chol(vcv))
	ones<-rep(1,times=Ntax)
	secondpart<-ones%*%invvcv%*%data
	firstpart<-ones%*%invvcv%*%ones
	ancstate<-chol2inv(chol(firstpart))%*%secondpart
	ancstate
}

brownieneglnL<-function(data, ancstatevector,vcv,sigmasquared) {
	V<-sigmasquared[1]*vcv
	Ntax<-dim(V)[1]
	invV<-chol2inv(chol(V))
	detV<-det(V)
	neglnL<-(-1)*((-0.5*(data-ancstatevector)%*%invV%*%(data-ancstatevector))-log(sqrt(((2*pi)^Ntax)*detV)))
	neglnL
}

brownievcv<-function(tree) {
	#takes a phylo4 or phylo4d object
	treedataframe<-as(tree,"data.frame")
	tiplabels<-labels(tree,which=c("tip"))
	vcv<-array(0,dim=c(nTips(tree),nTips(tree)))
	for(leftid in 1:nTips(tree)) {
		for(topid in leftid:nTips(tree)) {
			chosennode<-MRCA(tree,tiplabels[leftid],tiplabels[topid])
			if (leftid==topid) {
				chosennode<-getnodes(tree,tiplabels[leftid]) #this is to fix an error in phylobase: MRCA of node with itself is not itself
			}
			edgelength<-0
			while(!is.na(subset(treedataframe, treedataframe$node==chosennode)$branch.length)) {
				edgelength<-edgelength+subset(treedataframe, treedataframe$node==chosennode)$branch.length
				chosennode<-ancestor(tree,chosennode)
			}
			vcv[leftid,topid]<-edgelength
			vcv[topid,leftid]<-edgelength
		}
	}
	vcv
}