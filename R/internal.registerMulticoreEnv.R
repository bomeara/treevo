# internal functions for figuring out which multicore approach to use, depending on what the user has installed among suggested packages


registerMulticoreEnv<-function(nCore){	
    hasDoParallel<-requireNamespace("doParallel", quietly = TRUE)
    hasDoMC<-requireNamespace("doMC", quietly = TRUE)			
	#
	if(hasDoMC){
		cluster<-doMC::registerDoMC(nCore)
		attr(cluster,"doParallel")<-FALSE
	}else{
		if(hasDoParallel){
			cluster<-doParallel::registerDoParallel(nCore)
			attr(cluster,"doParallel")<-TRUE
		}else{
			stop("Argument multicore requires either suggested package 'doMC' or package 'doParallel' installed")
			}
		}
	return(cluster)
	}	

setupMulticore<-function(multicore,nSim,coreLimit){
	# set nCores to 1 as a placeholder
	nCores<-1
	#
	if (multicore) {
		if (is.na(coreLimit)){
			nCores<-c(nSim,parallel::detectCores()) #getDoParWorkers()
		}else{
			# count number of cores, subtract one for system processes
			nCoresExist<-parallel::detectCores()-1
			# don't allow it to take more nCores than exist
			nCores<-c(nSim,coreLimit,nCoresExist)
			}
		nCores<-min(nCores)
		}
	#
	if(nCores>1){
		platform<-.Platform$GUI
		if(!identical(platform,"RTerm")){
			warning(paste0("Your platform appears to be ",platform,
				"\n Multicore processes are ideally done in RTerm, preferably via 'R CMD BATCH'"))
			}
		cluster<-registerMulticoreEnv(nCores)
		attr(cluster,"multicore")<-TRUE
	}else{
		# if not multicore, need to setup registerDoSEQ so that foreach doesn't complain
		foreach::registerDoSEQ()
		cluster<-list()
		attr(cluster,"multicore")<-FALSE
		}
	attr(cluster,"nCores")<-nCores
	return(cluster)
	}
	
stopMulticore<-function(cluster){
	if(attr(cluster,"multicore")){
		#parallel::stopCluster(cluster)
		if(attr(cluster,"doParallel")){	
			doParallel::stopImplicitCluster()
			}
		}
	#
    foreach::registerDoSEQ()
	}
	
	
	#endFun<-function(){
	#	parallel::stopCluster(nCore)
	#	foreach::registerDoSEQ()
	#	}
