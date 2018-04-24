# internal functions for figuring out which multicore approach to use, depending on what the user has installed among suggested packages


registerMulticoreEnv<-function(nCore){	
    hasDoParallel<-requireNamespace("doParallel", quietly = TRUE)
    hasDoMC<-requireNamespace("doMC", quietly = TRUE)			
	#
	if(hasDoMC){
		cluster<-doMC::registerDoMC(nCore)	#set number of cores back to 1
	}else{
		if(hasDoParallel){
			cluster<-doParallel::registerDoParallel(nCore)
		}else{
			stop("Argument multicore requires either suggested package 'doMC' or package 'doParallel' installed")
			}
		}
	#endFun<-function(){
	#	parallel::stopCluster(nCore)
	#	foreach::registerDoSEQ()
	#	}
	}	

setupMulticore<-function(multicore,nSim,coreLimit){
	# set cores to 1 as a placeholder
	cores<-1
	#
	if (multicore) {
		if (is.na(coreLimit)){
			cores<-c(nSim,parallel::detectCores()) #getDoParWorkers()
		}else{
			# don't allow it to take more cores than exist
			cores<-c(nSim,coreLimit,parallel::detectCores())
			}
		cores<-min(cores)
		}
	#
	if(cores>1){
		platform<-.Platform$GUI
		if(platform!="Rterm"){
			warning(paste0("Your platform appears to be ",platform,
				"\n Multicore processes are ideally done in Rterm, preferably in BATCH"))
			}
		registerMulticoreEnv(cores)
	}else{
		# if not multicore, need to setup registerDoSEQ so that foreach doesn't complain
		foreach::registerDoSEQ()
		}
	return(cores)
	}