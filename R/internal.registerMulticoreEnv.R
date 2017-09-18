# internal function for figuring out which multicore approach to use, depending on what the user has installed among suggested packages


registerMulticoreEnv<-function(nCore){	
    hasDoParallel<-requireNamespace("doParallel", quietly = TRUE)
    hasDoMC<-requireNamespace("doMC", quietly = TRUE)			 
	#
	if(hasDoMC){
		doMC::registerDoMC(nCore)	#set number of cores back to 1
	}else{
		if(hasDoParallel){
			doParallel::registerDoParallel(nCore)
		}else{
			stop("Argument multicore requires either suggested package 'doMC' or package 'doParallel' installed")
			}
		}
	}	