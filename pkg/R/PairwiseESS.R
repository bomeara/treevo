setwd("/Users/BarbBanbury/Desktop/newton_BackUp/Newton-11.7.11/Models/Anolis-paper/TestTimeStep2")
#setwd("/Users/Barb/Desktop/TestTime/TestTimeStep2")
myFiles<-system("ls a*",intern=TRUE)
all.a<-vector("list", length(myFiles))
time<-vector()
timeSteps<-vector()

for (i in 1:length(myFiles)) {
	load(myFiles[i])
	all.a[[i]]<-a[[3]]
	time[[i]]<-a$time.per.gen[1]
	timeSteps[[i]]<-a$input.data[4]  #add back with new doRun

}

set2<-vector("list", 2)
set2[[1]]<-all.a[[1]]
set2[[2]]<-all.a[[2]]

set3<-vector("list", 3)
set3[[1]]<-all.a[[1]]
set3[[2]]<-all.a[[2]]
set3[[3]]<-all.a[[3]]


PairwiseESS<-function(particleDataFrame) {
	library(coda)
	# particleDataFrame can be single or a list of particleDataFrames (1:n)
	x<-particleDataFrame
	
	source("/Users/BarbBanbury/Desktop/treevo/pkg/R/ESSpairs.R")
	if(class(x)=="data.frame"){
		warning("ESS on a single run should be high; consider combining several runs")
		data1<-subset(x[which(x[,6]>0),], generation==max(x[,1])) 
		ESS<-effectiveSize(data1[,7:dim(data1)[2]])
	}
	
	if(class(x)=="list"){
		ESS<-vector("list")
		data1<-vector("list")
		pairs<-as.matrix(ESSpairs(length(x)))
		for (list in 1:length(x)) {
			data1[[list]]<-subset(x[[list]][which(x[[list]][,6]>0),], generation==max(x[[list]]$generation)) 
		}	
		for (combination in 1:dim(pairs)[2]){
			runNames<-c()
			data2<-c()
			runNames<-c()
			for (run in 1:dim(pairs)[1]){
				if (pairs[run, combination]==1){
					data2<-rbind(data2, data1[[run]])
					runNames<-append(runNames, run)
				}	
			}	
			ESS[[combination]]<-effectiveSize(data2[,7:dim(data2)[2]])
			names(ESS)[[combination]]<-paste(runNames, sep="", collapse=".")
		}
	}

	return(ESS)
}

PairwiseESS(set3)->l

PairwiseESS(set2)->l

effectiveSize(l)


			#all<-rbind(all, x[[list]])		





