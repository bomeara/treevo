#class for particles in ABC-PRC
setClass(
	Class="abcparticle",
	representation=representation(
		id="numeric",
		generation="numeric",
		weight="numeric",
		distance="numeric",
		startingStates="numeric",
		intrinsicValues="numeric",
		extrinsicValues="numeric",
		kernel="numeric"
	)
)

setGeneric("id",function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id","abcparticle",function(object) {
	object@id
})

setReplaceMethod("id", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@id <- value
	object
}
)

setGeneric("setId", function(object, value) {
	standardGeneric("setId")
})

setMethod("setId","abcparticle",function(object,value) {
	object@id<-value
	object
})

setGeneric("getId", function(object) {
	standardGeneric("getId")
})

setMethod("getId","abcparticle",function(object) {
	object@id
})

setGeneric("setWeight", function(object, value) {
	standardGeneric("setWeight")
})

setMethod("setWeight","abcparticle",function(object,value) {
	object@weight<-value
	object
})

setGeneric("getWeight", function(object) {
	standardGeneric("getWeight")
})

setMethod("getWeight","abcparticle",function(object) {
	object@weight
})

setGeneric("generation", function(object) {
	standardGeneric("generation")
})

setGeneric("generation<-", function(object, value) {
	standardGeneric("generation<-")
})


setMethod("generation","abcparticle",function(object) {
	object@generation
})

setReplaceMethod("generation", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@generation <- value
	object
}
)

setGeneric("weight", function(object) {
	standardGeneric("weight")
})

setGeneric("weight<-", function(object, value) {
	standardGeneric("weight<-")
})


setMethod("weight","abcparticle",function(object) {
	object@weight
})

setReplaceMethod("weight", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@weight <- value
	object
}
)

setGeneric("distance", function(object) {
	standardGeneric("distance")
})

setGeneric("distance<-", function(object, value) {
	standardGeneric("distance<-")
})


setMethod("distance","abcparticle",function(object) {
	object@distance
})

setReplaceMethod("distance", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@distance <- value
	object
}
)

setGeneric("setDistance", function(object, value) {
	standardGeneric("setDistance")
})

setMethod("setDistance","abcparticle",function(object,value) {
	object@distance<-value
	object
})

setGeneric("kernel", function(object) {
	standardGeneric("kernel")
})

setGeneric("kernel<-", function(object, value) {
	standardGeneric("kernel<-")
})


setMethod("kernel","abcparticle",function(object) {
	object@kernel
})

setReplaceMethod("kernel", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@kernel <- value
	object
}
)


setGeneric("startingStates", function(object) {
	standardGeneric("startingStates")
})

setGeneric("startingStates<-", function(object, value) {
	standardGeneric("startingStates<-")
})


setMethod("startingStates","abcparticle",function(object) {
	object@startingStates
})

setReplaceMethod("startingStates", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@startingStates <- value
	object
}
)

setGeneric("intrinsicValues", function(object) {
	standardGeneric("intrinsicValues")
})

setGeneric("intrinsicValues<-", function(object, value) {
	standardGeneric("intrinsicValues<-")
})


setMethod("intrinsicValues","abcparticle",function(object) {
	object@intrinsicValues
})

setReplaceMethod("intrinsicValues", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@intrinsicValues <- value
	object
}
)

setGeneric("extrinsicValues", function(object) {
	standardGeneric("extrinsicValues")
})

setGeneric("extrinsicValues<-", function(object, value) {
	standardGeneric("extrinsicValues<-")
})


setMethod("extrinsicValues","abcparticle",function(object) {
	object@extrinsicValues
})

setReplaceMethod("extrinsicValues", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@extrinsicValues <- value
	object
}
)

setGeneric("initializeStatesFromMatrices",function(x, startingMatrix, intrinsicMatrix, extrinsicMatrix){standardGeneric("initializeStatesFromMatrices")})

setMethod("initializeStatesFromMatrices",signature="abcparticle",definition=function(x, startingMatrix, intrinsicMatrix, extrinsicMatrix) {
		x@startingStates<-vector(mode="numeric",length=dim(startingMatrix)[2])
		x@intrinsicValues<-vector(mode="numeric",length=dim(intrinsicMatrix)[2])
		x@extrinsicValues <-vector(mode="numeric",length=dim(extrinsicMatrix)[2])
		for (j in 1:length(startingStates)) {
			x@startingStates[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
		}
		for (j in 1:length(intrinsicValues)) {
			x@intrinsicValues[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
			
		}
		for (j in 1:length(extrinsicValues)) {
			x@extrinsicValues[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
		}
		#print(x)
		x
		}
)

setGeneric("mutateStates",function(x, startingMatrix, intrinsicMatrix, extrinsicMatrix, standardDevFactor){standardGeneric("mutateStates")})

setMethod("mutateStates",signature="abcparticle",definition=function(x, startingMatrix, intrinsicMatrix, extrinsicMatrix, standardDevFactor) {
	#dput(x)
	#dput(x@startingStates)
	#typeof(x)
	#typeof(x@startingStates)
		replacementVector<-rep(NA, length(x@startingStates))
		for (j in 1:length(x@startingStates)) {
			newvalue<-Inf
			meantouse=x@startingStates[j]
			sdtouse=standardDevFactor*(max(startingMatrix[,j])-min(startingMatrix[,j]))
			while (newvalue<min(startingMatrix[,j]) || newvalue>max(startingMatrix[,j])) {
				newvalue<-rnorm(n=1,mean= meantouse ,sd= sdtouse)
				replacementVector[j]=newvalue
			}
		}
		x@startingStates<-replacementVector
		
		replacementVector<-rep(NA, length(x@intrinsicValues))
		for (j in 1:length(x@intrinsicValues)) {
			newvalue<-Inf
			meantouse=x@intrinsicValues[j]
			sdtouse=standardDevFactor*(max(intrinsicMatrix[,j])-min(intrinsicMatrix[,j]))
			while (newvalue<min(intrinsicMatrix[,j]) || newvalue>max(intrinsicMatrix[,j])) {
				newvalue<-rnorm(n=1,mean= meantouse ,sd= sdtouse)
				replacementVector[j]=newvalue
			}
		}
		x@intrinsicValues <-replacementVector

		replacementVector<-rep(NA, length(x@extrinsicValues))
		for (j in 1:length(x@extrinsicValues)) {
			newvalue<-Inf
			meantouse=x@extrinsicValues[j]
			sdtouse=standardDevFactor*(max(extrinsicMatrix[,j])-min(extrinsicMatrix[,j]))
			while (newvalue<min(extrinsicMatrix[,j]) || newvalue>max(extrinsicMatrix[,j])) {
				newvalue<-rnorm(n=1,mean= meantouse ,sd= sdtouse)
				replacementVector[j]=newvalue
			}
		}
		x@extrinsicValues <-replacementVector
		
		#kernel(x)<-lnTransitionProb
		x
		}
)

setGeneric("transformStates",function(x, sdVector ){standardGeneric("transformStates")})

setMethod("transformStates",signature="abcparticle",definition=function(x, sdVector) {
	positioncount<-1
	for (j in 1:length(startingStates)) {
			x@startingStates[j]=rnorm(n=1,mean=startingStates[j],sd=sdVector[positioncount])
			positioncount<-positioncount+1
	}
	for (j in 1:length(intrinsicValues)) {
		x@intrinsicValues[j]=rnorm(n=1,mean=intrinsicValues[j],sd=sdVector[positioncount])
		positioncount<-positioncount+1
		
		}
		for (j in 1:length(extrinsicValues)) {
		x@extrinsicValues[j]=rnorm(n=1,mean=extrinsicValues[j],sd=sdVector[positioncount])
			positioncount<-positioncount+1
		}
		x
	})
	
setGeneric("computeABCDistance",function(x, summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep){standardGeneric("computeABCDistance")})
	
setMethod("computeABCDistance",signature="abcparticle",definition=function(x, summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	x@distance<-dist(matrix(c(summaryFn(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,x@startingStates,x@intrinsicValues,x@extrinsicValues,timeStep),phy)), originalSummary),nrow=2,byrow=T))[1]
	#print("computeABCDistance")
	#print(x)
	x
	})

setGeneric("simulateTips",function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep){standardGeneric("simulateTips")})
	
setMethod("simulateTips",signature="abcparticle",definition=function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,x@startingStates,x@intrinsicValues,x@extrinsicValues,timeStep),phy)
	return(newtips)
	})	
	

	