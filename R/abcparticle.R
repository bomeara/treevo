#class for particles in ABC-PRC
setClass(
	Class="abcparticle",
	representation=representation(
		id="numeric",
		generation="numeric",
		weight="numeric",
		distance="numeric",
		startingStates="numeric",
		intrinsicStates="numeric",
		extrinsicStates="numeric"
	)
)

setGeneric("id", function(object) {
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

setGeneric("intrinsicStates", function(object) {
	standardGeneric("intrinsicStates")
})

setGeneric("intrinsicStates<-", function(object, value) {
	standardGeneric("intrinsicStates<-")
})


setMethod("intrinsicStates","abcparticle",function(object) {
	object@intrinsicStates
})

setReplaceMethod("intrinsicStates", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@intrinsicStates <- value
	object
}
)

setGeneric("extrinsicStates", function(object) {
	standardGeneric("extrinsicStates")
})

setGeneric("extrinsicStates<-", function(object, value) {
	standardGeneric("extrinsicStates<-")
})


setMethod("extrinsicStates","abcparticle",function(object) {
	object@extrinsicStates
})

setReplaceMethod("extrinsicStates", signature(object="abcparticle", value="numeric") ,
function(object, value) {
	object@extrinsicStates <- value
	object
}
)

setGeneric("initializeStatesFromMatrices",function(x, startingMatrix, intrinsicMatrix, ExtrinsicMatrix){standardGeneric("initializeStatesFromMatrices")})

setMethod("initializeStatesFromMatrices",signature="abcparticle",definition=function(x, startingMatrix, intrinsicMatrix, ExtrinsicMatrix) {
		x@startingStates<-rep(NA,dim(startingMatrix)[2])
		x@intrinsicValues<-rep(NA,dim(intrinsicMatrix)[2])
		x@extrinsicValues<-rep(NA,dim(extrinsicMatrix)[2])
		for (j in 1:length(startingStates)) {
			x@startingStates[j]=runif(n=1,min=min(startingMatrix[,j]),max=max(startingMatrix[,j]))
		}
		for (j in 1:length(intrinsicValues)) {
			x@intrinsicValues[j]=runif(n=1,min=min(intrinsicMatrix[,j]),max=max(intrinsicMatrix[,j]))
			
		}
		for (j in 1:length(extrinsicValues)) {
			x@extrinsicValues[j]=runif(n=1,min=min(extrinsicMatrix[,j]),max=max(extrinsicMatrix[,j]))
		}
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
	})
	
setGeneric("computeDistance",function(x, summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep){standardGeneric("computeDistance")})
	
setMethod("computeDistance",signature="abcparticle",definition=function(x, summaryFn, originalSummary, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	x@distance<-dist(matrix(c(summaryFn(phy, convertTaxonFrameToGeigerData(doSimulation(splits,intrinsicFn,extrinsicFn,x@startingStates,x@intrinsicStates,x@extrinsicStates,timeStep),phy)), originalSummary),nrow=2,byrow=T))[1]
	})

	