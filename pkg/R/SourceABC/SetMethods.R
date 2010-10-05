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

setGeneric("id", function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id", "abcparticle", function(object) {
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

setMethod("setId", "abcparticle", function(object, value) {
	object@id<-value
	object
})

setGeneric("getId", function(object) {
	standardGeneric("getId")
})

setMethod("getId", "abcparticle", function(object) {
	object@id
})

setGeneric("setWeight", function(object, value) {
	standardGeneric("setWeight")
})

setMethod("setWeight", "abcparticle", function(object, value) {
	object@weight<-value
	object
})

setGeneric("getWeight", function(object) {
	standardGeneric("getWeight")
})

setMethod("getWeight", "abcparticle", function(object) {
	object@weight
})

setGeneric("generation", function(object) {
	standardGeneric("generation")
})

setGeneric("generation<-", function(object, value) {
	standardGeneric("generation<-")
})


setMethod("generation", "abcparticle", function(object) {
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


setMethod("weight", "abcparticle", function(object) {
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


setMethod("distance", "abcparticle", function(object) {
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

setMethod("setDistance", "abcparticle", function(object, value) {
	object@distance<-value
	object
})

setGeneric("kernel", function(object) {
	standardGeneric("kernel")
})

setGeneric("kernel<-", function(object, value) {
	standardGeneric("kernel<-")
})


setMethod("kernel", "abcparticle", function(object) {
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


setMethod("startingStates", "abcparticle", function(object) {
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


setMethod("intrinsicValues", "abcparticle", function(object) {
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


setMethod("extrinsicValues", "abcparticle", function(object) {
	object@extrinsicValues
})

setReplaceMethod("extrinsicValues", signature(object="abcparticle", value="numeric") , 
function(object, value) {
	object@extrinsicValues <- value
	object
}
)

setGeneric("initializeStatesFromMatrices", function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns){standardGeneric("initializeStatesFromMatrices")})

setMethod("initializeStatesFromMatrices", signature="abcparticle", definition=function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns) {
		#cat("startingPriorsValues\n")
		#print(startingPriorsValues)
		#cat("\nintrinsicPriorsValues\n")
		#print(intrinsicPriorsValues)
		#cat("\nextrinsicPriorsValues\n")
		#print(extrinsicPriorsValues)
		#cat("\ndim(startingPriorsValues)[2]=", dim(startingPriorsValues)[2], " dim(intrinsicPriorsValues)[2]=", dim(intrinsicPriorsValues)[2], " dim(extrinsicPriorsValues)[2])=", dim(extrinsicPriorsValues)[2], "\n")
		x@startingStates<-vector(mode="numeric", length=dim(startingPriorsValues)[2])
		x@intrinsicValues<-vector(mode="numeric", length=dim(intrinsicPriorsValues)[2])
		x@extrinsicValues <-vector(mode="numeric", length=dim(extrinsicPriorsValues)[2])
		#cat("\n initializeStatesFromMatrices\n")
		for (j in 1:dim(startingPriorsValues)[2]) {
			x@startingStates[j]=pullFromPrior(priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
			#cat("starting states j=", j, " x@startingStates[j] =", x@startingStates[j], "\n")
		}
		for (j in 1:dim(intrinsicPriorsValues)[2]) {
			x@intrinsicValues[j]=pullFromPrior(priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
			#cat("intrinsicValues j=", j, " x@ intrinsicValues[j] =", x@intrinsicValues[j], "\n")

		}
		for (j in 1:dim(extrinsicPriorsValues)[2]) {
			x@extrinsicValues[j]=pullFromPrior(priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
			#cat("extrinsicValues j=", j, " x@ extrinsicValues[j] =", x@extrinsicValues[j], "\n")
		}
		#print(x)
		x
		}
)



setGeneric("mutateStates", function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor){standardGeneric("mutateStates")})

setMethod("mutateStates", signature="abcparticle", definition=function(x, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor) {
	#dput(x)
	#dput(x@startingStates)
	#typeof(x)
	#typeof(x@startingStates)
		replacementVector<-rep(NA, length(x@startingStates))
		for (j in 1:length(x@startingStates)) {
			replacementVector[j]<-mutateState(startingState=x@startingStates[j], standardDevFactor=standardDevFactor, priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
		}
		x@startingStates<-replacementVector
		
		replacementVector<-rep(NA, length(x@intrinsicValues))
		for (j in 1:length(x@intrinsicValues)) {
			replacementVector[j]<-mutateState(startingState=x@intrinsicValues[j], standardDevFactor=standardDevFactor, priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
		}
		x@intrinsicValues<-replacementVector

		replacementVector<-rep(NA, length(x@extrinsicValues))
		for (j in 1:length(x@extrinsicValues)) {
			replacementVector[j]<-mutateState(startingState=x@extrinsicValues[j], standardDevFactor=standardDevFactor, priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
		}
		x@extrinsicValues<-replacementVector
		
		#kernel(x)<-lnTransitionProb
		x
		}
)


	

setGeneric("simulateTips", function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep){standardGeneric("simulateTips")})
	
setMethod("simulateTips", signature="abcparticle", definition=function(x, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, x@startingStates, x@intrinsicValues, x@extrinsicValues, timeStep), phy)
	return(newtips)
	})	
	
setClass(
Class="abctaxon", 
representation=representation(
id="numeric", 
name="character", 
states="numeric", 
nextstates="numeric", 
timeSinceSpeciation="numeric"
), 
prototype(id=0, name="taxon0", timeSinceSpeciation=0)
)

setGeneric("id", function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id", "abctaxon", function(object) {
	object@id
})

setReplaceMethod("id", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@id <- value
	object
}
)

setGeneric("name", function(object) {
	standardGeneric("name")
})

setGeneric("name<-", function(object, value) {
	standardGeneric("name<-")
})


setMethod("name", "abctaxon", function(object) {
	object@name
})

setReplaceMethod("name", signature(object="abctaxon", value="character") , 
function(object, value) {
	object@name <- value
	object
}
)



setGeneric("states", function(object) {
	standardGeneric("states")
})

setGeneric("states<-", function(object, value) {
	standardGeneric("states<-")
})


setMethod("states", "abctaxon", function(object) {
	object@states
})

setReplaceMethod("states", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@states <- value
	object
}
)



setGeneric("nextstates", function(object) {
	standardGeneric("nextstates")
})

setGeneric("nextstates<-", function(object, value) {
	standardGeneric("nextstates<-")
})


setMethod("nextstates", "abctaxon", function(object) {
	object@nextstates
})

setReplaceMethod("nextstates", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@nextstates <- value
	object
}
)


setGeneric("timeSinceSpeciation", function(object) {
	standardGeneric("timeSinceSpeciation")
})

setGeneric("timeSinceSpeciation<-", function(object, value) {
	standardGeneric("timeSinceSpeciation<-")
})


setMethod("timeSinceSpeciation", "abctaxon", function(object) {
	object@timeSinceSpeciation
})

setReplaceMethod("timeSinceSpeciation", signature(object="abctaxon", value="numeric") , 
function(object, value) {
	object@timeSinceSpeciation <- value
	object
}
)
