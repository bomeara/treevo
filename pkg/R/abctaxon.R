setClass(
Class="abctaxon",
representation=representation(
id="numeric",
name="character",
states="numeric",
nextstates="numeric",
timeSinceSpeciation="numeric"
),
prototype(id=0,name="taxon0",timeSinceSpeciation=0)
)

setGeneric("id", function(object) {
	standardGeneric("id")
})

setGeneric("id<-", function(object, value) {
	standardGeneric("id<-")
})


setMethod("id","abctaxon",function(object) {
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


setMethod("name","abctaxon",function(object) {
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


setMethod("states","abctaxon",function(object) {
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


setMethod("nextstates","abctaxon",function(object) {
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


setMethod("timeSinceSpeciation","abctaxon",function(object) {
	object@timeSinceSpeciation
})

setReplaceMethod("timeSinceSpeciation", signature(object="abctaxon", value="numeric") ,
function(object, value) {
	object@timeSinceSpeciation <- value
	object
}
)

setGeneric("updateStates", function(object) {
	standardGeneric("updateStates")
})


setMethod("updateStates", signature(object='abctaxon'),
function(object) {
	object@states<-object@nextstates
	object
})
