# function for getting prior matrix for doRun_prc
getPriorMatrix <- function(
    startingPriorsValues, 
    intrinsicPriorsValues, 
    extrinsicPriorsValues, 
    startingPriorsFns, 
    intrinsicPriorsFns, 
    extrinsicPriorsFns, 
    numberParametersTotal){
    #
    #
    #create PriorMatrix
    namesForPriorMatrix <- c()
    PriorMatrix <- matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow = 1, ncol = numberParametersTotal)
    for (a in 1:length(startingPriorsValues)){
        namesForPriorMatrix <- c(paste0("starting_", a, sep = ""))
        }
    #
    for (b in 1:length(intrinsicPriorsValues)){
        namesForPriorMatrix <- append(namesForPriorMatrix, paste0("intrinsic_", b, sep = ""))
        }
    #
    #message(extrinsicPriorsValues)
    for (c in 1:length(extrinsicPriorsValues)){
        namesForPriorMatrix  <- append(namesForPriorMatrix, paste0("extrinsic_", c, sep = ""))
        }
    #
    PriorMatrix <- rbind(
        PriorMatrix, 
        cbind(startingPriorsValues, 
            intrinsicPriorsValues, 
            extrinsicPriorsValues
            )
        )
    colnames(PriorMatrix) <- namesForPriorMatrix
    rownames(PriorMatrix) <- c("shape", "value1", "value2")
    return(PriorMatrix)
    }
