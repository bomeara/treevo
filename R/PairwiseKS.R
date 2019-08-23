#' Pairwise Kolmogorov-Smirnov test
#' 
#' This function calculates Kolmogorov-Smirnov on results.
#' 
#' 

#' @param particleDataFrameList An object of type list,
#' composed of particleDataFrames from separate analyses.

#' @return Returns a matrix with Kolmogorov-Smirnov values of all pairwise runs.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#' 
#' data(simRunExample)
#' 
#' pdfList <- list(
#'    Brownian = resultsBMExample[[1]]$particleDataFrame,
#'    Bounded = resultsBoundExample[[1]]$particleDataFrame
#'    )
#' 
#' pairwiseKS(particleDataFrameList = pdfList)


#' @name pairwiseKS
#' @rdname pairwiseKS
#' @export
pairwiseKS <- function(particleDataFrameList) {
    #Combine doRun$particleDataFrame results and performs Kolmogorov-Smirnov tests
    # particleDataFrame should be a list of particleDataFrames (1:n)

    #source("/Users/Barb/Desktop/treevo/pkg/R/pairings.R")
        
    if(class(particleDataFrameList) == "data.frame" | class(particleDataFrameList) != "list" ){
                    # | length(particleDataFrameList) != 2
        warning("KS test requires a list, composed of particleDataFrameList objects from two ABC runs")
    }
    
    x <- list()
    for (list in 1:length(particleDataFrameList)) {
        nonZeroFrame <- particleDataFrameList[[list]][which(particleDataFrameList[[list]][, 6]>0), ]
        x[[list]] <- subset(nonZeroFrame, nonZeroFrame$generation == max(nonZeroFrame$generation))
        }    

    KSMatrixList <- vector("list", dim(x[[1]])[2]-6)
    names(KSMatrixList) <- names(x[[1]][7:dim(x[[1]])[2]])
    for (m in 1:length(KSMatrixList)){
        param <- 6+m
        KSMatrix <- matrix(nrow = length(x), ncol = length(x))
        if(is.null(names(particleDataFrameList))){
            colnames(KSMatrix) <- rownames(KSMatrix) <- paste0("run", 1:length(x))
        }else{
            colnames(KSMatrix) <- rownames(KSMatrix) <- names(particleDataFrameList)
            }
        #colnames(KSMatrix) <- paste0("run", 1:length(x))
        for (i in 1:dim(KSMatrix)[1]){
            for (j in 1:dim(KSMatrix)[2]){
                KSMatrix[i, j] <- suppressWarnings(round(ks.test(x[[i]][, param],
                    x[[j]][, param])$p.value, digits = 4)) #make diag
                }
            }
        KSMatrixList[[m]] <- KSMatrix
        }

    return(KSMatrixList)
}
