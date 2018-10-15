#   Pull value from prior
#
#   This function retreives random value from variable prior functions
#
#
#   @param priorValues Matrix with ncol = number of states (characters) at root
#   and nrow = 2 (two parameters to pass to prior distribution)
#
#   @param priorFn Prior shape
#
#   @return Returns mean
#
#   @author Brian O'Meara and Barb Banbury
#
#	@references O'Meara and Banbury, unpublished


pullFromPrior <- function(priorValues, priorFn) {
    #fixed, uniform, normal, lognormal, gamma, exponential
    x <- NA
    priorFn <- match.arg(arg = priorFn,
		choices = c("fixed", "uniform", "normal",
			"lognormal", "gamma", "exponential"),
		several.ok = FALSE)
    if (priorFn == "fixed") {
		if(length(priorValues) != 1){
			stop(" the 'fixed' prior requires one parameter, not more, not less")
			}
		x <- priorValues[1]
    }else{
		if (priorFn == "uniform") {
			if(length(priorValues) != 2){
				stop(" the 'uniform' prior requires two parameters, not more, not less")
				}
			x <- runif(n = 1, min = min(priorValues), max = max(priorValues))
		}else{
			if (priorFn == "normal") {
				if(length(priorValues) != 2){
					stop(" the 'normal' prior requires two parameters, not more, not less")
					}
				x <- rnorm(n = 1, mean = priorValues[1], sd = priorValues[2])
			}else{
				if (priorFn == "lognormal") {
					if(length(priorValues) != 2){
						stop(" the 'lognormal' prior requires two parameters, not more, not less")
						}					
					x <- rlnorm(n = 1, meanlog = priorValues[1], sdlog = priorValues[2])
				}else{
					if (priorFn == "gamma") {
						if(length(priorValues) != 2){
							stop(" the 'gamma' prior requires two parameters, not more, not less")
							}					
						x <- rgamma(n = 1, shape = priorValues[1], scale = priorValues[2])
					}else{
						if (priorFn == "exponential") {
							if(length(priorValues) != 1){
								stop(" the 'exponential' prior requires one parameter, not more, not less")
								}
							x <- rexp(n = 1, rate = priorValues[1])
						}else{
							stop(priorFn, " was not a recognized prior function")
							}
						}
					}
				}
			}
		}
    return(x)
	}
