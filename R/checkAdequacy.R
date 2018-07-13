
#'

#*' @details

#*' @inheritParams

#*' @param

#*' @return

#*' @aliases

#*' @seealso

#*' @author David W. Bapst

#*' @references

#*' @examples

#*' @name
#*' @rdname
#*' @export


#############################################################################################
# notes from conversation with Smits (05-09-18)
# Me: so I'm doing approximate bayesian computation and the question is, what do I want to show to the reader
# Me: I want to show posterior parameter estimates from real data, and show that they are very different from parameter estimates made under other models, or under the same model but with simulated data, for scenarios with a small number of models
# what i do is make the same series of posterior predictive checks and demonstrate how your prefered model better recapitulates the data it was fit to
# like....how well do they all reproduce to ecdf or the density of the original data? if your model does a better job of doing that, then it is straight up a better model
# it also goes beyond parameter estimates and towards the model describing the data
# ECDF - empirical cumulative distribution function = the ranked order accumulation curve
# http://stat.ethz.ch/R-manual/R-devel/library/stats/html/ecdf.html
# bayes wants to describe more than just the expected value. it is greedy and wants to describe the whole posterior
# the posterior predictive distribution describes all data sets that are consistent with the model, given the original input information
# if the PPD doesn't look like the empirical data, then the model is not describing your data
# ecdf is a cool way of summarizing the entire dataset graphically
#####################################################################################
# from 06-21-18
# okay so the general sketch is particles from the posterior, simulate under this set of parameters N times, 
	# and compare the original ECDF for each parameter to the simulated
# my other idea:
# draw parameter estimates from some posterior particle, simulate under those parameters, 
	#  then test if 'true' generating parameters are actually within the 95% HDF of the simulated posterior
	# deals with how we don't really understand how adequate the models are for giving unbiased estimates of parameters

# checkAdequacy              # sixAnalysesTwoModels
	

checkAdequacy <- function(){
	}



#I mean, writing a function that just takes arguments: tree, params, etc. and returns results would be good

#a lot of this is (dataset) and six corresponding analyses

#fit model A, model B to real data, then simulate under model A, model B and fits both model A and B to both
#(where A is usually BM)


