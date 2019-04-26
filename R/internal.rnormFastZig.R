# an alternative replacement to rpgm::rpgm.rnorm
	# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	# 04-25-19

rnormFastZig <- function(n, mean, sd){
	# an alternative replacement to rpgm::rpgm.rnorm
		# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	(RcppZiggurat::znorm(n) * sd) + mean
	}