# an alternative replacement to rpgm::rpgm.rnorm
	# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	# 04-25-19

rnormFastZig <- function(n, mean = 0, sd = 1){
	# an alternative replacement to rpgm::rpgm.rnorm
		# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	(RcppZiggurat::zrnorm(n) * sd) + mean
	}

