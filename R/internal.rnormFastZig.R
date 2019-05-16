# an alternative replacement to rpgm::rpgm.rnorm
	# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	# 04-25-19

rnormFastZig <- function(nZig, meanZig = 0, sdZig = 1){
	# an alternative replacement to rpgm::rpgm.rnorm
		# uses a C++ implementation of the Zigguratt algorithm
		# RcppZiggurat::znorm
	(RcppZiggurat::zrnorm(nZig) * sdZig) + meanZig
	}

