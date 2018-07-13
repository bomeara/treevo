test_that("highestDensityInterval works", {
set.seed(444)

# let's imagine we have some variable with
	# an extreme bimodal distribution
z <- sample(c(rnorm(50, 1, 2), rnorm(100, 50, 3)))
hist(z)

# now let's say we want to know the what sort of values
# are reasonably consistent with this distribution

# for example, let's say we wanted the ranges within
# which 80% of our data occurs

# one way to do this would be a quantile
	# two tailed 80% quantiles
quantile(z, probs = c(0.1, 0.9))

# that seems overly broad - there's essentially no density
# in the central valley - but we want to exclude values found there!
# A value of 30 doesn't match this data sample, right??

# the problem is that all quantile methods are essentially based on
# the empirical cumulative distribution function - which is monotonic
# (as any cumulutative function should be), and thus
# quantiles can only be a single interval

# A different approach is to use density from stats
density(z)
# we could then take the density estimates for
# particular intervals, rank-order them, and
# then cumulative sample until we reach
# our desired probability density (alpha)

# let's try that
alpha <- 0.8
zDensOut <- density(z)
zDensity <- zDensOut$y/sum(zDensOut$y)
inHPD <- cumsum(-sort(-zDensity)) <= alpha
# now reorder
inHPD <- inHPD[order(order(-zDensity))]
colDens <- rep(1, length(zDensity))
colDens[inHPD] <- 2
# and let's plot it, with colors
plot(zDensOut$x, zDensity, col = colDens)

# and we can do all that (except the plotting)
	# with highestDensityInterval
highestDensityInterval(z, alpha = 0.8)

})

test_that("testMultivarOutlierHDR works", {
# simulate two correlated variables
set.seed(444)
x <- rnorm(100, 1, 1)
y <- (x*1.5)+rnorm(100)

# find the highest density intervals for each variable
pIntX <- highestDensityInterval(x, alpha = 0.8)
pIntY <- highestDensityInterval(y, alpha = 0.8)

# These define a box-like region that poorly
# describes the actual distribution of
# the data in multivariate space.

# Let's see this ourselves...
xLims <- c(min(c(x, pIntX)), max(c(x, pIntX)))
yLims <- c(min(c(y, pIntY)), max(c(y, pIntY)))
plot(x, y, xlim = xLims, ylim = yLims)
rect(pIntX[1], pIntY[1], pIntX[2], pIntY[2])

# So, that doesn't look good.
# Let's imagine we wanted to test if some outlier
# was within that box:

outlier <- c(2, -1)
points(x = outlier[1], y = outlier[2], col = 2)

# we can use testMultivarOutlierHDR with pca = FALSE
# to do all of the above without visually checking
testMultivarOutlierHDR(dataMatrix = cbind(x, y), 
	outlier = outlier, alpha = 0.8, pca = FALSE)
	
# Should that really be considered to be within
# the 80% density region of this data cloud?

#####

# let's try it with PCA

pcaRes <- princomp(cbind(x, y))
PC <- pcaRes$scores

pIntPC1 <- highestDensityInterval(PC[, 1], alpha = 0.8)
pIntPC2 <- highestDensityInterval(PC[, 2], alpha = 0.8)

# plot it
xLims <- c(min(c(PC[, 1], pIntPC1)), max(c(PC[, 1], pIntPC1)))
yLims <- c(min(c(PC[, 2], pIntPC2)), max(c(PC[, 2], pIntPC2)))
plot(PC[, 1], PC[, 2], xlim = xLims, ylim = yLims)
rect(pIntPC1[1], pIntPC2[1], pIntPC1[2], pIntPC2[2])

# That looks a lot better, isnt' filled with lots of
# white space not supported by the observed data.

# And now let's apply testMultivarOutlierHDR, with pca = TRUE
testMultivarOutlierHDR(dataMatrix = cbind(x, y), 
	outlier = outlier, alpha = 0.8, pca = TRUE)

#####################

# Example with four complex variables
	# including correlated and multimodal variables

x <- rnorm(100, 1, 1)
y <- (x*0.8)+rnorm(100)
z <- sample(c(rnorm(50, 3, 2), 
  rnorm(50, 30, 3)))
a <- sample(c(rnorm(50, 3, 2), 
  rnorm(50, 10, 3)))+x^2

#plot(x, y)
#plot(x, z)
#plot(x, a)
data <- cbind(x, y, z, a)

# actual outlier, but maybe not obvious if PCA isn't applied
outlier <- c(2, 0.6, 10, 8)
# this one should appear to be an outlier (but only if PCA is applied)
testMultivarOutlierHDR(dataMatrix = data,
  outlier = outlier, alpha = 0.8)
testMultivarOutlierHDR(dataMatrix = data, 
  outlier = outlier, alpha = 0.8, pca = FALSE)

# this one should be within the 80% area
outlier <- c(1, 0, 30, 5)
testMultivarOutlierHDR(dataMatrix = data, 
  outlier = outlier, alpha = 0.8)
testMultivarOutlierHDR(dataMatrix = data, 
  outlier = outlier, alpha = 0.8, pca = FALSE)

# this one should be an obvious outlier no matter what
outlier <- c(3, -2, 20, 18)
# this one should be outside the 80% area
testMultivarOutlierHDR(dataMatrix = data, 
  outlier = outlier, alpha = 0.8)
testMultivarOutlierHDR(dataMatrix = data, 
  outlier = outlier, alpha = 0.8, pca = FALSE)
})

test_that("summarizePosteriors works",{
# example with output from doRun_prc
data(simRunExample)
#'
# alpha = 0.95
summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.95)

# you might be tempted to use alphas like 95%, but with bayesian statistics
# we often don't sample the distribution well enough to know
# its shape to exceeding detail. alpha = 0.8 may be more reasonable.
summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.8)

# or even better, for coverage purposes, maybe 0.5
summarizePosterior(results[[1]]$particleDataFrame, alpha = 0.5)
})