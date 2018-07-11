#' 3D ABCplots
#' 
#' Plot posterior density distribution for each generation in 3d plot window
#' 
#' This opens a new interactive 3d plotting window and plots the posterior
#' density distribution of accepted particles from each generation.  Several
#' options are available to add to the plot: plotting particles by weight or
#' distance, plotting particle parantage, and plotting the real parameter
#' values (if known).
#' 
#' As of version 0.6.0, rejected particles are not saved for outputting by the parallelized algorithm,
#' and thus they are no longer displayed by this function, unlike previous versions.


#' @note
#' This function requires access functions \code{triangulate} and the \code{as} method for
#' class \code{gpc.poly} from package \code{gpclib}. As of 05-30-18, this package
#' was not available from CRAN as a Windows binary, and thus this function is likely
#' unavailable to many (if not all) Windows users.
#' 
#' This function also requires the package \code{rgl}, which is usually easier to
#' obtain than package \code{gpclib}, but may not be buildable on some UNIX workstations.
#' 

#* @param particleDataFrame A \code{particleDataFrame} object, as found among the output from \code{\link{doRun}} functions.

#' @param parameter Column number of parameter of interest from
#' \code{particleDataFrame}.

#' @param show.particles Option to show particles on 3d plot as "none" or as a
#' function of "weights" or "distance".

#' @param plot.parent Option to plot lines on the floor of the 3d plot to show
#' particle parantage.

#' @param realParam Option to display real parameter value as a solid line,
#' also must give actual value for this (realParamValues).  Note: this should
#' only be done with simulated data where real param values are recorded.

#' @param realParamValues Value for \code{realParam}.

#' @author Barb Banbury

# @references O'Meara and Banbury, unpublished

# @keywords plotABC_3D

#' @examples
#' 
#' # need to check for required suggested packages
#' if(requireNamespace("gpclib", quietly = TRUE) & requireNamespace("rgl", quietly = TRUE)){
#' 
#'  data(simRunExample)
#'  plotABC_3D(particleDataFrame = results[[1]]$particleDataFrame, parameter = 7,
#'      show.particles = "none", plot.parent = FALSE, realParam = FALSE, realParamValues = NA)
#' 
#' }

#' @name plotABC_3D
#' @rdname plotABC_3D
#' @export
plotABC_3D <- function(particleDataFrame, parameter, show.particles = "none",
	plot.parent = FALSE, realParam = FALSE, realParamValues = NA) {

	# check if gpclib exists, if not - FAIL
	has_gpclib <- requireNamespace("gpclib", quietly = TRUE)
	if(!has_gpclib){stop(
		"This function cannot be run without package gpclib available (Note: Windows binaries of gpclib were not unavailable as of 09-12-17).")
		}
	
	# check if rgl exists, if not - FAIL
	has_rgl <- requireNamespace("rgl", quietly = TRUE)
	if(!has_rgl){stop(
		"This function cannot be run without package rgl available.")
		}
	
	#ugh
	#library(gpclib)
	#library(rgl)
	#library(geometry)
	
	# ugh ugh
	#generation <- NULL #to appease R CMD CHECK
	# yes??? I think this is right, not sure
	generation <- particleDataFrame$generation
		
	x <- particleDataFrame	
	param.position <- parameter
	pdfMaxGen <- particleDataFrame[particleDataFrame$weight>0  & particleDataFrame$generation == max(particleDataFrame$generation),]
	nParticles <- dim(pdfMaxGen)[1]
	nparams <- dim(x)[2]-6
	#
	q <- vector() #vector of x vals
	r <- vector() #vector of y vals
	s <- vector() #generation each x-y coord is found
	#
	rangeParam <- (max(particleDataFrame[particleDataFrame$weight>0,param.position])
		-min(particleDataFrame[particleDataFrame$weight>0,param.position]))
	#
	if (rangeParam  !=  0) {
			v <- vector("list", max(particleDataFrame$generation))
			for (i in 1:max(particleDataFrame$generation)){
				which.gen <- (i+1)-1 # the hell is this?
				pdfParam <- particleDataFrame[(particleDataFrame$weight>0 & particleDataFrame$generation == i),param.position]
				pdfSix <- particleDataFrame[particleDataFrame$weight>0 & particleDataFrame$generation == i,6]
				v[[i]] <- density(pdfParam,weights = nParticles*pdfSix/sum(nParticles*pdfSix))
				q <- c(q, v[[i]]$x)
				r <- c(r, v[[i]]$y)
				#s <- c(s, v[[i]]$x)	# return a$generation which v[[i]]
				#
				s <- c(s,rep(which.gen,times = length(v[[i]]$x)))
				T <- as.matrix(cbind(q, r, s))
				}
					
			x <- T[,1]
			y <- T[,2]
			z <- T[,3]
			
			#
			# plotting with rgl
			rgl::open3d()  #make bigger window
			#bg3d("color)  #gives background color for plot
			rgl::plot3d(x, y, z, col = "black", box = FALSE, type = "n", xlab = "", ylab = "", zlab = "",
				zlim = c(0, max(particleDataFrame$generation)), ylim = c(0, max(y)))
			#message(paste("HERE"))
			rgl::rgl.viewpoint(35, 1, 90)  #sets viewpoint for initial plot
			rgl::title3d(colnames(x)[param.position], col = 'red', pos = c(NA, -2, max(z)))
			#text3d(x = min(x), y = mean(y), z = max(z), text = "Density" col = 'blue')
			#title3d("Starting States", col = 'purple', pos = c(NA, 0, max(which.gen)))
			#
			
			for (i in 1:max(s)){
				ngen <- (i+1)-1
				triangles <- gpclib::triangulate(as(cbind(x[which(z == i)],y[which(z == i)]), "gpc.poly"))
				xyzGen <- data.frame(zG = z[which(z == i)], xG = x[which(z == i)], yG = y[which(z == i)])
				zfit <- predict(lm(zG ~ xG + yG, xyzGen), newdata = data.frame(xG = triangles[,1], yG = triangles[,2]))
				opacity <- 0.8*(ngen/length(v))
				rgl::rgl.material(color = "black", alpha = opacity, lit = FALSE)
				#message(dim(triangles))
				#message(dim(zfit))
				if(length(zfit) != (dim(triangles)[1])){
					stop("Error, predict() not properly returning vector for same number of input")
				}else{
					rgl::triangles3d(cbind(triangles, zfit), col = "red")
					}
				#readline(prompt = "hit enter ")
				if (realParam) {
					rgl::rgl.material(color = "blue", lwd = 2)
					rgl::lines3d(x = c(realParamValues[1], realParamValues[1]), y = c(0, 0), z = c(min(s), max(s)))	
				}
			}
		}		
		
		else {
			warning(paste("You are attempting to plot",colnames(x)[param.position],", which is a fixed value")) #return fixed val
			}
		
			
		show.particles <- match.arg(arg = show.particles, choices = c("none", "weights", "distance"),several.ok = FALSE)
		#kept <- subset(particleDataFrame[which(particleDataFrame$id>0),])[,]	
		#reject <- subset(particleDataFrame[which(particleDataFrame$id<0),])[,]
		#short.kept <- subset(kept[which(kept$generation>1),])[,]
		#
		kept <- particleDataFrame[particleDataFrame$id>0,]
		short.kept <- kept[kept$generation>1,]

		if (show.particles == "none"){
			message("Note: Currently not plotting particles.  To plot particles modify the show.particles =  argument.")
			}
		else if (show.particles == "weights") {
			for (j in 1:(dim(kept)[1])) {
				circle.size <- (kept[j, 6]/max(kept[,6]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))*10
				rgl::spheres3d(x = kept[j, param.position], y = -1*(.05*(max(y)-min(y))), z = kept[j, 1], radius = circle.size, type = "s", col = "black", add = TRUE, aspect = TRUE)
			}	
		}		
		else if (show.particles == "distance") {
			for (j in 1:(dim(kept)[1])) {
				circle.size <- (kept[j, 5]/max(kept[,5]))*(.1*(max(kept[,param.position])-min(kept[,param.position])))*15
				rgl::spheres3d(x = kept[j, param.position], y = -1*(.05*(max(y)-min(y))), z = kept[j, 1], radius = circle.size, type = "s", col = "black", add = TRUE, aspect = TRUE)
			}	
		}	
		if (plot.parent){
			for (k in 1:(dim(short.kept)[1])) {
			prev.gen <- subset(kept[which(kept$generation == short.kept[k, 1]-1),])[,]  #works to retreive prev gen
			rgl::rgl.material(color = "black", alpha = 0.5)
			rgl::lines3d(x = c(short.kept[k, param.position], prev.gen[short.kept[k,]$parentid, param.position]), y = c(0, 0), z = c(short.kept[k, 1], short.kept[k, 1]-1))		}	
		}	
	}
