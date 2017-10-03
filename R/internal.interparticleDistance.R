#   Calculate interparticle distance
#   
#   This function calculates interparticle distance between two vectors
#   
#   IPD is used by several internal TreEvo function to compare two
#   distributions.
#   
#   @param x Vector of parameter values
#   @param y Vector of parameter values
#   @param abs If TRUE, uses absolute value for distance
#   @return Returns a matrix of distances between pairwise points in the vectors
#   @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished

interparticleDistance<-function(x,y,abs=TRUE) {
   distances<-matrix(nrow=length(x),ncol=length(y))
   for (i in sequence(length(x))) {
      for (j in sequence(length(y))) {
        distances[i,j]<-x[i]-y[j] 
      }
   }
   if (abs==TRUE) {
      distances<-abs(distances) 
   }
   return(distances)
}
