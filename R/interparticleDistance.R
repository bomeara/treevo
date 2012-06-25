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