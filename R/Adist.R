#'Linking unique water level measurements to actual
#'water level measurements
#'
#'Adist links unique water level measurements (\strong{w'}) to actual
#'water level measurements (w) such that \strong{w}=\strong{Aw'}.
#'from the measurements.
#'@param w Numerical vector containing water level measurements
#'@return
#'\itemize{
#'\item A: Matrix \strong{A} linking unique water level measurements (\strong{w'}) to actual
#'water level measurements (w) such that \strong{w}=\strong{Aw'}
#'\item dist:  Matrix of distances between unique water level measurements
#'\item n:     Number of unique measurements
#'\item N:     Number of measurements
#'}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
create_A <- function(RC){
  A=matrix(0,nrow=RC$n,ncol=RC$n_unique)
  A[1,1]=1
  i=1
  for(ii in 2:RC$n){
    if(RC$w[ii]==RC$w[ii-1]){
          A[ii,i]=1
    }else{
      i=i+1
      A[ii,i]=1
    }
  }
  return(A)
}

