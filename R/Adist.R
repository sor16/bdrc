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
Adist <- function(w){
N=length(w)
O=t(unique(w))
n=length(O)

A=matrix(0,nrow=N,ncol=n)

A[1,1]=1

e=1
for(ee in 2:N){

  if(w[ee]==w[ee-1]){

        A[ee,e]=1

  }else{
    e=e+1
    A[ee,e]=1

  }
}
W=O
 for(ee in 2:n){
   W=rbind(W,O)
 }

dist=abs(W-t(W))
return(list("dist"=dist,"A"=A,"n"=n,"N"=N,"O"=O))
}
