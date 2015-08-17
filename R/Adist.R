#'Distance matrix for unique water level measurements
#'
#'Adist extracts unique elements of water level measurements and creates a distance matrix
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

library(Matrix)



w=as.matrix(w)
N=length(w[,1])
O=unique(w[,1])

n_e=length(O)


A=matrix(0,nrow=N,ncol=n_e)

A[1,1]=1

e=1
for(ee in 2:N){
  if(w[ee]==w[ee-1]){

A[ee,e]=1

  }else{
    e=e+1
A[ee,e]=1
#O[ee]=w[ee]
  }
  }



O=t(O)
#A=diag(1,nrow=N,ncol=N)
w2=O

n=length(O)

W=O

for(ee in 2:n){
  W=rbind(W,O)
}
  dist=abs(W-t(W))
return(list("dist"=dist,"A"=A,"n"=n,"N"=N,"O"=O))
}
