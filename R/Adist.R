#equivalent to
#unique(wq[,1])


Adist <- function(w){

library(Matrix)
  # %Zdist creates the
  #
  # %Input w: water level measurements
  #
  # %Outputs
  # %Z:     Matrix \mathbold{Z} linking unique water level measurements (mathbold{w}') to actual
  # %       water level measurements (w) such that \mathbold{w}=mathbold{Zw}'
  # %dist:  Matrix of distances between unique water level measurements
  # %       dist_{ij}=|w_{i}'-w_{j}'|
  # %n:     Number of unique measurements
  # %N:     Number of measurements


w=as.matrix(w)
N=length(w[,1])
O=unique(w[,1])

n_e=length(O)


A=matrix(0,nrow=N,ncol=n_e)

A[1,1]=1

e=1
for(ee in 2:N){
  if( w[ee]==w[ee-1]){

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
return(list("dist"=dist,"A"=A,"n"=n,"N"=N))
}
