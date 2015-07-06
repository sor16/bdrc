Adistclean <- function(W){

    W_unique=unique(W)
    N=length(W)
    n=length(W_unique)
    A=matrix(0,nrow=N,ncol=n)
    #Initialize A and add first value
    A[1,1]=1
    e=1
    for(i in 2:N){
        if( W[i]==W[i-1]){

            A[i,e]=1

        }
        else{
            e=e+1
            A[i,e]=1
        }
    }

    dist=abs(outer(W_unique,W_unique,FUN="-"))
    return(list("dist"=dist,"A"=A,"O","N"=N,"n"=n))
}
