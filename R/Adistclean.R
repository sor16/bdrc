Adistclean <- function(W){

    W_unique=unique(W)
    N=length(W)
    n=length(W_unique)
    A=matrix(0,nrow=N,ncol=n)
    #Initialize A and add first value
    A[1,1]=1
    col=1
    for(row in 2:N){
        if( W[row]==W[row-1]){

            A[row,col]=1

        }
        else{
            col=col+1
            A[row,col]=1
        }
    }

    dist=abs(outer(W_unique,W_unique,FUN="-"))
    return(list("dist"=dist,"A"=A,"O","N"=N,"n"=n))
}
