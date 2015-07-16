Adistclean <- function(W){

    #W_unique=unique(W)
    N=length(W)
    distvect=abs(c(1000,W[-length(W)])-W)
    O=W[which(distvect>0.03)]
    n=length(O)
    A=matrix(0,nrow=N,ncol=n)
    #Initialize A and add first value
    A[1,1]=1
    col=1
    for(row in 2:N){
        if(abs(W[row]-W[row-1]<0.01)){

            A[row,col]=1

        }else{
            col=col+1
            A[row,col]=1
        }
    }

    dist=abs(outer(O,O,FUN="-"))
    return(list("dist"=dist,"A"=A,"O"=O,"N"=N,"n"=n))
}
