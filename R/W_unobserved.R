W_unobserved <- function(W_unique,min=NULL,max=NULL){
    W_u=NULL
    W_u_tild=NULL
    #distance between subsequent elements in vector with additional dummy point 1000
    min=round(min,2)
    max=round(max,2)
    w=100*W_unique
    distvect=abs(w-c(w[2:length(w)],1000))
    #add datapoints to corresponding distances to see range of distance
    distwithdata=rbind(w,distvect,c(w[2:length(w)],1000))
    distfilter=distwithdata[,distvect>4]

    if(!is.null(ncol(distfilter))){
        #remove dummy distance
        distfilter=distfilter[,-ncol(distfilter)]
        #make sequence from the ranges with length.out equal to corresponding elelement in distvect
        W_u=0.01*unlist(apply(distfilter,2,FUN=function(x){setdiff(seq(x[1],x[3],length.out=round(x[2])),c(x[1],x[3]))
        }))
    }
    if(!is.null(min)|!is.null(max)){
        minseq=setdiff(seq(min,min(W_unique),by=0.01),c(min(W_unique)))
        maxseq=setdiff(seq(round(max(W_unique),2),max,by=0.01),c(round(max(W_unique),2)))
        W_spline=c(rep(min(W_unique),length(minseq)),W_u,rep(max(W_unique),length(maxseq)))
        W_u=c(minseq,W_u,maxseq)
        W_u_tild=W_spline-min(W_unique)
    }
    else{
        W_u_tild=W_u-min(W_unique)
    }
    return(list("W_u"=W_u,"W_u_tild"=W_u_tild))
}
