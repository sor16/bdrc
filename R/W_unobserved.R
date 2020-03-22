#'Unobserved stages
#'
#'W_unobserved returns the stages that are needed to make an equally spaced grid of stages from data of stages.
#'
#'@param W_unique vector containing unique stages from river data.
#'@param min minimum stage of rating curve.
#'@param max maximum stage of rating curve.
#'@return W_unobserved returns a list of vectors, W_u and W_u_tild. W_u is a vector of unobserved stage values
#' needed to make an equally spaced grid of stages. W_u_tild is a vector which is calculated by W_u-min(W_unique) needed to input into B_splines.
#' The unobserved stages are lower or higher than that of the data, take the same value in W_u_tild as the minimum value and maximum value of the
#' data respectively. This is done to ensure constant variance below and above observed data.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
W_unobserved <- function(W_unique,min=NA,max=NA){
    W_u=NULL
    W_u_tild=NULL
    w=100*W_unique #work in cm
    filter=5
    #distance between subsequent elements in vector with additional dummy point 1000
    distvect=abs(w-c(w[2:length(w)],1000))
    #add datapoints to corresponding distances to see range of distance
    distwithdata=rbind(w,distvect,c(w[2:length(w)],1000))
    distfilter=distwithdata[,distvect>filter]
    #remove dummy distance
    distfilter=as.matrix(distfilter[,-ncol(distfilter)])
    if(ncol(distfilter)!=0){
        #make sequence from the ranges with length.out equal to corresponding elelement in distvect
        W_u=0.01*unlist(apply(distfilter,2,FUN=function(x){setdiff(seq(x[1],x[3],length.out=2+round(x[2]/filter)),c(x[1],x[3]))
        }))
    }
    if(!is.na(min)|!is.na(max)){
        min=ceiling(min*10)/10
        max=ceiling(max*10)/10
        minseq=setdiff(seq(min,min(W_unique),by=0.05),c(min(W_unique)))
        maxseq=setdiff(seq(max(W_unique),max,length.out=2+ceiling(20*(max-max(W_unique)))),max(W_unique))
        W_spline=c(rep(min(W_unique),length(minseq)),W_u,rep(max(W_unique),length(maxseq)))
        W_u=c(minseq,W_u,maxseq)
        W_u_tild=W_spline-min(W_unique)
    }else{
        W_u_tild=W_u-min(W_unique)
    }
    return(list("W_u"=W_u,"W_u_tild"=W_u_tild))
}
