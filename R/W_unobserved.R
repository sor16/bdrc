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
W_unobserved <- function(RC,w_min=NA,w_max=NA){
    w_u=NULL
    w=100*c(RC$O) #work in cm
    max_w_diff=5
    #distance between subsequent elements in vector with additional dummy point 1000
    distvect=abs(w-c(w[2:length(w)],1000))
    #add datapoints to corresponding distances to see range of distance
    distwithdata=rbind(w,distvect,c(w[2:length(w)],1000))
    distfilter=distwithdata[,distvect>max_w_diff]
    #remove dummy distance
    distfilter=as.matrix(distfilter[,-ncol(distfilter)])
    if(ncol(distfilter)!=0){
        #make sequence from the ranges with length.out equal to corresponding elelement in distvect
        w_u=0.01*unlist(apply(distfilter,2,FUN=function(x){setdiff(seq(x[1],x[3],length.out=2+ceiling(x[2]/max_w_diff)),c(x[1],x[3]))
        }))
    }
    w_before_data=setdiff(seq(w_min,RC$w_min,by=0.05),c(RC$w_min))
    w_after_data=setdiff(seq(RC$w_max,w_max,length.out=2+ceiling(20*(w_max-RC$w_max))),RC$w_max)
    w_u=c(w_before_data,w_u,w_after_data)
    return(w_u)
}
