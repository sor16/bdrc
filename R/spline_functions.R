spline_functions <- function(ZZZ,tau,dx,k,M){
  # function [XX] = spline_functions(ZZZ,tau,dx,k,M)
  # 
  # XX = zeros(k+M,length(ZZZ));
  XX = matrix(0,nrow=(k+M),ncol=length(ZZZ))
  # 
  # % i = 1
  XX[1,] = (1/dx^3)*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1]);
  
  # % i = 2
  XX[2,] = (1/dx^3)*(ZZZ-tau[2])*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
    (1/2/dx^3)*(tau[M+2]-ZZZ)*(ZZZ-tau[3])*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
    (1/4/dx^3)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
    (1/4/dx^3)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])
  
  # % i = 3
  XX[3,] = (1/2/dx^3)*(ZZZ-tau[3])*(ZZZ-tau[3])*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
    (1/4/dx^3)*(ZZZ-tau[3])*(tau[M+2]-ZZZ)*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+ 
    (1/4/dx^3)*(ZZZ-tau[3])*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
    (1/6/dx^3)*(tau[M+3]-ZZZ)*(ZZZ-tau[M])*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
    (1/6/dx^3)*(tau[M+3]-ZZZ)*(ZZZ-tau[M])*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
    (1/6/dx^3)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(ZZZ-tau[M+1])*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
    (1/6/dx^3)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(tau[M+2]<=ZZZ)*(ZZZ<tau[M+3])
  
  # % i = 4,...,k + 1
#   for (kk in M:(k + 1)){
#     XX[kk,] = (1/6/dx^3)*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(tau[kk]<=ZZZ)*(ZZZ<tau[kk+1])+
#       (1/6/dx^3)*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(tau[kk+2]-ZZZ)*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+
#       (1/6/dx^3)*(ZZZ-tau[kk])*(tau[kk+3]-ZZZ)*(ZZZ-tau[kk+1])*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+ 
#       (1/6/dx^3)*(ZZZ-tau[kk])*(tau[kk+3]-ZZZ)*(tau[kk+3]-ZZZ)*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
#       (1/6/dx^3)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+1])*(ZZZ-tau[kk+1])*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+ 
#       (1/6/dx^3)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+1])*(tau[kk+3]-ZZZ)*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
#       (1/6/dx^3)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+2])*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
#       (1/6/dx^3)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(tau[kk+3]<=ZZZ)*(ZZZ<tau[kk+4]) 
#   }
#   
#Axel/check-spline-algorithm, in Matlab this for loop does not run. Change to 4:-1:3 in Matlab?
  # % i = k + 2
  XX[k+2,] =  -(1/6/dx^3)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(tau[k+2]<=ZZZ)*(ZZZ<tau[k+3])-
    (1/6/dx^3)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+4])*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
    (1/6/dx^3)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
    (1/6/dx^3)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+5])*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-       
    (1/4/dx^3)*(ZZZ-tau[k+6])*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
    (1/4/dx^3)*(ZZZ-tau[k+6])*(tau[k+3]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
    (1/2/dx^3)*(ZZZ-tau[k+6])*(ZZZ-tau[k+6])*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])
  
  # % i = k + 3
  XX[k+3,] = - (1/4/dx^3)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
    (1/4/dx^3)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
    (1/2/dx^3)*(tau[k+3]-ZZZ)*(ZZZ-tau[k+6])*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
    (1/dx^3)*(ZZZ-tau[k+7])*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])
  
  # % i = k + 4
  XX[k+4,] = -(1/dx^3)*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<=tau[k+5])
  
  XX = t(XX)
  
}