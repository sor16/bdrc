B_splines <- function(ZZ){

#A script to test the B-splines.

#The number of equally spaced interior knots.
kx=2

#Delta x and Delta y.
dx=1/(kx+1)

#The order of the splines.
M = 4

#Determine the number of functions.
Nx = kx + M
#
#The epsilon-knots
# epsilon_x = dx*[0:(kx+1)];
epsilon_x = dx*seq(0,kx+1,by=1)


#the tau-knots.
# tau_x = zeros(1,kx+2*M);
# tau_x(1:M) = epsilon_x(1)*ones(1,M);
# tau_x(M+1:kx+M) = epsilon_x(2:kx+1);
# tau_x(kx+M+1:kx+2*M) = epsilon_x(kx+2)*ones(1,M);
tau_x = matrix(0,nrow=1,ncol=(kx+2*M))
tau_x[1:M] = epsilon_x[1]*matrix(1,nrow=1,ncol=M)
tau_x[(M+1):(kx+M)]=epsilon_x[2:(kx+1)]
tau_x[(kx+M+1):(kx+2*M)]=epsilon_x[kx+2]*matrix(1,nrow=1,ncol=M)

#Vector with values of x and y.
lx = length(ZZ)

#Compute the x-splines and the y-splines.
#[XX] = spline_functions(ZZ,tau_x,dx,kx,M);

    # function [XX] = spline_functions(ZZ,tau_x,dx,kx,M)
    #
    # XX = zeros(kx+M,length(ZZ));
    XX = matrix(0,nrow=(kx+M),ncol=length(ZZ))
    #
    # % i = 1
    XX[1,] = (1/dx^3)*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1]);

    # % i = 2
    XX[2,] = (1/dx^3)*(ZZ-tau_x[2])*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/2/dx^3)*(tau_x[M+2]-ZZ)*(ZZ-tau_x[3])*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/4/dx^3)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/4/dx^3)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])

    # % i = 3
    XX[3,] = (1/2/dx^3)*(ZZ-tau_x[3])*(ZZ-tau_x[3])*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/4/dx^3)*(ZZ-tau_x[3])*(tau_x[M+2]-ZZ)*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/4/dx^3)*(ZZ-tau_x[3])*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
        (1/6/dx^3)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M])*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
        (1/6/dx^3)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M])*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
        (1/6/dx^3)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M+1])*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
        (1/6/dx^3)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(tau_x[M+2]<=ZZ)*(ZZ<tau_x[M+3])

    # % i = 4,...,kx + 1
    #   for (kk in M:(kx + 1)){
    #     XX[kk,] = (1/6/dx^3)*(ZZ-tau_x[kk])*(ZZ-tau_x[kk])*(ZZ-tau_x[kk])*(tau_x[kk]<=ZZ)*(ZZ<tau_x[kk+1])+
    #       (1/6/dx^3)*(ZZ-tau_x[kk])*(ZZ-tau_x[kk])*(tau_x[kk+2]-ZZ)*(tau_x[kk+1]<=ZZ)*(ZZ<tau_x[kk+2])+
    #       (1/6/dx^3)*(ZZ-tau_x[kk])*(tau_x[kk+3]-ZZ)*(ZZ-tau_x[kk+1])*(tau_x[kk+1]<=ZZ)*(ZZ<tau_x[kk+2])+
    #       (1/6/dx^3)*(ZZ-tau_x[kk])*(tau_x[kk+3]-ZZ)*(tau_x[kk+3]-ZZ)*(tau_x[kk+2]<=ZZ)*(ZZ<tau_x[kk+3])+
    #       (1/6/dx^3)*(tau_x[kk+4]-ZZ)*(ZZ-tau_x[kk+1])*(ZZ-tau_x[kk+1])*(tau_x[kk+1]<=ZZ)*(ZZ<tau_x[kk+2])+
    #       (1/6/dx^3)*(tau_x[kk+4]-ZZ)*(ZZ-tau_x[kk+1])*(tau_x[kk+3]-ZZ)*(tau_x[kk+2]<=ZZ)*(ZZ<tau_x[kk+3])+
    #       (1/6/dx^3)*(tau_x[kk+4]-ZZ)*(tau_x[kk+4]-ZZ)*(ZZ-tau_x[kk+2])*(tau_x[kk+2]<=ZZ)*(ZZ<tau_x[kk+3])+
    #       (1/6/dx^3)*(tau_x[kk+4]-ZZ)*(tau_x[kk+4]-ZZ)*(tau_x[kk+4]-ZZ)*(tau_x[kk+3]<=ZZ)*(ZZ<tau_x[kk+4])
    #   }
    #
    #Axel/check-spline-algorithm, in Matlab this for loop does not run. Change to 4:-1:3 in Matlab?
    # % i = kx + 2
    XX[kx+2,] =  -(1/6/dx^3)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]<=ZZ)*(ZZ<tau_x[kx+3])-
        (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+4])*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
        (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
        (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+5])*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
        (1/4/dx^3)*(ZZ-tau_x[kx+6])*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
        (1/4/dx^3)*(ZZ-tau_x[kx+6])*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
        (1/2/dx^3)*(ZZ-tau_x[kx+6])*(ZZ-tau_x[kx+6])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])

    # % i = kx + 3
    XX[kx+3,] = - (1/4/dx^3)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
        (1/4/dx^3)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
        (1/2/dx^3)*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+6])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
        (1/dx^3)*(ZZ-tau_x[kx+7])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])

    # % i = kx + 4
    XX[kx+4,] = -(1/dx^3)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<=tau_x[kx+5])

    XX = t(XX)

  return(XX)
}
