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
XX <- spline_functions(ZZ,tau_x,dx,kx,M)

  return(XX)
}
