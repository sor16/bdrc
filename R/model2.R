#Model2BH
library(stats)
model2<-function(wq,RC){
    Nit=20000;

    RC$mu_a=3.20
    RC$mu_b=2.29
    RC$sig_a=sqrt(1.21)
    RC$sig_b=sqrt(0.48)
    RC$p_ab=-0.61
    RC$mu_c=1.9000
    RC$nugget=10^-8
    RC$mu_sb=0.5
    RC$mu_pb=0.5
    RC$tau_pb2=0.25^2
    RC$s=3
    RC$v=5
    RC$y=rbind(as.matrix(log(wq[,2])),0)
    RC$w=as.matrix(0.01*wq[,1])
    RC$w_tild=RC$w-min(RC$w)

    Adist1 <- Adist(RC$w)
    RC$A=Adist1$A
    RC$dist=Adist1$dist
    RC$n=Adist1$n
    RC$N=Adist1$N

    RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
    RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
    RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n))) #Setja i RC

    RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
    RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))
    RC$m1=matrix(0,nrow=2,ncol=RC$n)
    RC$m2=matrix(0,nrow=RC$n,ncol=2)
    theta.init=rep(0,9)

    Dens = function(th) {-Densevalm22(th,RC)$p}

    Densmin=optim(par=theta.init,Dens,method="BFGS",hessian=TRUE)

    t_m =Densmin$par
    H=Densmin$hessian

    t1=matrix(0,9,Nit)
    t2=matrix(0,9,Nit)
    t3=matrix(0,9,Nit)
    t4=matrix(0,9,Nit)


    xsiz=max(dim(mu));
    x1=matrix(0,xsiz,Nit)
    x2=matrix(0,xsiz,Nit)
    x3=matrix(0,xsiz,Nit)
    x4=matrix(0,xsiz,Nit)


    for(j in 1:4){
      t_old=t_m
      t=matrix(0,nrow=9,ncol=Nit)
      x=matrix(0,nrow=xsiz,ncol=Nit)
      yp=matrix(0,nrow=RC$N,ncol=Nit)
      ypo=matrix(0,nrow=RC$N,ncol=Nit)
      varr=matrix(0,nrow=RC$N,ncol=Nit)
      D=matrix(0,nrow=1,ncol=Nit)



      Dens<-Densevalm22(t_old,RC)
      p_old=Dens$p
      x_old=Dens$x
      yp_old=Dens$yp
      ypo_old=Dens$ypo
      D_old=Dens$D
      varr_old=Dens$varr

      for(i in 1:Nit){
        t_new=t_old+solve(t(LH),as.matrix(rnorm(9,0,1)))

        Densnew<-Densevalm22(t_new,RC)
        p_new=Densnew$p
        x_new=Densnew$x
        yp_new=Densnew$yp
        ypo_new=Densnew$ypo
        D_new=Densnew$D
        varr_new=Densnew$varr

        logR=p_new-p_old

        if (logR>log(runif(1))){
          t_old=t_new
          x_old=x_new
          p_old=p_new
          yp_old=yp_new
          ypo_old=ypo_new
          D_old=D_new
          varr_old=varr_new
        }

        t[,i]=t_old
        yp[,i]=yp_old
        ypo[,i]=ypo_old

        D[1,i]=D_old
        varr[,i]=varr_old
      }

      if(j==1){
        t1=t
        yp1=yp
        ypo1=ypo
        D1=D
        varr1=varr
      } else if(j==2){
        t2=t
        yp2=yp
        ypo2=ypo
        D2=D
        varr2=varr
      } else if(j==3){
        t3=t
        yp3=yp
        ypo3=ypo
        D3=D
        varr3=varr
      } else if(j==4){
        t4=t
        yp4=yp
        ypo4=ypo
        D4=D
        varr4=varr
      }
    }
}
