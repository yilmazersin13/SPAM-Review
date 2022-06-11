a<-1234578
library(mgcv)
set.seed(a)
require(SimDesign)
SimFunctions()
des<-createDesign(sample_size = c(30, 60, 100, 150),
                  method = c("LL", "KS","SS"),no.cov=c(2,4))
Generate  <- function(condition, fixed_objects = NULL) {
  n       <- condition$sample_size
  method  <- condition$method
  criteria<- condition$criteria
  no.cov  <- condition$no.cov
  if(no.cov == 2){
    t<-0
    x<-matrix(c(runif(n),runif(n)),n,2)
    for (j in 1:n){
      t[j]<-(j-0.5)/n
    }
    z<-seq(-2,2,length.out = n)
    f1<-1-48*t+218*t^2-315*t^3+145*t^4
    f2<-sin(2*z)+2*exp(-16*z^2)
    beta<-c(-1,2)
    error<-rnorm(n,sd=0.5)
    y<-x%*%beta+f1+f2+error
    nc<-matrix(c(t,z),n,2)
    allf<-matrix(c(f1,f2),n,2)
    dat<-new.env()
    dat$x<-x
    dat$y<-y
    dat$nc<-nc
    dat$beta<-beta
    dat$allf<-allf
  } 
  else if(no.cov == 4){
    t2<-0
    t<-0
    x<-matrix(c(runif(n),runif(n)),n,2)
    for (j in 1:n){
      t[j] <-(j-0.5)/n
      t2[j]<-5*(j-0.5)/n
    }
    z<-seq(-2,2,length.out = n)
    z2<-seq(-2,2,length.out = n)
    f1<-1-48*t+218*t^2-315*t^3+145*t^4
    f2<-sin(2*z)+2*exp(-16*z^2)
    f3<-t2*(sin(t2))^2
    f4<-z2+2*exp(-15*z2^2)
    beta<-c(-1,2)
    error<-rnorm(n,sd=0.5)
    y<-x%*%beta+f1+f2+f3+f4+error
    nc<-matrix(c(t,t2,z,z2),n,4)
    allf<-matrix(c(f1,f2,f3,f4),n,4)
    dat<-new.env()
    dat$x<-x
    dat$y<-y
    dat$nc<-nc
    dat$beta<-beta
    dat$allf<-allf
  }
  dat
}
Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  SPAMLL<-function(x,y,nc,allf){
    library(optR)
    library(psych)
    library(pracma)
    size<-dim(nc)
    q<-size[2]
    #-------Improved AIC (AICc)----------
    aiccfunc<-function(y,yhat,H){
      p <- 2
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-1+log10(norm(y-yhat)^2/n)+(2*(p+1)/(n-tr(H)-2))
      return(score)
    }
    #---------------GCV-----------------
    gcvfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-(1/n)*(t(y-yhat)%*%(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
      return(score)
    }
    #---------------GCV-----------------
    cpfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      MSE<-mean((y-yhat)^2)
      ssq<-var(y-yhat)
      score<-MSE+(2*ssq*tr(H))
      return(score)
    }
    #----------------------------------
    K<-function(z){
      res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
      return(res)
    }
    #---------------------------------
    Smat<-function(z,bw){
      K<-function(z){
        res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
        return(res)
      }
      n<-length(z)
      za<-seq(min(z),max(z),length.out = n)
      W  <- matrix(0,n,n)
      Sr  <- matrix(0,n,n)
      for (i in 1:n){
        d  <- matrix(c(1,0),2,1)
        Zr <- matrix(c(matrix(1,n,1),za-z[i]),n,2)
        u<-(za-z[i])/bw
        Wdiag <- K(u)
        Wr<-diag(Wdiag)
        Sr[,i]<-t(d)%*%solve(t(Zr)%*%Wr%*%Zr,tol=1e-100)%*%t(Zr)%*%Wr
      }
      SLL <- Sr
      return(SLL)
    }
    if (q==2){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,2]
      
      
      #------------------------------------
      #LLR estimator for SPAM Part (Non-iterative)
      #Selection of Bandwidth-------------------------------------------------------
      h<-seq(0.005,0.3,length.out = n)
      gcv<-0
      aic<-0
      for (i in 1:n){
        SLL<-Smat(z1,bw=h[i])
        xtil <- (diag(n)-SLL)%*%x
        ytil <- (diag(n)-SLL)%*%y
        
        s.betaLL <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
        s.fLL    <- SLL%*%(y-x%*%s.betaLL)
        s.yhatLL <- x%*%s.betaLL+s.fLL
        s.H      <- xtil%*%solve(t(xtil)%*%xtil)%*%t(xtil)
        aic[i]<-aiccfunc(y,s.yhatLL,s.H)
      }
      bw<-h[which.min(aic)]
      # ESTIMATION AFTER SELECTION OF BW----------------------------------------------
      SLL<-Smat(z1,bw)
      
      xtil <- (diag(n)-SLL)%*%x
      ytil <- (diag(n)-SLL)%*%y
      
      s.betaLL <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
      s.fLL <- SLL%*%(y-x%*%s.betaLL)
      s.yhatLL <- x%*%s.betaLL+s.fLL
      
      #plot(s.yhatLL,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(y,type="l",ylim=c(min(y),max(y)),col="red")
      
      #------------------------------------------------------------------------------
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.005,0.3,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.05){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2)
          e[i]           <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          tol<-e[i]
          f1hat<-f1p[,(i+1)]
          i<-i+1
        }
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]+0.02
      bw_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-Smat(z1,bw_aic)
      SLLf1_gcv<-Smat(z1,bw_gcv)
      SLLf1_cp<-Smat(z1,bw_cp)
      #--------------------------------------
      #LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-f1p[,(i+1)]
        i<-i+1
      }
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv   <-f1p[,(i+1)]
        i<-i+1
      }
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-f1p[,(i+1)]
        i<-i+1
      }
      
      
      
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.01,0.3,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.05){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i])
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1)
          e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          tol<-e[i]
          f2hat<-f2p[,(i+1)]
          i<-i+1
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]+0.02
      bw_cp<-h[which.min(cp)]
      
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-Smat(z2,bw_aic)
      SLLf2_gcv<-Smat(z2,bw_gcv)
      SLLf2_cp<-Smat(z2,bw_cp)
      #AIC-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-f2p[,(i+1)]
        i<-i+1
      }
      #GCV-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-f2p[,(i+1)]
        i<-i+1
      }
      #cp-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-f2p[,(i+1)]
        i<-i+1
      }
      
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic),n,2)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv),n,2)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp),n,2)
    }
    ###############################################################################
    #FOR q=4-----------------------------------------------------------------------
    if (q==4){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      f3<-allf[,3]
      f4<-allf[,4]
      
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,3]
      z3<-nc[,2]
      z4<-nc[,4]
      #------------------------------------
      #LLR estimator for SPAM Part (Non-iterative)
      #Selection of Bandwidth---------------------------------------------------------
      h<-seq(0.01,1,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (i in 1:n){
        SLL<-Smat(z1,bw=h[i])
        xtil <- (diag(n)-SLL)%*%x
        ytil <- (diag(n)-SLL)%*%y
        
        s.betaLL <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
        s.fLL <- SLL%*%(y-x%*%s.betaLL)
        s.yhatLL <- x%*%s.betaLL+s.fLL
        s.H      <- xtil%*%solve(t(xtil)%*%xtil)%*%t(xtil)
        aic[i]<-aiccfunc(y,s.yhatLL,s.H)
      }
      bw<-h[which.min(aic)]
      # ESTIMATION AFTER SELECTION OF BW----------------------------------------------
      SLL<-Smat(z1,bw)
      
      xtil <- (diag(n)-SLL)%*%x
      ytil <- (diag(n)-SLL)%*%y
      
      s.betaLL <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
      s.fLL <- SLL%*%(y-x%*%s.betaLL)
      s.yhatLL <- x%*%s.betaLL+s.fLL
      
      plot(s.yhatLL,type="l",ylim=c(min(y),max(y)))
      par(new=TRUE)
      plot(y,type="l",ylim=c(min(y),max(y)),col="red")
      
      #---------------------------------------------------------------------------
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.01,1,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL<-solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2-f3-f4)
          e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          tol<-e[i]
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          f1hat<-(f1p[,(i+1)])
          i<-i+1 
        }
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw1_aic<-h[which.min(aic)]
      bw1_gcv<-h[which.min(gcv)]+0.02
      bw1_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-Smat(z1,bw1_aic)
      SLLf1_gcv<-Smat(z1,bw1_gcv)
      SLLf1_cp<-Smat(z1,bw1_cp)
      
      #--------------------------------------------------------------------------
      #LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      #GCV----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      #Cp----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.01,1,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f3-f4)
          e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          tol<-e[i]
          f2hat<-(f2p[,(i+1)])
          i<-i+1 
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
        
      }
      bw2_aic<-h[which.min(aic)]
      bw2_gcv<-h[which.min(gcv)]+0.025
      bw2_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-Smat(z2,bw2_aic)
      SLLf2_gcv<-Smat(z2,bw2_gcv)
      SLLf2_cp<-Smat(z2,bw2_cp)
      
      #AICc---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-(f2p[,(i+1)])
        i<-i+1 
      }
      #GCV---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-(f2p[,(i+1)])
        i<-i+1 
      }
      #Cp---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-(f2p[,(i+1)])
        i<-i+1 
      }
      #BW SELECTION WITH f3 FUNCTION----------------------------------------------
      h<-seq(0.05,1,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z3,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
          f3p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f4)
          e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
          xt <- (diag(n)-SLL)%*%x
          Hs3 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          tol<-e[i]
          f3hat<-(f3p[,(i+1)])
          i<-i+1 
        }
        aic[j]<-aiccfunc(f3,f3hat,Hs3)
        gcv[j]<-gcvfunc(f3,f3hat,Hs3)
        cp[j]<-cpfunc(f3,f3hat,Hs3)
      }
      bw3_aic<-h[which.min(aic)]
      bw3_gcv<-h[which.min(gcv)]+0.03
      bw3_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f3")
      SLLf3_aic<-Smat(z3,bw3_aic)
      SLLf3_gcv<-Smat(z3,bw3_gcv)
      SLLf3_cp<-Smat(z3,bw3_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_aic%*%(y-x%*%betaLL_aic-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_aic<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_gcv<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_cp%*%(y-x%*%betaLL_cp-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_cp<-(f3p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      #BW SELECTION WITH f4 FUNCTION----------------------------------------------
      h<-seq(0.05,0.1,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z4,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
          f4p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f3)
          e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
          xt <- (diag(n)-SLL)%*%x
          Hs4 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
          tol<-e[i]
          f4hat<-(f4p[,(i+1)])
          i<-i+1 
        }
        aic[j]<-aiccfunc(f4,f4hat,Hs4)
        gcv[j]<-gcvfunc(f4,f4hat,Hs4)
        cp[j]<-cpfunc(f4,f4hat,Hs4)
      }
      bw4_aic<-h[which.min(aic)]
      bw4_gcv<-h[which.min(gcv)]+0.02
      bw4_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f4")
      SLLf4_aic<-Smat(z4,bw4_aic)
      SLLf4_gcv<-Smat(z4,bw4_gcv)
      SLLf4_cp<-Smat(z4,bw4_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_aic%*%(y-x%*%betaLL_aic-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_aic<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #GCV-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_gcv<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #Cp-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_cp%*%(y-x%*%betaLL_cp-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_cp<-(f4p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic+f3hat_aic+f4hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv+f3hat_gcv+f4hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp+f3hat_cp+f4hat_cp
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic,f3hat_aic,f4hat_aic),n,4)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv,f3hat_gcv,f4hat_gcv),n,4)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp,f3hat_cp,f4hat_cp),n,4)
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    res<-new.env()
    res$betahat_noniter<-s.betaLL
    res$fitted_noniter<-s.yhatLL
    res$sum_fhat_noniter<-s.fLL
    #-----------------------------------------------------------------------------
    res$beta_aic   <-betaLL_aic
    res$fhat_aic   <-FHAT_aic
    res$beta_gcv   <-betaLL_gcv
    res$fhat_gcv   <-FHAT_gcv
    res$beta_cp    <-betaLL_cp
    res$fhat_cp    <-FHAT_cp
    res$fitted_aic <-yhatBF_aic
    res$fitted_gcv <-yhatBF_gcv
    res$fitted_cp  <-yhatBF_cp
    res$Smat       <-SLL 
    res$H          <-s.H 
    
    return(res)
  }     #REVISED
  SPAMKS<-function(x,y,nc,allf){
    #x: matrix of parametric covariates
    #y: response variable
    #nc: Matrix of nonparamtric covariate
    library(optR)
    library(psych)
    library(pracma)
    size<-dim(nc)
    q<-size[2]
    #-------Improved AIC (AICc)----------
    aiccfunc<-function(y,yhat,H){
      p <- 2
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-1+log10(norm(y-yhat)^2/n)+(2*(p+1)/(n-tr(H)-2))
      return(score)
    }
    #---------------GCV-----------------
    gcvfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-(1/n)*(t(y-yhat)%*%(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
      return(score)
    }
    #---------------GCV-----------------
    cpfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      MSE<-mean((y-yhat)^2)
      ssq<-var(y-yhat)
      score<-MSE+(2*ssq*tr(H))
      return(score)
    }
    #----------------------------------
    K<-function(z){
      res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
      return(res)
    }
    #---------------------------------
    WKS<-function(z,bw){
      K<-function(z){
        res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
        return(res)
      }
      n<-length(z)
      W  <- matrix(0,n,n)
      zk<-seq(min(z)-0.1,max(z)+0.1,length.out = n)
      for (i in 1:n){
        u<-(z[i]-zk)/bw
        Ku <- K(u)
        W[,i]<-Ku/sum(Ku)
      }
      WKS <- W
      return(WKS)
    }
    if (q==2){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,2]
      
      
      #------------------------------------------------------------------------------
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.1,0.3,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-WKS(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.05){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2)
          e[i]           <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          tol<-e[i]
          f1hat<-f1p[,(i+1)]
          i<-i+1
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        yhat<-x%*%betaLL+f2p[,1]+f1hat
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]+0.02
      bw_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-WKS(z1,bw_aic)
      SLLf1_gcv<-WKS(z1,bw_gcv)
      SLLf1_cp<-WKS(z1,bw_cp)
      #--------------------------------------
      #KS estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-f1p[,(i+1)]
        i<-i+1
      }
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv<-f1p[,(i+1)]
        i<-i+1
      }
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.05){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-f1p[,(i+1)]
        i<-i+1
      }
      
      
      
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.05,0.5,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-WKS(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.05){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1)
          e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          tol<-e[i]
          f2hat<-f2p[,(i+1)]
          i<-i+1
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]+0.03
      bw_cp<-h[which.min(cp)]
      
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-WKS(z2,bw_aic)
      SLLf2_gcv<-WKS(z2,bw_gcv)
      SLLf2_cp<-WKS(z2,bw_cp)
      #AIC-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-f2p[,(i+1)]
        i<-i+1
      }
      #GCV-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-f2p[,(i+1)]
        i<-i+1
      }
      #cp-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.05){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-f2p[,(i+1)]
        i<-i+1
      }
      H <- Hs2
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic),n,2)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv),n,2)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp),n,2)
    }
    ###############################################################################
    #FOR q=4-----------------------------------------------------------------------
    if (q==4){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      f3<-allf[,3]
      f4<-allf[,4]
      
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,3]
      z3<-nc[,2]
      z4<-nc[,4]
      
      #---------------------------------------------------------------------------
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.01,0.5,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-WKS(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2-f3-f4)
          e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          tol<-e[i]
          f1hat<-(f1p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw1_aic<-h[which.min(aic)]
      bw1_gcv<-h[which.min(gcv)]
      bw1_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-WKS(z1,bw1_aic)
      SLLf1_gcv<-WKS(z1,bw1_gcv)
      SLLf1_cp<-WKS(z1,bw1_cp)
      
      #--------------------------------------------------------------------------
      #LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-(f1p[,(i+1)])
        i<-i+1 
      }
      #GCV----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv<-(f1p[,(i+1)])
        i<-i+1 
      }
      #Cp----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-(f1p[,(i+1)])
        i<-i+1 
      }
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.05,0.5,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-WKS(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f2p[,i]-f1-f3-f4)
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f3-f4)
          e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          tol<-e[i]
          f2hat<-(f2p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
        
      }
      bw2_aic<-h[which.min(aic)]
      bw2_gcv<-h[which.min(gcv)]
      bw2_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-WKS(z2,bw2_aic)
      SLLf2_gcv<-WKS(z2,bw2_gcv)
      SLLf2_cp<-WKS(z2,bw2_cp)
      
      #AICc---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f2p[,i]-f1-f3-f4)
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-(f2p[,(i+1)])
        i<-i+1 
      }
      #GCV---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f2p[,i]-f1-f3-f4)
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-(f2p[,(i+1)])
        i<-i+1 
      }
      #Cp---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f2p[,i]-f1-f3-f4)
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-(f2p[,(i+1)])
        i<-i+1 
      }
      #BW SELECTION WITH f3 FUNCTION----------------------------------------------
      h<-seq(0.05,1,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-WKS(z3,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f3p[,i]-f1-f2-f4)
          f3p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f4)
          e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
          tol<-e[i]
          f3hat<-(f3p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs3 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f3,f3hat,Hs3)
        gcv[j]<-gcvfunc(f3,f3hat,Hs3)
        cp[j]<-cpfunc(f3,f3hat,Hs3)
      }
      bw3_aic<-h[which.min(aic)]
      bw3_gcv<-h[which.min(gcv)]
      bw3_cp<-h[which.min(cp)]
      
      #plot(aic,type="l", main="AIC scores for bandwidth for f3")
      SLLf3_aic<-WKS(z3,bw3_aic)
      SLLf3_gcv<-WKS(z3,bw3_gcv)
      SLLf3_cp<-WKS(z3,bw3_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f3p[,i]-f1-f2-f4)
        f3p[,(i+1)] <- SLLf3_aic%*%(y-x%*%betaLL_aic-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_aic<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f3p[,i]-f1-f2-f4)
        f3p[,(i+1)] <- SLLf3_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_gcv<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f3p[,i]-f1-f2-f4)
        f3p[,(i+1)] <- SLLf3_cp%*%(y-x%*%betaLL_cp-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_cp<-(f3p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      #BW SELECTION WITH f4 FUNCTION----------------------------------------------
      h<-seq(0.05,0.1,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-WKS(z4,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f4p[,i]-f1-f2-f3)
          f4p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f3)
          e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
          tol<-e[i]
          f4hat<-(f4p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs4 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f4,f4hat,Hs4)
        gcv[j]<-gcvfunc(f4,f4hat,Hs4)
        cp[j]<-cpfunc(f4,f4hat,Hs4)
      }
      bw4_aic<-h[which.min(aic)]
      bw4_gcv<-h[which.min(gcv)]
      bw4_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f4")
      SLLf4_aic<-WKS(z4,bw4_aic)
      SLLf4_gcv<-WKS(z4,bw4_gcv)
      SLLf4_cp<-WKS(z4,bw4_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f4p[,i]-f1-f2-f3)
        f4p[,(i+1)] <- SLLf4_aic%*%(y-x%*%betaLL_aic-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_aic<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #GCV-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f4p[,i]-f1-f2-f3)
        f4p[,(i+1)] <- SLLf4_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_gcv<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #Cp-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f4p[,i]-f1-f2-f3)
        f4p[,(i+1)] <- SLLf4_cp%*%(y-x%*%betaLL_cp-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_cp<-(f4p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      H <- Hs4
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic+f3hat_aic+f4hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv+f3hat_gcv+f4hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp+f3hat_cp+f4hat_cp
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic,f3hat_aic,f4hat_aic),n,4)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv,f3hat_gcv,f4hat_gcv),n,4)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp,f3hat_cp,f4hat_cp),n,4)
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    res<-new.env()
    #Results obtained from backfitting algorithm
    res$beta_aic      <-betaLL_aic
    res$fhat_aic      <-FHAT_aic
    res$fhat_gcv      <-FHAT_gcv
    res$fhat_cp       <-FHAT_cp
    res$fitted_aic    <-yhatBF_aic
    res$fitted_gcv    <-yhatBF_gcv
    res$fitted_cp     <-yhatBF_cp
    res$Smat          <-SLL 
    res$H             <-H 
    return(res)
  }     #REVISED
  SPAMSS<-function(x,y,nc,allf){
    #y <- scale(y)
    #x <- scale(x)
    #nc <- scale(nc)
    #allf <- scale(allf)
    #x: matrix of parametric covariates
    #y: response variable
    #nc: Matrix of nonparamtric covariate
    library(optR)
    library(psych)
    library(pracma)
    size<-dim(nc)
    q<-size[2]
    #-------Improved AIC (AICc)----------
    aiccfunc<-function(y,yhat,H){
      p <- 2
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-1+log10(norm(y-yhat)^2/n)+(2*(p+1)/(n-tr(H)-2))
      return(score)
    }
    #---------------GCV-----------------
    gcvfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      score<-(1/n)*(t(y-yhat)%*%(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
      return(score)
    }
    #---------------GCV-----------------
    cpfunc<-function(y,yhat,H){
      y<-matrix(c(y))
      yhat<-matrix(c(yhat))
      n<-length(y)
      MSE<-mean((y-yhat)^2)
      ssq<-var(y-yhat)
      score<-MSE+(2*ssq*tr(H))
      return(score)
    }
    #----------------------------------
    K<-function(z){
      res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
      return(res)
    }
    #---------------------------------
    Smat<-function(z,bw){
      h<-0
      n<-length(z)
      for (b in 1:(n-1)) {
        h[b]<-z[b+1]-z[b]
      }
      Q<-matrix(0,n,(n-2))
      for (i in 1:(n-2)) {
        for (j in 2:n) {
          if (i==(j-1)) {
            Q[j-1,i]<-(1/h[j-1])
          }
          if (i==j) {
            Q[j-1,i]<- (-(1/h[j-1])+(1/h[j]))
          }
          if (i==(j+1)) {
            Q[j-1,i]<-(1/h[j])
          }
          if (abs(i-j)>=2) {
            Q[j-1,i]<-0
          }
        }
      }
      R<-matrix(0,n-2,n-2)
      for (i in 2:(n-1)){
        for (j in 2:(n-1)) {
          if (i==j) {
            R[j-1,i-1]<-1/3*(h[j-1]+h[j])
          }
          if (i==(j-1)) {
            R[j-1,i-1]<-1/6*h[j]
          }
          if (i==(j+1)){
            R[j-1,i-1]<-1/6*h[j]
          }
          if (abs(i-j)>=2) {
            R[j-1,i-1]<-0
          }
        }
      }
      K<-matrix(0,n,n)
      K<-(Q%*%solve(R)%*%t(Q))
      SLL<-solve(diag(n)+bw*K)
      return(SLL)
    }
    if (q==2){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,2]
      
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.000000000005,0.000000005,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2)
          e[i]           <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          tol<-e[i]
          f1hat<-f1p[,(i+1)]
          i<-i+1
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]
      bw_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-Smat(z1,bw_aic)
      SLLf1_gcv<-Smat(z1,bw_gcv)
      SLLf1_cp<-Smat(z1,bw_cp)
      #--------------------------------------
      #LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-f1p[,(i+1)]
        i<-i+1
      }
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv<-f1p[,(i+1)]
        i<-i+1
      }
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1p[,i])-f2)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2)
        e[i]        <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-f1p[,(i+1)]
        i<-i+1
      }
      
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.000000000005,0.00000005,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i])
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1)
          e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          tol<-e[i]
          f2hat<-f2p[,(i+1)]
          i<-i+1
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
      }
      bw_aic<-h[which.min(aic)]
      bw_gcv<-h[which.min(gcv)]
      bw_cp<-h[which.min(cp)]
      
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-Smat(z2,bw_aic)
      SLLf2_gcv<-Smat(z2,bw_gcv)
      SLLf2_cp<-Smat(z2,bw_cp)
      #AIC-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-f2p[,(i+1)]
        i<-i+1
      }
      #GCV-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-f2p[,(i+1)]
        i<-i+1
      }
      #cp-------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-(f1)-f2p[,i])
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1)
        e[i]        <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-f2p[,(i+1)]
        i<-i+1
      }
      H <- Hs2
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic),n,2)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv),n,2)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp),n,2)
    }
    ###############################################################################
    #FOR q=4-----------------------------------------------------------------------
    if (q==4){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      f3<-allf[,3]
      f4<-allf[,4]
      
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,3]
      z3<-nc[,2]
      z4<-nc[,4]
      
      
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      h<-seq(0.00000000005,0.00000005,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z1,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL<-solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
          f1p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f2-f3-f4)
          e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
          tol<-e[i]
          f1hat<-(f1p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f1,f1hat,Hs)
        gcv[j]<-gcvfunc(f1,f1hat,Hs)
        cp[j]<-cpfunc(f1,f1hat,Hs)
      }
      bw1_aic<-h[which.min(aic)]
      bw1_gcv<-h[which.min(gcv)]#+0.00000001
      bw1_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f1")
      SLLf1_aic<-Smat(z1,bw1_aic)
      SLLf1_gcv<-Smat(z1,bw1_gcv)
      SLLf1_cp<-Smat(z1,bw1_cp)
      
      #--------------------------------------------------------------------------
      #LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_aic%*%(y-x%*%betaLL_aic-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_aic<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      #GCV----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_gcv%*%(y-x%*%betaLL_gcv-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_gcv<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      #Cp----------------------------------------------------------------------- 
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1p[,i]-f2-f3-f4)
        f1p[,(i+1)] <- SLLf1_cp%*%(y-x%*%betaLL_cp-f2-f3-f4)
        e[i]  <- mean(abs((f1p[,(i+1)]-f1p[,i])))
        tol<-e[i]
        f1hat_cp<-(f1p[,(i+1)])
        i<-i+1 
      }
      
      ################################################################################
      #BW SELECTION WITH f2 FUNCTION----------------------------------------------
      h<-seq(0.000000000005,0.00000005,length.out = n)
      gcv<-0
      aic<-0
      cp<-0
      for (j in 1:n){
        SLL<-Smat(z2,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        i<-1
        tol<-99
        e<-0
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
          f2p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f3-f4)
          e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
          tol<-e[i]
          f2hat<-(f2p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs2 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f2,f2hat,Hs2)
        gcv[j]<-gcvfunc(f2,f2hat,Hs2)
        cp[j]<-cpfunc(f2,f2hat,Hs2)
        
      }
      bw2_aic<-h[which.min(aic)]
      bw2_gcv<-h[which.min(gcv)]#+0.00000001
      bw2_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f2")
      SLLf2_aic<-Smat(z2,bw2_aic)
      SLLf2_gcv<-Smat(z2,bw2_gcv)
      SLLf2_cp<-Smat(z2,bw2_cp)
      
      #AICc---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_aic%*%(y-x%*%betaLL_aic-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_aic<-(f2p[,(i+1)])
        i<-i+1 
      }
      #GCV---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_gcv%*%(y-x%*%betaLL_gcv-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_gcv<-(f2p[,(i+1)])
        i<-i+1 
      }
      #Cp---------------------------------------------------------------------------
      i<-1
      tol<-99
      e<-0
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2p[,i]-f3-f4)
        f2p[,(i+1)] <- SLLf2_cp%*%(y-x%*%betaLL_cp-f1-f3-f4)
        e[i]  <- mean(abs((f2p[,(i+1)]-f2p[,i])))
        tol<-e[i]
        f2hat_cp<-(f2p[,(i+1)])
        i<-i+1 
      }
      #BW SELECTION WITH f3 FUNCTION----------------------------------------------
      h<-seq(0.0000000005,0.00000005,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z3,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
          f3p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f4)
          e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
          tol<-e[i]
          f3hat<-(f3p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs3 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f3,f3hat,Hs3)
        gcv[j]<-gcvfunc(f3,f3hat,Hs3)
        cp[j]<-cpfunc(f3,f3hat,Hs3)
      }
      bw3_aic<-h[which.min(aic)]
      bw3_gcv<-h[which.min(gcv)]#+0.000000015
      bw3_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f3")
      SLLf3_aic<-Smat(z3,bw3_aic)
      SLLf3_gcv<-Smat(z3,bw3_gcv)
      SLLf3_cp<-Smat(z3,bw3_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_aic%*%(y-x%*%betaLL_aic-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_aic<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #GCV--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_gcv      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_gcv<-(f3p[,(i+1)])
        i<-i+1 
      }
      
      #Cp--------------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3p[,i]-f4)
        f3p[,(i+1)] <- SLLf3_cp%*%(y-x%*%betaLL_cp-f1-f2-f4)
        e[i]  <- mean(abs((f3p[,(i+1)]-f3p[,i])))
        tol<-e[i]
        f3hat_cp<-(f3p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      #BW SELECTION WITH f4 FUNCTION----------------------------------------------
      h<-seq(0.00000000005,0.00000005,length.out = n)
      gcv<-0
      aic<-0
      cp <-0
      for (j in 1:n){
        SLL<-Smat(z4,bw=h[j])
        p<-2
        f0<-matrix(0,n,p)
        f1p<-matrix(0,n,500)
        f2p<-matrix(0,n,500)
        f3p<-matrix(0,n,500)
        f4p<-matrix(0,n,500)
        tol<-99
        e<-0
        i<-1
        while (tol>0.005){
          betaLL      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
          f4p[,(i+1)] <- SLL%*%(y-x%*%betaLL-f1-f2-f3)
          e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
          tol<-e[i]
          f4hat<-(f4p[,(i+1)])
          i<-i+1 
          xt <- (diag(n)-SLL)%*%x
          Hs4 <- xt%*%solve(t(xt)%*%xt)%*%t(xt)
        }
        aic[j]<-aiccfunc(f4,f4hat,Hs4)
        gcv[j]<-gcvfunc(f4,f4hat,Hs4)
        cp[j]<-cpfunc(f4,f4hat,Hs4)
      }
      bw4_aic<-h[which.min(aic)]
      bw4_gcv<-h[which.min(gcv)]#+0.00000002
      bw4_cp<-h[which.min(cp)]
      #plot(aic,type="l", main="AIC scores for bandwidth for f4")
      SLLf4_aic<-Smat(z4,bw4_aic)
      SLLf4_gcv<-Smat(z4,bw4_gcv)
      SLLf4_cp<-Smat(z4,bw4_cp)
      #--------------------------------------
      #(f3) LLR estimator for SPAM Part (Iterative-Backfitting)
      #AIC-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_aic      <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_aic%*%(y-x%*%betaLL_aic-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_aic<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #GCV-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      
      while (tol>0.005){
        betaLL_gcv  <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_gcv%*%(y-x%*%betaLL_gcv-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_gcv<-(f4p[,(i+1)])
        i<-i+1 
      }
      
      #Cp-----------------------------------------------------------------------
      p<-2
      f0<-matrix(0,n,p)
      f1p<-matrix(0,n,500)
      f2p<-matrix(0,n,500)
      f3p<-matrix(0,n,500)
      f4p<-matrix(0,n,500)
      tol<-99
      e<-0
      i<-1
      while (tol>0.005){
        betaLL_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-f1-f2-f3-f4p[,i])
        f4p[,(i+1)] <- SLLf4_cp%*%(y-x%*%betaLL_cp-f1-f2-f3)
        e[i]  <- mean(abs((f4p[,(i+1)]-f4p[,i])))
        tol<-e[i]
        f4hat_cp<-(f4p[,(i+1)])
        i<-i+1 
      }
      #-------------------------------------------------------------------------------
      H <- Hs4
      yhatBF_aic <- x%*%betaLL_aic+f1hat_aic+f2hat_aic+f3hat_aic+f4hat_aic
      yhatBF_gcv <- x%*%betaLL_gcv+f1hat_gcv+f2hat_gcv+f3hat_gcv+f4hat_gcv
      yhatBF_cp <- x%*%betaLL_cp+f1hat_cp+f2hat_cp+f3hat_cp+f4hat_cp
      
      FHAT_aic<-matrix(c(f1hat_aic,f2hat_aic,f3hat_aic,f4hat_aic),n,4)
      FHAT_gcv<-matrix(c(f1hat_gcv,f2hat_gcv,f3hat_gcv,f4hat_gcv),n,4)
      FHAT_cp<-matrix(c(f1hat_cp,f2hat_cp,f3hat_cp,f4hat_cp),n,4)
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(f1hat_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(f2hat_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(f3hat_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(f4hat_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    res<-new.env()
    #Results obtained from backfitting algorithm
    res$beta_aic      <-betaLL_aic
    res$fhat_aic      <-FHAT_aic
    res$fhat_gcv      <-FHAT_gcv
    res$fhat_cp       <-FHAT_cp
    res$fitted_aic    <-yhatBF_aic
    res$fitted_gcv    <-yhatBF_gcv
    res$fitted_cp     <-yhatBF_cp
    res$Smat          <-SLL 
    res$H             <-H 
    return(res)
  }     #REVISED
  
 
  LLR_estimates<-SPAMLL(dat$x,dat$y,dat$nc,dat$allf)
  KS_estimates<-SPAMKS(dat$x,dat$y,dat$nc,dat$allf)
  SS_estimates<-SPAMSS(dat$x,dat$y,dat$nc,dat$allf)
  #-------------------------------------------------
  #FUNCTIONS OF EVALUATION METRICS------------------------------------------------
  #Functions
  #arg: model object from SPAM functions
  biasfunc<-function(arg,x,fhat){
    sumf<-rowSums(fhat)
    n <- length(sumf)
    xtil<-(diag(n)-arg$Smat)%*%x 
    ftil<-(diag(n)-arg$Smat)%*%sumf
    bias<-solve(t(x)%*%xtil)%*%t(x)%*%ftil
    return(bias)
  }
  varfunc <-function(arg,x,y,yhat,H){
    n    <-length(y)
    xtil<-(diag(n)-arg$Smat)%*%x
    sig2<-sum((y-yhat)^2)/((tr(diag(n)-H))^2)
    var  <-sig2*solve(t(x)%*%xtil)%*%t(x)%*%((diag(n)-arg$Smat)^2)%*%x%*%solve(t(x)%*%xtil)
    return(var)
  }
  smdefunc<-function(var,bias){
    smde<-sum(diag(var)+bias^2)
    return(smde)
  }
  REfunc<-function(smde1,smde2){
    re<-smde1/smde2
    return(re)
  }
  msefunc<-function(arg,allf,q,fhat){
    sim <- 1000
    MSE<-0
    for (i in 1:q){
      MSE[i]<-sum((allf[,i]-fhat[,i])^2)/sim
    }
    return(MSE)
  }
  amsefunc<-function(arg,allf,q,fhat){
    sim <- 1000
    sumfhat<-rowSums(fhat)
    sumf   <-rowSums(allf)
    amse   <-(sum((sumfhat-sumf)^2)/q)/sim 
    return(amse)
  }
  ###############################################################################
  #BIASES------------------------------------------------------------------------
  bias.LL_AIC <- biasfunc(LLR_estimates,dat$x,LLR_estimates$fhat_aic)
  bias.LL_GCV <- biasfunc(LLR_estimates,dat$x,LLR_estimates$fhat_gcv)
  bias.LL_Cp <- biasfunc(LLR_estimates,dat$x,LLR_estimates$fhat_cp)
  #----------------------------------------------------------------------
  bias.KS_AIC <- biasfunc(KS_estimates,dat$x,KS_estimates$fhat_aic)
  bias.KS_GCV <- biasfunc(KS_estimates,dat$x,KS_estimates$fhat_gcv)
  bias.KS_Cp <- biasfunc(KS_estimates,dat$x,KS_estimates$fhat_cp)
  #----------------------------------------------------------------------
  bias.SS_AIC <- biasfunc(SS_estimates,dat$x,SS_estimates$fhat_aic)
  bias.SS_GCV <- biasfunc(SS_estimates,dat$x,SS_estimates$fhat_gcv)
  bias.SS_Cp <- biasfunc(SS_estimates,dat$x,SS_estimates$fhat_cp)
  #-----------------------------------------------------------------------
  size<-dim(dat$nc)
  q<-size[2]
  p<-2
  #VARIANCES---------------------------------------------------------------------
  var.LL_AIC <- varfunc(LLR_estimates,dat$x,dat$y,LLR_estimates$fitted_aic,LLR_estimates$H)
  var.LL_GCV <- varfunc(LLR_estimates,dat$x,dat$y,LLR_estimates$fitted_gcv,LLR_estimates$H)
  var.LL_Cp <- varfunc(LLR_estimates,dat$x,dat$y,LLR_estimates$fitted_cp,LLR_estimates$H)
  #-----------------------------------------------------------------------------
  var.KS_AIC <- varfunc(KS_estimates,dat$x,dat$y,KS_estimates$fitted_aic,KS_estimates$H)
  var.KS_GCV <- varfunc(KS_estimates,dat$x,dat$y,KS_estimates$fitted_gcv,KS_estimates$H)
  var.KS_Cp <- varfunc(KS_estimates,dat$x,dat$y,KS_estimates$fitted_cp,KS_estimates$H)
  #------------------------------------------------------------------------------
  var.SS_AIC <- varfunc(SS_estimates,dat$x,dat$y,SS_estimates$fitted_aic,SS_estimates$H)
  var.SS_GCV <- varfunc(SS_estimates,dat$x,dat$y,SS_estimates$fitted_gcv,SS_estimates$H)
  var.SS_Cp <- varfunc(SS_estimates,dat$x,dat$y,SS_estimates$fitted_cp,SS_estimates$H)
  #SMDES-------------------------------------------------------------------------
  SMDE.LL_AIC<-smdefunc(var.LL_AIC,bias.LL_AIC)
  SMDE.LL_GCV<-smdefunc(var.LL_GCV,bias.LL_GCV)
  SMDE.LL_Cp<-smdefunc(var.LL_Cp,bias.LL_Cp)
#--------------------------------------------------------------------------------
  SMDE.KS_AIC<-smdefunc(var.KS_AIC,bias.KS_AIC)
  SMDE.KS_GCV<-smdefunc(var.KS_GCV,bias.KS_GCV)
  SMDE.KS_Cp<-smdefunc(var.KS_Cp,bias.KS_Cp)
#---------------------------------------------------------------------------------
  SMDE.SS_AIC<-smdefunc(var.SS_AIC,bias.SS_AIC)
  SMDE.SS_GCV<-smdefunc(var.SS_GCV,bias.SS_GCV)
  SMDE.SS_Cp<-smdefunc(var.SS_Cp,bias.SS_Cp)
#Relative Efficiencies (RE)-----------------------------------------------------
  RE.LL_AIC <-c(REfunc(SMDE.LL_AIC,SMDE.KS_AIC),REfunc(SMDE.LL_AIC,SMDE.SS_AIC))    #LL/KS, LL/SS
  RE.LL_GCV <-c(REfunc(SMDE.LL_GCV,SMDE.KS_GCV),REfunc(SMDE.LL_GCV,SMDE.SS_GCV))    #LL/KS, LL/SS
  RE.LL_Cp  <-c(REfunc(SMDE.LL_Cp,SMDE.KS_Cp),REfunc(SMDE.LL_Cp,SMDE.SS_Cp))        #LL/KS, LL/SS
#---------------------------------------------------------------------------------  
  RE.KS_AIC <-c(REfunc(SMDE.KS_AIC,SMDE.LL_AIC),REfunc(SMDE.KS_AIC,SMDE.SS_AIC))    #KS/LL, KS/SS
  RE.KS_GCV <-c(REfunc(SMDE.KS_GCV,SMDE.LL_GCV),REfunc(SMDE.KS_GCV,SMDE.SS_GCV))    #KS/LL, KS/SS
  RE.KS_Cp  <-c(REfunc(SMDE.KS_Cp,SMDE.LL_Cp),REfunc(SMDE.KS_Cp,SMDE.SS_Cp))        #KS/LL, KS/SS
#---------------------------------------------------------------------------------  
  RE.SS_AIC  <-c(REfunc(SMDE.SS_AIC,SMDE.LL_AIC),REfunc(SMDE.SS_AIC,SMDE.KS_AIC))   #SS/LL, SS/KS
  RE.SS_GCV  <-c(REfunc(SMDE.SS_GCV,SMDE.LL_GCV),REfunc(SMDE.SS_GCV,SMDE.KS_GCV))   #SS/LL, SS/KS
  RE.SS_Cp   <-c(REfunc(SMDE.SS_Cp,SMDE.LL_Cp),REfunc(SMDE.SS_Cp,SMDE.KS_Cp))       #SS/LL, SS/KS
#MSE----------------------------------------------------------------------------
  MSE.LL_AIC <-msefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_aic)
  MSE.LL_GCV <-msefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_gcv)
  MSE.LL_Cp  <-msefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_cp)
#--------------------------------------------------------------------------------
  MSE.KS_AIC <-msefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_aic)
  MSE.KS_GCV <-msefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_gcv)
  MSE.KS_Cp  <-msefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_cp)
#--------------------------------------------------------------------------------
  MSE.SS_AIC <-msefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_aic)
  MSE.SS_GCV <-msefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_gcv)
  MSE.SS_Cp  <-msefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_cp)
#AMSE----------------------------------------------------------------------------
  AMSE.LL_AIC  <-amsefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_aic)
  AMSE.LL_GCV  <-amsefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_gcv)
  AMSE.LL_Cp  <-amsefunc(LLR_estimates,dat$allf,q,LLR_estimates$fhat_cp)
#--------------------------------------------------------------------------------
  AMSE.KS_AIC  <-amsefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_aic)
  AMSE.KS_GCV  <-amsefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_gcv)
  AMSE.KS_Cp   <-amsefunc(KS_estimates,dat$allf,q,KS_estimates$fhat_cp)
#--------------------------------------------------------------------------------
  AMSE.SS_AIC  <-amsefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_aic)
  AMSE.SS_GCV  <-amsefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_gcv)
  AMSE.SS_Cp   <-amsefunc(SS_estimates,dat$allf,q,SS_estimates$fhat_cp)
  
  
  ret1 <- c(BiasLL.AIC=bias.LL_AIC,BiasLL.GCV=bias.LL_GCV, BiasLL.Cp=bias.LL_Cp, BiasKS.AIC=bias.KS_AIC,BiasKS.GCV=bias.KS_GCV,BiasKS.Cp=bias.KS_Cp,BiasSS.AIC=bias.SS_AIC,BiasSS.GCV=bias.SS_GCV,BiasSS.Cp=bias.SS_Cp,varLL.AIC=var.LL_AIC,varLL.GCV=var.LL_GCV,varLL.Cp=var.LL_Cp,varKS.AIC=var.KS_AIC,varKS.GCV=var.KS_GCV,varKS.Cp=var.KS_Cp,varSS.AIC=var.SS_AIC, varSS.GCV=var.SS_GCV,varSS.Cp=var.SS_Cp)
  ret2 <-c(SMDELL.AIC=SMDE.LL_AIC,SMDELL.GCV=SMDE.LL_GCV,SMDELL.Cp=SMDE.LL_Cp,SMDEKS.AIC=SMDE.KS_AIC,SMDEKS.GCV=SMDE.KS_GCV,SMDEKS.Cp=SMDE.KS_Cp,SMDESS.AIC=SMDE.SS_AIC,SMDESS.GCV=SMDE.SS_GCV,SMDESS.Cp=SMDE.SS_Cp) 
  ret3 <- c(RELL.AIC=RE.LL_AIC,RELL.GCV=RE.LL_GCV,RELL.Cp=RE.LL_Cp,REKS.AIC=RE.KS_AIC,REKS.GCV=RE.KS_GCV,REKS.Cp=RE.KS_Cp,RESS.AIC=RE.SS_AIC,RESS.GCV=RE.SS_GCV,RESS.Cp=RE.SS_Cp)
  ret4 <- c(MSELL.AIC=MSE.LL_AIC,MSELL.GCV=MSE.LL_GCV,MSELL.Cp=MSE.LL_Cp,MSEKS.AIC=MSE.KS_AIC,MSEKS.GCV=MSE.KS_GCV,MSEKS.Cp=MSE.KS_Cp,MSESS.AIC=MSE.SS_AIC,MSESS.GCV=MSE.SS_GCV,MSESS.Cp=MSE.SS_Cp)
  ret5 <- c(AMSELL.AIC=AMSE.LL_AIC,AMSELL.GCV=AMSE.LL_GCV,AMSELL.Cp=AMSE.LL_Cp,AMSEKS.AIC=AMSE.KS_AIC,AMSEKS.GCV=AMSE.KS_GCV,AMSEKS.Cp=AMSE.KS_Cp,AMSESS.AIC=AMSE.SS_AIC,AMSESS.GCV=AMSE.SS_GCV,AMSESS.Cp=AMSE.SS_Cp)
 
  ret <-c(ret1,ret2,ret3,ret4,ret5)
  ret
#------------------------------------------------------------------------------------
  
  }

res <- runSimulation(des, replications = 1000, generate=Generate, 
                     analyse=Analyse)
