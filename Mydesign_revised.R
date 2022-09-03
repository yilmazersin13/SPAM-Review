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
    #----------------------------------
    K<-function(z){
      res<-(1/sqrt(2*pi))*exp(-0.5*z^2)
      return(res)
    }
    #---------------------------------
    
    if (q==2){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-allf[,1]
      f2<-allf[,2]
      n<-length(y)
      z1<-nc[,1]
      z2<-nc[,2]
      #BACKFITTING PROCEDURE----------------------------------------------------
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2))
      #SMOOTHING MATRIX FOR KS--------------------------------------------------
      Smat<-function(z,bw){
        library(condSURV)
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
      #-----------------------------------------------------------------------------
      selectionLL <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
        Smat<-function(z,bw){
          library(condSURV)
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.05,0.9,length.out=50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      
      fhat_aic <- matrix(0,n,2)
      fhat_gcv <- matrix(0,n,2)
      fhat_cp  <- matrix(0,n,2)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      for (k in 1:2){
        tp <- selectionLL(x,nc[,k],(y))
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        tol    <- 0.05
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+fhat_aic[,-k])) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-fhat_aic[,-k])
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+fhat_gcv[,-k])) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-fhat_gcv[,-k])
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+fhat_cp[,-k])) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-fhat_cp[,-k])
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #  plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #  par(new=TRUE)
      #  plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #  par(new=TRUE)
      #  plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #  par(new=TRUE)
      #  plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      # plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #  par(new=TRUE)
      #  plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #  par(new=TRUE)
      #  plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #  par(new=TRUE)
      #  plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic),n,2)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv),n,2)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp),n,2)
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
      
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2))
      f03 <- fitted(lm(y~z3))
      f04 <- fitted(lm(y~z4))
      
      #SMOOTHING MATRIX FOR KS--------------------------------------------------
      Smat<-function(z,bw){
        library(condSURV)
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
      #-----------------------------------------------------------------------------
      selectionLL <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
        Smat<-function(z,bw){
          library(condSURV)
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.05,0.9,length.out=50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      fhat3_aic <- f03 #matrix(0,n,1)
      fhat4_aic <- f04 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      fhat3_gcv <- f03
      fhat4_gcv <- f04
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      fhat3_cp <- f03
      fhat4_cp <- f04
      
      fhat_aic <- matrix(0,n,4)
      fhat_gcv <- matrix(0,n,4)
      fhat_cp  <- matrix(0,n,4)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      fh3_aic <- matrix(0,n,iter)
      fh4_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      fh3_gcv <- matrix(0,n,iter)
      fh4_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      fh3_cp <- matrix(0,n,iter)
      fh4_cp <- matrix(0,n,iter)
      ones <- matrix(1,n,1)
      x <- matrix(c(ones,x),n,3)
      for (k in 1:4){
        tp <- selectionLL(x,nc[,k],y)
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        ones <- matrix(1,n,1)
        tol    <- 0.01
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+rowSums(allf[,-k]))) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-rowSums(allf[,-k]))
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+rowSums(allf[,-k]))) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-rowSums(allf[,-k]))
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+rowSums(allf[,-k]))) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-rowSums(allf[,-k]))
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (k==3){
            fhat3_aic <- fhat_aic[,k]
            fhat3_gcv <- fhat_gcv[,k]
            fhat3_cp  <- fhat_cp[,k]
            
            fh3_aic[,i]   <- fhat3_aic
            fh3_gcv[,i]   <- fhat3_gcv
            fh3_cp[,i]    <- fhat3_cp
          }
          if (k==4){
            fhat4_aic <- fhat_aic[,k]
            fhat4_gcv <- fhat_gcv[,k]
            fhat4_cp  <- fhat_cp[,k]
            
            fh4_aic[,i]   <- fhat4_aic
            fh4_gcv[,i]   <- fhat4_gcv
            fh4_cp[,i]    <- fhat4_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      #-------------------------------------------------------------------------------
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic+fhat3_aic+fhat4_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv+fhat3_gcv+fhat4_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp+fhat3_cp+fhat4_cp
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic,fhat3_aic,fhat4_aic),n,4)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv,fhat3_gcv,fhat4_gcv),n,4)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp,fhat3_cp,fhat4_cp),n,4)
      
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    
    
    res<-new.env()
    #----------------------------------------------
    res$beta_aic   <-betai_aic
    res$fhat_aic   <-FHAT_aic
    res$beta_gcv   <-betai_gcv
    res$fhat_gcv   <-FHAT_gcv
    res$beta_cp    <-betai_cp
    res$fhat_cp    <-FHAT_cp
    res$fitted_aic <-yhatBF_aic
    res$fitted_gcv <-yhatBF_gcv
    res$fitted_cp  <-yhatBF_cp
    res$Smat.aic   <-S_aic
    res$Smat.gcv   <-S_gcv
    res$Smat.cp    <-S_cp
    res$H.aic      <-H.aic
    res$H.gcv      <-H.gcv
    res$H.cp       <-H.cp
    
    return(res)
  } 
  SPAMKS<-function(x,y,nc,allf){
    #x: matrix of parametric covariates
    #y: response variable
    #nc: Matrix of nonparamtric covariate
    library(optR)
    library(psych)
    library(pracma)
    size<-dim(nc)
    q<-size[2]
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
      
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2)) 
      #------------------------------------------------------------------------
      Smat<-function(z,bw){
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
      selectionKS <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
        Smat<-function(z,bw){
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.1,0.9,length.out=50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      
      fhat_aic <- matrix(0,n,2)
      fhat_gcv <- matrix(0,n,2)
      fhat_cp  <- matrix(0,n,2)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      for (k in 1:2){
        tp <- selectionKS(x,nc[,k],(y))
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        tol    <- 0.05
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+fhat_aic[,-k])) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-fhat_aic[,-k])
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+fhat_gcv[,-k])) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-fhat_gcv[,-k])
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+fhat_cp[,-k])) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-fhat_cp[,-k])
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic),n,2)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv),n,2)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp),n,2)
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
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2))
      f03 <- fitted(lm(y~z3))
      f04 <- fitted(lm(y~z4))
      
      #SMOOTHING MATRIX FOR KS--------------------------------------------------
      Smat<-function(z,bw){
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
      #-----------------------------------------------------------------------------
      selectionKS <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
        Smat<-function(z,bw){
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.1,0.9,length.out=50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      fhat3_aic <- f03 #matrix(0,n,1)
      fhat4_aic <- f04 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      fhat3_gcv <- f03
      fhat4_gcv <- f04
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      fhat3_cp <- f03
      fhat4_cp <- f04
      
      fhat_aic <- matrix(0,n,4)
      fhat_gcv <- matrix(0,n,4)
      fhat_cp  <- matrix(0,n,4)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      fh3_aic <- matrix(0,n,iter)
      fh4_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      fh3_gcv <- matrix(0,n,iter)
      fh4_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      fh3_cp <- matrix(0,n,iter)
      fh4_cp <- matrix(0,n,iter)
      ones <- matrix(1,n,1)
      x <- matrix(c(ones,x),n,3)
      for (k in 1:4){
        tp <- selectionKS(x,nc[,k],y)
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        ones <- matrix(1,n,1)
        tol    <- 0.01
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+rowSums(allf[,-k]))) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-rowSums(allf[,-k]))
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+rowSums(allf[,-k]))) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-rowSums(allf[,-k]))
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+rowSums(allf[,-k]))) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-rowSums(allf[,-k]))
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (k==3){
            fhat3_aic <- fhat_aic[,k]
            fhat3_gcv <- fhat_gcv[,k]
            fhat3_cp  <- fhat_cp[,k]
            
            fh3_aic[,i]   <- fhat3_aic
            fh3_gcv[,i]   <- fhat3_gcv
            fh3_cp[,i]    <- fhat3_cp
          }
          if (k==4){
            fhat4_aic <- fhat_aic[,k]
            fhat4_gcv <- fhat_gcv[,k]
            fhat4_cp  <- fhat_cp[,k]
            
            fh4_aic[,i]   <- fhat4_aic
            fh4_gcv[,i]   <- fhat4_gcv
            fh4_cp[,i]    <- fhat4_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      #-------------------------------------------------------------------------------
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic+fhat3_aic+fhat4_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv+fhat3_gcv+fhat4_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp+fhat3_cp+fhat4_cp
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic,fhat3_aic,fhat4_aic),n,4)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv,fhat3_gcv,fhat4_gcv),n,4)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp,fhat3_cp,fhat4_cp),n,4)
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    res<-new.env()
    #----------------------------------------------
    res$beta_aic   <-betai_aic
    res$fhat_aic   <-FHAT_aic
    res$beta_gcv   <-betai_gcv
    res$fhat_gcv   <-FHAT_gcv
    res$beta_cp    <-betai_cp
    res$fhat_cp    <-FHAT_cp
    res$fitted_aic <-yhatBF_aic
    res$fitted_gcv <-yhatBF_gcv
    res$fitted_cp  <-yhatBF_cp
    res$Smat.aic   <-S_aic
    res$Smat.gcv   <-S_gcv
    res$Smat.cp    <-S_cp
    res$H.aic      <-H.aic
    res$H.gcv      <-H.gcv
    res$H.cp       <-H.cp
    return(res)
  }
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
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2)) 
      #------------------------------------------------------------------------
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
      selectionSS <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.000000000005,0.000000005,length.out = 50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      
      fhat_aic <- matrix(0,n,2)
      fhat_gcv <- matrix(0,n,2)
      fhat_cp  <- matrix(0,n,2)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      for (k in 1:2){
        tp <- selectionSS(x,nc[,k],(y))
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        tol    <- 0.05
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+fhat_aic[,-k])) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-fhat_aic[,-k])
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+fhat_gcv[,-k])) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-fhat_gcv[,-k])
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+fhat_cp[,-k])) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-fhat_cp[,-k])
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic),n,2)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv),n,2)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp),n,2)
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
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1))
      f02 <- fitted(lm(y~z2))
      f03 <- fitted(lm(y~z3))
      f04 <- fitted(lm(y~z4))
      
      #SMOOTHING MATRIX FOR KS--------------------------------------------------
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
      selectionSS <- function(x,z,y){
        aic <- 0
        gcv <- 0
        cp  <- 0
        aiccfunc<-function(y,yhat,H){
          p <- tr(H)
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
          return(score)
        }
        #---------------GCV-----------------
        gcvfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
          return(score)
        }
        #---------------GCV-----------------
        cpfunc<-function(y,yhat,H){
          y<-matrix(c(y))
          yhat<-matrix(c(yhat))
          n<-length(y)
          MSE<-norm(((diag(n)-H)%*%y)^2)
          ssq<-var(y-yhat)
          DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
          score<-MSE+(2*ssq*DF)
          return(score)
        }
        #-----------------------------------------------------------------------------
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
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.000000000005,0.000000005,length.out = 50)
        for (i in 1:50){
          W    <- Smat(z,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil)%*%t(x)%*%(diag(n)-W)^2
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:50){
          if (aic[i2]==min(aic)){
            lam_aic <- tp_seq[i2]
          }
          if (gcv[i2]==min(gcv)){
            lam_gcv <- tp_seq[i2]
          }
          if (cp[i2]==min(cp)){
            lam_cp <- tp_seq[i2]
          }
        }
        res <- new.env()
        res$lam.aic <- lam_aic
        res$lam.gcv <- lam_gcv
        res$lam.cp  <- lam_cp
        return(res)
      }
      #BACKFITTING BEGINS-------------------------------------------------------
      #-------------------------------------------------------------------------
      alpha0 <- mean(y)
      fhat1_aic <- f01 #matrix(0,n,1)
      fhat2_aic <- f02 #matrix(0,n,1)
      fhat3_aic <- f03 #matrix(0,n,1)
      fhat4_aic <- f04 #matrix(0,n,1)
      
      fhat1_gcv <- f01
      fhat2_gcv <- f02
      fhat3_gcv <- f03
      fhat4_gcv <- f04
      
      fhat1_cp <- f01
      fhat2_cp <- f02
      fhat3_cp <- f03
      fhat4_cp <- f04
      
      fhat_aic <- matrix(0,n,4)
      fhat_gcv <- matrix(0,n,4)
      fhat_cp  <- matrix(0,n,4)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      fh3_aic <- matrix(0,n,iter)
      fh4_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      fh3_gcv <- matrix(0,n,iter)
      fh4_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      fh3_cp <- matrix(0,n,iter)
      fh4_cp <- matrix(0,n,iter)
      ones <- matrix(1,n,1)
      x <- matrix(c(ones,x),n,3)
      for (k in 1:4){
        tp <- selectionSS(x,nc[,k],y)
        tp.aic <- tp$lam.aic
        tp.gcv <- tp$lam.gcv
        tp.cp  <- tp$lam.cp
        S_aic  <- Smat(nc[,k],tp.aic) 
        S_gcv  <- Smat(nc[,k],tp.gcv)
        S_cp   <- Smat(nc[,k],tp.cp)
        
        xtil.aic <- (diag(n)-S_aic)%*%x
        xtil.gcv <- (diag(n)-S_gcv)%*%x
        xtil.cp  <- (diag(n)-S_cp)%*%x
        
        H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)^2
        H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)^2
        H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)^2
        
        ones <- matrix(1,n,1)
        tol    <- 0.01
        ctol   <- 99
        i <- 1
        while (ctol>=tol){
          
          betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+rowSums(allf[,-k]))) 
          fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-rowSums(allf[,-k]))
          
          betai_gcv   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_gcv[,k]+rowSums(allf[,-k]))) 
          fhat_gcv[,k]<- S_gcv%*%(y-alpha0-(x%*%betai_gcv)-rowSums(allf[,-k]))
          
          betai_cp   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_cp[,k]+rowSums(allf[,-k]))) 
          fhat_cp[,k]<- S_cp%*%(y-alpha0-(x%*%betai_cp)-rowSums(allf[,-k]))
          
          if (k==1){
            fhat1_aic <- fhat_aic[,k]
            fhat1_gcv <- fhat_gcv[,k]
            fhat1_cp  <- fhat_cp[,k]
            
            fh1_aic[,i]   <- fhat1_aic
            fh1_gcv[,i]   <- fhat1_gcv
            fh1_cp[,i]    <- fhat1_cp
          }
          if (k==2){
            fhat2_aic <- fhat_aic[,k]
            fhat2_gcv <- fhat_gcv[,k]
            fhat2_cp  <- fhat_cp[,k]
            
            fh2_aic[,i]   <- fhat2_aic
            fh2_gcv[,i]   <- fhat2_gcv
            fh2_cp[,i]    <- fhat2_cp
          }
          if (k==3){
            fhat3_aic <- fhat_aic[,k]
            fhat3_gcv <- fhat_gcv[,k]
            fhat3_cp  <- fhat_cp[,k]
            
            fh3_aic[,i]   <- fhat3_aic
            fh3_gcv[,i]   <- fhat3_gcv
            fh3_cp[,i]    <- fhat3_cp
          }
          if (k==4){
            fhat4_aic <- fhat_aic[,k]
            fhat4_gcv <- fhat_gcv[,k]
            fhat4_cp  <- fhat_cp[,k]
            
            fh4_aic[,i]   <- fhat4_aic
            fh4_gcv[,i]   <- fhat4_gcv
            fh4_cp[,i]    <- fhat4_cp
          }
          if (i>1){
            tol2      <-(mean(abs(fh2_cp[,(i-1)]-fhat2_cp)))
            tol1      <-(mean(abs(fh1_gcv[,(i-1)]-fhat1_gcv)))
            ctol <- (tol1+tol2)/2
          }
          i <- i+1
          if (i==iter){
            break
          }
        }
      }
      #-------------------------------------------------------------------------------
      
      yhatBF_aic <- x%*%betai_aic+fhat1_aic+fhat2_aic+fhat3_aic+fhat4_aic
      yhatBF_gcv <- x%*%betai_gcv+fhat1_gcv+fhat2_gcv+fhat3_gcv+fhat4_gcv
      yhatBF_cp <- x%*%betai_cp+fhat1_cp+fhat2_cp+fhat3_cp+fhat4_cp
      
      FHAT_aic<-matrix(c(fhat1_aic,fhat2_aic,fhat3_aic,fhat4_aic),n,4)
      FHAT_gcv<-matrix(c(fhat1_gcv,fhat2_gcv,fhat3_gcv,fhat4_gcv),n,4)
      FHAT_cp<-matrix(c(fhat1_cp,fhat2_cp,fhat3_cp,fhat4_cp),n,4)
      
      #plot(f1,pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_aic,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_gcv,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      #par(new=TRUE)
      #plot(fhat1_cp,col="red",pch=19,type="l",ylim=c(min(f1),max(f1)))
      
      #plot(f2,pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_aic,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_gcv,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      #par(new=TRUE)
      #plot(fhat2_cp,col="red",pch=19,type="l",ylim=c(min(f2),max(f2)))
      
      #plot(f3,pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_aic,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_gcv,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      #par(new=TRUE)
      #plot(fhat3_cp,col="red",pch=19,type="l",ylim=c(min(f3),max(f3)))
      
      #plot(f4,pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_aic,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_gcv,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      #par(new=TRUE)
      #plot(fhat4_cp,col="red",pch=19,type="l",ylim=c(min(f4),max(f4)))
      
      #plot(y,pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_aic,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_gcv,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(yhatBF_cp,col="red",pch=19,type="l",ylim=c(min(y),max(y)))
    }
    res<-new.env()
    #----------------------------------------------
    res$beta_aic   <-betai_aic
    res$fhat_aic   <-FHAT_aic
    res$beta_gcv   <-betai_gcv
    res$fhat_gcv   <-FHAT_gcv
    res$beta_cp    <-betai_cp
    res$fhat_cp    <-FHAT_cp
    res$fitted_aic <-yhatBF_aic
    res$fitted_gcv <-yhatBF_gcv
    res$fitted_cp  <-yhatBF_cp
    res$Smat.aic   <-S_aic
    res$Smat.gcv   <-S_gcv
    res$Smat.cp    <-S_cp
    res$H.aic      <-H.aic
    res$H.gcv      <-H.gcv
    res$H.cp       <-H.cp
    
    return(res)
  }   
  
  LLR_estimates<-SPAMLL(dat$x,dat$y,dat$nc,dat$allf)
  KS_estimates<-SPAMKS(dat$x,dat$y,dat$nc,dat$allf)
  SS_estimates<-SPAMSS(dat$x,dat$y,dat$nc,dat$allf)
  #-------------------------------------------------
  #FUNCTIONS OF EVALUATION METRICS------------------------------------------------
  #Functions
  #arg: model object from SPAM functions
  biasfunc<-function(Smat,x,fhat){
    sumf<-rowSums(fhat)
    n <- length(sumf)
    xtil<-(diag(n)-Smat)%*%x 
    ftil<-(diag(n)-Smat)%*%sumf
    bias<-solve(t(x)%*%xtil)%*%t(x)%*%ftil
    return(bias)
  }
  varfunc <-function(Smat,x,y,yhat,H){
    n    <-length(y)
    xtil<-(diag(n)-Smat)%*%x
    sig2<-sum((y-yhat)^2)/((tr(diag(n)-H))^2)
    var  <-sig2*solve(t(x)%*%xtil)%*%t(x)%*%((diag(n)-Smat)^2)%*%x%*%solve(t(x)%*%xtil)
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
  msefunc<-function(allf,q,fhat){
    sim <- 1000
    MSE<-0
    for (i in 1:q){
      MSE[i]<-sum((allf[,i]-fhat[,i])^2)/sim
    }
    return(MSE)
  }
  amsefunc<-function(allf,q,fhat){
    sim <- 1000
    sumfhat<-rowSums(fhat)
    sumf   <-rowSums(allf)
    amse   <-(sum((sumfhat-sumf)^2)/q)/sim 
    return(amse)
  }
  ###############################################################################
  #BIASES------------------------------------------------------------------------
  bias.LL_AIC <- biasfunc(LLR_estimates$Smat.aic,dat$x,LLR_estimates$fhat_aic)
  bias.LL_GCV <- biasfunc(LLR_estimates$Smat.gcv,dat$x,LLR_estimates$fhat_gcv)
  bias.LL_Cp <- biasfunc(LLR_estimates$Smat.cp,dat$x,LLR_estimates$fhat_cp)
  #----------------------------------------------------------------------
  bias.KS_AIC <- biasfunc(KS_estimates$Smat.aic,dat$x,KS_estimates$fhat_aic)
  bias.KS_GCV <- biasfunc(KS_estimates$Smat.gcv,dat$x,KS_estimates$fhat_gcv)
  bias.KS_Cp <- biasfunc(KS_estimates$Smat.cp,dat$x,KS_estimates$fhat_cp)
  #----------------------------------------------------------------------
  bias.SS_AIC <- biasfunc(SS_estimates$Smat.aic,dat$x,SS_estimates$fhat_aic)
  bias.SS_GCV <- biasfunc(SS_estimates$Smat.gcv,dat$x,SS_estimates$fhat_gcv)
  bias.SS_Cp <- biasfunc(SS_estimates$Smat.cp,dat$x,SS_estimates$fhat_cp)
  #-----------------------------------------------------------------------
  size<-dim(dat$nc)
  q<-size[2]
  p<-2
  #VARIANCES---------------------------------------------------------------------
  var.LL_AIC <- varfunc(LLR_estimates$Smat.aic,dat$x,dat$y,LLR_estimates$fitted_aic,LLR_estimates$H.aic)
  var.LL_GCV <- varfunc(LLR_estimates$Smat.gcv,dat$x,dat$y,LLR_estimates$fitted_gcv,LLR_estimates$H.gcv)
  var.LL_Cp <- varfunc(LLR_estimates$Smat.cp,dat$x,dat$y,LLR_estimates$fitted_cp,LLR_estimates$H.cp)
  #-----------------------------------------------------------------------------
  var.KS_AIC <- varfunc(KS_estimates$Smat.aic,dat$x,dat$y,KS_estimates$fitted_aic,KS_estimates$H.aic)
  var.KS_GCV <- varfunc(KS_estimates$Smat.gcv,dat$x,dat$y,KS_estimates$fitted_gcv,KS_estimates$H.gcv)
  var.KS_Cp <- varfunc(KS_estimates$Smat.cp,dat$x,dat$y,KS_estimates$fitted_cp,KS_estimates$H.cp)
  #------------------------------------------------------------------------------
  var.SS_AIC <- varfunc(SS_estimates$Smat.aic,dat$x,dat$y,SS_estimates$fitted_aic,SS_estimates$H.aic)
  var.SS_GCV <- varfunc(SS_estimates$Smat.gcv,dat$x,dat$y,SS_estimates$fitted_gcv,SS_estimates$H.gcv)
  var.SS_Cp <- varfunc(SS_estimates$Smat.cp,dat$x,dat$y,SS_estimates$fitted_cp,SS_estimates$H.cp)
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
  MSE.LL_AIC <-msefunc(dat$allf,q,LLR_estimates$fhat_aic)
  MSE.LL_GCV <-msefunc(dat$allf,q,LLR_estimates$fhat_gcv)
  MSE.LL_Cp  <-msefunc(dat$allf,q,LLR_estimates$fhat_cp)
  #--------------------------------------------------------------------------------
  MSE.KS_AIC <-msefunc(dat$allf,q,KS_estimates$fhat_aic)
  MSE.KS_GCV <-msefunc(dat$allf,q,KS_estimates$fhat_gcv)
  MSE.KS_Cp  <-msefunc(dat$allf,q,KS_estimates$fhat_cp)
  #--------------------------------------------------------------------------------
  MSE.SS_AIC <-msefunc(dat$allf,q,SS_estimates$fhat_aic)
  MSE.SS_GCV <-msefunc(dat$allf,q,SS_estimates$fhat_gcv)
  MSE.SS_Cp  <-msefunc(dat$allf,q,SS_estimates$fhat_cp)
  #AMSE----------------------------------------------------------------------------
  AMSE.LL_AIC  <-amsefunc(dat$allf,q,LLR_estimates$fhat_aic)
  AMSE.LL_GCV  <-amsefunc(dat$allf,q,LLR_estimates$fhat_gcv)
  AMSE.LL_Cp  <-amsefunc(dat$allf,q,LLR_estimates$fhat_cp)
  #--------------------------------------------------------------------------------
  AMSE.KS_AIC  <-amsefunc(dat$allf,q,KS_estimates$fhat_aic)
  AMSE.KS_GCV  <-amsefunc(dat$allf,q,KS_estimates$fhat_gcv)
  AMSE.KS_Cp   <-amsefunc(dat$allf,q,KS_estimates$fhat_cp)
  #--------------------------------------------------------------------------------
  AMSE.SS_AIC  <-amsefunc(dat$allf,q,SS_estimates$fhat_aic)
  AMSE.SS_GCV  <-amsefunc(dat$allf,q,SS_estimates$fhat_gcv)
  AMSE.SS_Cp   <-amsefunc(dat$allf,q,SS_estimates$fhat_cp)
  
  ret1 <- c(BiasLL.AIC=bias.LL_AIC,BiasLL.GCV=bias.LL_GCV, BiasLL.Cp=bias.LL_Cp, BiasKS.AIC=bias.KS_AIC,BiasKS.GCV=bias.KS_GCV,BiasKS.Cp=bias.KS_Cp,BiasSS.AIC=bias.SS_AIC,BiasSS.GCV=bias.SS_GCV,BiasSS.Cp=bias.SS_Cp,varLL.AIC=var.LL_AIC,varLL.GCV=var.LL_GCV,varLL.Cp=var.LL_Cp,varKS.AIC=var.KS_AIC,varKS.GCV=var.KS_GCV,varKS.Cp=var.KS_Cp,varSS.AIC=var.SS_AIC, varSS.GCV=var.SS_GCV,varSS.Cp=var.SS_Cp)
  ret2 <-c(SMDELL.AIC=SMDE.LL_AIC,SMDELL.GCV=SMDE.LL_GCV,SMDELL.Cp=SMDE.LL_Cp,SMDEKS.AIC=SMDE.KS_AIC,SMDEKS.GCV=SMDE.KS_GCV,SMDEKS.Cp=SMDE.KS_Cp,SMDESS.AIC=SMDE.SS_AIC,SMDESS.GCV=SMDE.SS_GCV,SMDESS.Cp=SMDE.SS_Cp) 
  ret3 <- c(RELL.AIC=RE.LL_AIC,RELL.GCV=RE.LL_GCV,RELL.Cp=RE.LL_Cp,REKS.AIC=RE.KS_AIC,REKS.GCV=RE.KS_GCV,REKS.Cp=RE.KS_Cp,RESS.AIC=RE.SS_AIC,RESS.GCV=RE.SS_GCV,RESS.Cp=RE.SS_Cp)
  ret4 <- c(MSELL.AIC=MSE.LL_AIC,MSELL.GCV=MSE.LL_GCV,MSELL.Cp=MSE.LL_Cp,MSEKS.AIC=MSE.KS_AIC,MSEKS.GCV=MSE.KS_GCV,MSEKS.Cp=MSE.KS_Cp,MSESS.AIC=MSE.SS_AIC,MSESS.GCV=MSE.SS_GCV,MSESS.Cp=MSE.SS_Cp)
  ret5 <- c(AMSELL.AIC=AMSE.LL_AIC,AMSELL.GCV=AMSE.LL_GCV,AMSELL.Cp=AMSE.LL_Cp,AMSEKS.AIC=AMSE.KS_AIC,AMSEKS.GCV=AMSE.KS_GCV,AMSEKS.Cp=AMSE.KS_Cp,AMSESS.AIC=AMSE.SS_AIC,AMSESS.GCV=AMSE.SS_GCV,AMSESS.Cp=AMSE.SS_Cp)
  
  ret <-c(ret1,ret2,ret3,ret4,ret5)
  ret
  #------------------------------------------------------------------------------------
  
}

res <- runSimulation(des, replications = 2, generate=Generate, 
                     analyse=Analyse)
