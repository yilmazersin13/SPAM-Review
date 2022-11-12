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
      z1<-scale(nc[,1])     # After Revision 3
      z2<-scale(nc[,2])     # After Revision 3
      #BACKFITTING PROCEDURE----------------------------------------------------
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
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
      
      tol    <- 0.05
      ctol   <- 99
      i <- 1
      while (ctol>=tol){
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
      z1<-scale(nc[,1])
      z2<-scale(nc[,3])
      z3<-scale(nc[,2])
      z4<-scale(nc[,4])
      
      #------------------------------------
      #LLR estimator for SPAM Part (Non-iterative)
      #Selection of Bandwidth---------------------------------------------------------
      
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
      f03 <- fitted(lm(y~z3+x[,1]+x[,2]))
      f04 <- fitted(lm(y~z4+x[,1]+x[,2]))
      
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
      
      ones <- matrix(1,n,1)
      tol    <- 0.01
      ctol   <- 99
      i <- 1
      while (ctol>=tol){
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
      z1<-scale(nc[,1])
      z2<-scale(nc[,2])
      
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
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
      
      
      tol    <- 0.05
      ctol   <- 99
      i <- 1
      while (ctol>=tol){
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
      z1<-scale(nc[,1])
      z2<-scale(nc[,3])
      z3<-scale(nc[,2])
      z4<-scale(nc[,4])
      
      #---------------------------------------------------------------------------
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
      f03 <- fitted(lm(y~z3+x[,1]+x[,2]))
      f04 <- fitted(lm(y~z4+x[,1]+x[,2]))
      
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
      
      
      ones <- matrix(1,n,1)
      tol    <- 0.01
      ctol   <- 99
      i <- 1
      while (ctol>=tol){
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
    Smat = function(x, df,lam){
      n = length(x);
      A = matrix(0, n, n);
      for(i in 1:n){
        y = rep(0, n); y[i]=1;
        yi = smooth.spline(x, y, df=df,spar=lam)$y;
        A[,i]= yi;
      }
      return(A)
    } 
    if (q==2){
      #x: matrix of parametric covariates
      #y: response variable
      #nc: Matrix of nonparamtric covariate
      f1<-(allf[,1])
      f2<-scale(allf[,2])
      n<-length(y)
      
      z1<-scale(nc[,1])
      z2<-scale(nc[,2])
      
      #BW SELECTION WITH f1 FUNCTION----------------------------------------------
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
      
      #------------------------------------------------------------------------
      Smat = function(x, df,lam){
        n = length(x);
        A = matrix(0, n, n);
        for(i in 1:n){
          y = rep(0, n); y[i]=1;
          yi = smooth.spline(x, y, df=df,spar=lam)$y;
          A[,i]= yi;
        }
        return(A)
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
        Smat = function(x, df,lam){
          n = length(x);
          A = matrix(0, n, n);
          for(i in 1:n){
            y = rep(0, n); y[i]=1;
            yi = smooth.spline(x, y, df=df,spar=lam)$y;
            A[,i]= yi;
          }
          return(A)
        } 
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.7,0.85,length.out = 20)
        for (i in 1:20){
          W    <- Smat(z,2,tp_seq[i]) 
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
        for (i2 in 1:20){
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
      
      fhat_aic <- matrix(1,n,2)
      fhat_gcv <- matrix(1,n,2)
      fhat_cp  <- matrix(1,n,2)
      
      iter <- 100
      fh1_aic <- matrix(0,n,iter)
      fh2_aic <- matrix(0,n,iter)
      
      fh1_gcv <- matrix(0,n,iter)
      fh2_gcv <- matrix(0,n,iter)
      
      fh1_cp <- matrix(0,n,iter)
      fh2_cp <- matrix(0,n,iter)
      
      
      tol    <- 0.05
      ctol   <- 99
      i <- 1
      
      while (ctol>=tol){
        for (k in 1:2){
          tp <- selectionSS(x,nc[,k],(y))
          tp.aic <- tp$lam.aic
          tp.gcv <- tp$lam.gcv
          tp.cp  <- tp$lam.cp
          S_aic  <- Smat(nc[,k],2,tp.aic) 
          S_gcv  <- Smat(nc[,k],2,tp.gcv)
          S_cp   <- Smat(nc[,k],2,tp.cp)
          
          xtil.aic <- (diag(n)-S_aic)%*%x
          xtil.gcv <- (diag(n)-S_gcv)%*%x
          xtil.cp  <- (diag(n)-S_cp)%*%x
          
          H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic)%*%t(x)%*%(diag(n)-S_aic)
          H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv)%*%t(x)%*%(diag(n)-S_gcv)
          H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp)%*%t(x)%*%(diag(n)-S_cp)
          
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
      f1<-scale(allf[,1])
      f2<-scale(allf[,2])
      f3<-scale(allf[,3])
      f4<-scale(allf[,4])
      
      n<-length(y)
      z1<-scale(nc[,1])
      z2<-scale(nc[,3])
      z3<-scale(nc[,2])
      z4<-scale(nc[,4])
      
      
      #---------------------------------------------------------------------------
      #Initialization-----------------------------------------------------------
      f01 <- fitted(lm(y~z1+x[,1]+x[,2]))
      f02 <- fitted(lm(y~z2+x[,1]+x[,2]))
      f03 <- fitted(lm(y~z3+x[,1]+x[,2]))
      f04 <- fitted(lm(y~z4+x[,1]+x[,2]))
      
      #SMOOTHING MATRIX FOR KS--------------------------------------------------
      Smat = function(x, df,lam){
        n = length(x);
        A = matrix(0, n, n);
        for(i in 1:n){
          y = rep(0, n); y[i]=1;
          yi = smooth.spline(x, y, df=df,spar=lam)$y;
          A[,i]= yi;
        }
        return(A)
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
        Smat = function(x, df,lam){
          n = length(x);
          A = matrix(0, n, n);
          for(i in 1:n){
            y = rep(0, n); y[i]=1;
            yi = smooth.spline(x, y, df=df,spar=lam)$y;
            A[,i]= yi;
          }
          return(A)
        } 
        n      <- length(y)
        index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
        tp_seq <-seq(0.7,0.85,length.out = 20)
        for (i in 1:20){
          W    <- Smat(z,2,tp_seq[i]) 
          xtil <- (diag(n)-W)%*%x
          ytil <- (diag(n)-W)%*%y
          
          beta <- solve(t(xtil)%*%xtil,tol=1e-100)%*%t(xtil)%*%ytil
          fhat <- W%*%(y-x%*%beta)
          yhat <- x%*%beta+fhat
          H    <- W+xtil%*%solve(t(xtil)%*%xtil,tol=1e-100)%*%t(x)%*%(diag(n)-W)
          aic[i] <- aiccfunc(y,yhat,H)
          gcv[i] <- gcvfunc(y,yhat,H)
          cp[i]  <- cpfunc(y,yhat,H)
        }
        for (i2 in 1:20){
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
      
      tol    <- 0.05
      ctol   <- 99
      i <- 1
      while (ctol>=tol){
        for (k in 1:4){
          tp <- selectionSS(x,nc[,k],y)
          tp.aic <- tp$lam.aic
          tp.gcv <- tp$lam.gcv
          tp.cp  <- tp$lam.cp
          S_aic  <- Smat(nc[,k],2,tp.aic) 
          S_gcv  <- Smat(nc[,k],2,tp.gcv)
          S_cp   <- Smat(nc[,k],2,tp.cp)
          
          xtil.aic <- (diag(n)-S_aic)%*%x
          xtil.gcv <- (diag(n)-S_gcv)%*%x
          xtil.cp  <- (diag(n)-S_cp)%*%x
          
          H.aic <- S_aic+xtil.aic%*%solve(t(xtil.aic)%*%xtil.aic,tol=1e-100)%*%t(x)%*%(diag(n)-S_aic)
          H.gcv <- S_gcv+xtil.gcv%*%solve(t(xtil.gcv)%*%xtil.gcv,tol=1e-100)%*%t(x)%*%(diag(n)-S_gcv)
          H.cp  <- S_cp+xtil.cp%*%solve(t(xtil.cp)%*%xtil.cp,tol=1e-100)%*%t(x)%*%(diag(n)-S_cp)
          
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
          
        }
        i <- i+1
        if (i==iter){
          break
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
