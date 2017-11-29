
library(MASS)
library(car)
library('ggplot2')
library('pracma')
library('glmnet')

############################################


generate_data<-function(params){
  
   n<-params$n
   p<-params$p
   q<-params$q
   bsig<-params$bsig
   SNR<-params$SNR
   
  
  #### Population X (Signal) Covariance
  V= matrix(0,p,p)
  for (i in seq(p)){
    for (j in seq(p)){
      rr=round(runif(1, min=-1, max=1),1)
      V[i,j]=rr}}
  V = t(V)%*%V # makes it postive-semi definite
  X = round(mvrnorm(n,mu=matrix(0,p),V),4)
  
  ##Regenerate B's.. 
  
  # want correlated across y but independent across x.. 
  # so draw each predict seperately from a mv guassian.. # .. figure out later...
  B = matrix(0,p,q)
  for (pp in seq(p)){
    # each predictor x, has some effect.. random -1 to 1.
    # then y have similar responses depending on the spread.
    B[pp,]<-rnorm(q,runif(1,min=-1,max=1),bsig)
  }
  # normalize across q (not across p)
  for (qq in seq(q)){
    B[,qq]=(B[,qq]-mean(B[,qq]))/sd(B[,qq])
  }

  
  #B = matrix(round(mvrnorm(p*q,rep(shared_response,p*q),diag(p*q)*bsig),4),p,q)
  
  # signal to noise matrix
  F = t(B)%*%V%*%B
  
  #Calculate appropriate Noise
  sigma=mean(F)/SNR
  # R = Sigma%*%solve(F) # noise-signal matrix.. 1/R is signal to noise... geeks out if q>p..
  
  #### Population E
  Sigma = diag(q)*sigma
  E = round(mvrnorm(n,matrix(0,q),Sigma),4) # could also do diagnol
  
  # Responses 
  Y = X%*%B+ E
  

  X.train = X[1:(n/2),]
  Y.train = Y[1:(n/2),]  
  X.test = X[(n/2+1):n,]
  Y.test = Y[(n/2+1):n,]

  # Center  (do i need to?)
  for (y in seq(q)){
    Y.train[,y]=Y.train[,y]-mean(Y.train[,y])
    Y.test[,y]=Y.test[,y]-mean(Y.test[,y])
  }
  for (x in seq(p)){
    X.train[,x]=X.train[,x]-mean(X.train[,x])
    X.test[,x]=X.test[,x]-mean(X.test[,x])
  }

  data = list(X.train=X.train,Y.train=Y.train,Y.test=Y.test,X.test=X.test,B=B,V=V)
}




########################################################



ols.fit<-function(X.train,Y.train){
  B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  Yhat.train<- X.train%*%B.hat
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=1)
  return(fitted)
}



ols.cw.gcv.fit<-function(X.train,Y.train){
  
  
  B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  r = p/n
 
  tryCatch({
    Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
    e<-eigen(Q)
    That=e$vectors
    ## alternative (exact same)
    #cc=cancor(X.train,Y.train)#xcenter=FALSE,ycenter=FALSE
    #That=cc$ycoef ### I think X>Y.. here =
    #C2=diag(q)*(cc$cor**2) # squared
    ## shrink the canonical correlations based on generalized cross validation
    C2=diag(q)*e$values # the eigen values are the squared canonical correlations. 
    #if (p>q){
    #  diag(C2)[(p+1):q]=0 # oterwise it repeats ccors 
    #}
    d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
    Dhat = diag(q)*diag(d)
    Dhat[Dhat<0]=0 # take positive part as B-F recomment
    S = That%*%Dhat%*%solve(That) # project into canonical directions, and project back..
    succeed=1
    B.hat<-B.hat%*%S # shrink the B.hat estimates, though not sure if these are interpretable.. 
    Yhat.train<- X.train%*%B.hat 
    fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
    return(fitted)
  },error= function(c){
    Q=ginv(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%ginv(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
    e<-eigen(Q)
    That=Re(e$vectors) # I have no idea if this is ok
    C2=diag(q)*Re(e$values)
    #if (p>q){
    #  diag(C2)[(p+1):q]=0 # oterwise it repeats ccors 
    #}
    d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
    Dhat = diag(q)*diag(d)
    Dhat[Dhat<0]=0
    S = That%*%Dhat%*%ginv(That)
    succeed=0
    B.hat<-B.hat%*%S # shrink the B.hat estimates, though not sure if these are interpretable.. 
    Yhat.train<- X.train%*%B.hat 
    fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
    return(fitted)
  })
  
}


ols.cw.cv.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-50
  #d.options<-linspace(0,.999,options)
  #cverr.fold.d.q<-array(0,c(folds,options,q))

  #
  
  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain)%*%t(cv.Xtrain)%*%cv.Ytrain
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Xtrain)%*%solve(t(cv.Xtrain)%*%cv.Xtrain)%*%(t(cv.Xtrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)
    
    C2=diag(q)*Re(e$values)
    r=p/n
    d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
    Dhat.c = diag(q)*diag(d)
    

    # oh right each y is shrunk seperately so its fine.. to iterate through.. 
    for (qq in seq(q)){
      dd=1
      for (d in d.options){
        #Dhat.try<-diag(q) # make all others 1 # when he says quadratic.. 
        Dhat.try<-Dhat.c # try gcv as starting point for other's.. 
        Dhat.try[qq,qq]<- d # Dhat.try[qq,qq]-d # make this basis vector shrunk.. 
        S <- That.cv%*%Dhat.try%*%solvedThat.cv
        cverr.fold.d.q[i,dd,qq]<-mean((cv.Ytest-cv.Xtest%*%B.hat.cv.d%*%S)**2)
        dd=dd+1
      }
      
    }
  }
  
  # find d that gives minimum cv err
  C2=diag(q)*Re(e$values)
  r=p/n
  d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
  Dhat = diag(q)*diag(d)
  
  # subtract off best for each
  cverr.fold.d.mean<-apply(cverr.fold.d.q,c(2,3),mean)
  for (qq in seq(q)){
    Dhat[qq,qq]<-d.options[which.min(cverr.fold.d.mean[,qq])] #Dhat[qq,qq]-d.options[which.min(cverr.fold.d.mean[,qq])]
  }
  

  
  # refind basis with full set
  Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
  e<-eigen(Q)
  That=e$vectors

  S = That%*%Dhat%*%solve(That) 
  
  # get beta
  B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  B.hat<-B.hat%*%S 
  
  # get Yhat 
  Yhat.train<- X.train%*%B.hat 
  
  succeed=1
  
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed,d=d)
  return(fitted)
  
  
}


ols.cw.cv2.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-100
  d.options<-linspace(-1,1,options)
  s.options<-linspace(0.01,10,options)
  cverr.fold.d.q<-array(0,c(folds,options,options))
  
  
  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain)%*%t(cv.Xtrain)%*%cv.Ytrain
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Xtrain)%*%solve(t(cv.Xtrain)%*%cv.Xtrain)%*%(t(cv.Xtrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)
    C2=diag(q)*Re(e$values)
    r=p/n
    d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
    Dhat.c = diag(q)*diag(d)
  #  Dhat.c = diag(q)*diag(C2)  
  
    # oh right each y is shrunk seperately so its fine.. to iterate through.. 
    sss=1
    for (ss in s.options){
      dd=1
      for (d in d.options){
        Dhat.try<-Dhat.c-diag(q)*d
        Dhat.try[Dhat.try<0]=0
        S <- That.cv%*%Dhat.try%*%solvedThat.cv
        cverr.fold.d.q[i,dd,sss]<-mean((cv.Ytest-cv.Xtest%*%B.hat.cv.d%*%S)**2)
        dd=dd+1
      }
    sss=sss+1
    }
  }
  
  # find d that gives minimum cv err

  cverr.fold.d.mean<-apply(cverr.fold.d.q,c(2,3),mean)
  best<-which(cverr.fold.d.mean == min(cverr.fold.d.mean), arr.ind = TRUE)
  sbest<-s.options[best[2]]
  dbest<-d.options[best[1]]
  
  # refind basis with full set
  Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
  e<-eigen(Q)
  That=e$vectors
  C2=diag(q)*Re(e$values)
  r=p/n
  d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
  Dhat.c = diag(q)*diag(d)
  #Dhat.c = diag(q)*diag(C2)  
  Dhat.c = Dhat.c -diag(q)*dbest
  Dhat.c[Dhat.c<0]=0
  
  S = That%*%Dhat.c%*%solve(That) 
  
  # get beta
  B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  B.hat<-B.hat%*%S 
  
  # get Yhat 
  Yhat.train<- X.train%*%B.hat 
  
  succeed=1
  
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed,d=d)
  return(fitted)
  
  
}




ols.cw.cv3.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-50
  d.reg<-matrix(0,folds,q)

  
  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain)%*%t(cv.Xtrain)%*%cv.Ytrain
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Xtrain)%*%solve(t(cv.Xtrain)%*%cv.Xtrain)%*%(t(cv.Xtrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)

    Dhat.cv<-diag(q)
    YT<-cv.Ytest%*%That.cv
    YhatT<-cv.Xtest%*%B.hat.cv%*%That.cv
    for (qq in seq(q)){
      r<-lm(YT[,qq]~YhatT[,qq])
      d.reg[i,qq]<-r$coef[2]
    }
    
  }

d.reg.mean<-apply(d.reg,2,mean)
  

# refind basis with full set
Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
e<-eigen(Q)
That=e$vectors

Dhat<-diag(q)*diag(d.reg.mean)

Dhat[Dhat<0]=0

S = That%*%Dhat%*%solve(That) 

# get beta
B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
B.hat<-B.hat%*%S 

# get Yhat 
Yhat.train<- X.train%*%B.hat 

succeed=1

fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
return(fitted)


}




ols.cw.cv4.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-50

  
  YT.all<-matrix(0,n,q)
  YhatT.all<-matrix(0,n,q)
  

  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain)%*%t(cv.Xtrain)%*%cv.Ytrain
    
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Xtrain)%*%solve(t(cv.Xtrain)%*%cv.Xtrain)%*%(t(cv.Xtrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)
    
    Dhat.cv<-diag(q)
    YT<-cv.Ytest%*%That.cv
    YhatT<-cv.Xtest%*%B.hat.cv%*%That.cv
    
    YT.all[idx.tf,]<-YT
    YhatT.all[idx.tf,]<-YhatT
  }
  
  d.reg<-matrix(0,q)
  Dhat<-diag(q)
  for (qq in seq(q)){
    r<-lm(YT.all[,qq]~YhatT.all[,qq])
    d.reg[qq]<-r$coef[2]
    Dhat[qq,qq]<-d.reg[qq]
  }
  
  # refind basis with full set
  Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
  e<-eigen(Q)
  That=e$vectors

  
  Dhat[Dhat<0]=0
  
  S = That%*%Dhat%*%solve(That) 
  
  # get beta
  B.hat<- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  B.hat<-B.hat%*%S 
  
  # get Yhat 
  Yhat.train<- X.train%*%B.hat 
  
  succeed=1
  
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
  return(fitted)
  
  
}



ridge.cw.cv4.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-50
  
  
  YT.all<-matrix(0,n,q)
  YhatT.all<-matrix(0,n,q)
  
  lambdas<-matrix(0,q)
  for (qq in seq(dim(Y.train)[2])){
    cvfit = cv.glmnet(X.train, Y.train[,qq],nfolds=5,nlambda=20)
    lambdas[qq]<-cvfit$lambda.min

  }
  
  lambda<-mean(lambdas)
                
  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain+diag(p)*lambda)%*%t(cv.Xtrain)%*%cv.Ytrain
    
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Xtrain)%*%solve(t(cv.Xtrain)%*%cv.Xtrain)%*%(t(cv.Xtrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)
    
    Dhat.cv<-diag(q)
    YT<-cv.Ytest%*%That.cv
    YhatT<-cv.Xtest%*%B.hat.cv%*%That.cv
    
    YT.all[idx.tf,]<-YT
    YhatT.all[idx.tf,]<-YhatT
  }
  
  d.reg<-matrix(0,q)
  Dhat<-diag(q)
  for (qq in seq(q)){
    r<-lm(YT.all[,qq]~YhatT.all[,qq])
    d.reg[qq]<-r$coef[2]
    Dhat[qq,qq]<-d.reg[qq]
  }
  
  # refind basis with full set
  Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%X.train)%*%solve(t(X.train)%*%X.train)%*%(t(X.train)%*%Y.train)
  e<-eigen(Q)
  That=e$vectors
  
  
  Dhat[Dhat<0]=0
  
  S = That%*%Dhat%*%solve(That) 
  
  # get beta
  B.hat<- solve(t(X.train)%*%X.train+diag(p)*lambda)%*%t(X.train)%*%Y.train
  B.hat<-B.hat%*%S 
  
  # get Yhat 
  Yhat.train<- X.train%*%B.hat 
  
  succeed=1
  
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
  return(fitted)
  
                
}



ridge.cw.cv5.fit<-function(X.train,Y.train){
  
  
  p=dim(X.train)[2]
  n=dim(X.train)[1]
  q=dim(Y.train)[2]
  
  folds<-10
  options<-50
  
  
  YT.all<-matrix(0,n,q)
  YhatT.all<-matrix(0,n,q)
  
  lambdas<-matrix(0,q)
  for (qq in seq(dim(Y.train)[2])){
    cvfit = cv.glmnet(X.train, Y.train[,qq],nfolds=5,nlambda=20)
    lambdas[qq]<-cvfit$lambda.min
    
  }
  
  lambda<-mean(lambdas)
  
  idx<-seq(1,n,folds)
  for (i in seq(folds)){
    
    # get index # 
    if (i==folds){ # last fold
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:n]<-TRUE
    }else{
      idx.tf<-rep(FALSE,n)
      idx.tf[idx[i]:(idx[i+1]-1)]<-TRUE
    }
    
    # split up data #
    cv.Ytrain<-Y.train[!idx.tf,]
    cv.Ytest<-Y.train[idx.tf,]
    cv.Xtrain<-X.train[!idx.tf,]
    cv.Xtest<-X.train[idx.tf,]
    
    # for each fold (calculate beta, and T^)
    B.hat.cv<- solve(t(cv.Xtrain)%*%cv.Xtrain+diag(p)*lambda)%*%t(cv.Xtrain)%*%cv.Ytrain
    cv.Yhattrain<-cv.Xtrain%*%B.hat.cv
    
    # calculate canonical correlations with ridge #
    Q.cv=solve(t(cv.Ytrain)%*%cv.Ytrain)%*%(t(cv.Ytrain)%*%cv.Yhattrain)%*%solve(t(cv.Yhattrain)%*%cv.Yhattrain)%*%(t(cv.Yhattrain)%*%cv.Ytrain)
    e<-eigen(Q.cv)
    That.cv=e$vectors
    solvedThat.cv<-solve(That.cv)
    
    Dhat.cv<-diag(q)
    YT<-cv.Ytest%*%That.cv
    YhatT<-cv.Xtest%*%B.hat.cv%*%That.cv
    
    YT.all[idx.tf,]<-YT
    YhatT.all[idx.tf,]<-YhatT
  }
  
  d.reg<-matrix(0,q)
  Dhat<-diag(q)
  for (qq in seq(q)){
    r<-lm(YT.all[,qq]~YhatT.all[,qq])
    d.reg[qq]<-r$coef[2]
    Dhat[qq,qq]<-d.reg[qq]
  }
  
  # get beta
  B.hat<- solve(t(X.train)%*%X.train+diag(p)*lambda)%*%t(X.train)%*%Y.train
  Yhat.train<-X.train%*%B.hat
  
  
  # refind basis with full set
  Q=solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%Yhat.train)%*%solve(t(Yhat.train)%*%Yhat.train)%*%(t(Yhat.train)%*%Y.train)
  e<-eigen(Q)
  That=e$vectors
  
  
  Dhat[Dhat<0]=0
  
  S = That%*%Dhat%*%solve(That) 
  
 # shrink beta
  B.hat<-B.hat%*%S 
  
  # get Yhat 
  Yhat.train<- X.train%*%B.hat 
  
  succeed=1
  
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
  return(fitted)
  
  
}





ridge.fit<-function(X.train,Y.train){
  
  ## Ridge ##
  q = dim(Y.train)[2]
  p= dim(X.train)[2]
  n=dim(X.train)[1]
  Yhat.train.rr<-matrix(0,n,q)
  B.hat<-matrix(0,dim(X.train)[2],q)
  for (qq in seq(dim(Y.train)[2])){
    cvfit = cv.glmnet(X.train, Y.train[,qq],nfolds=5,nlambda=20)
    #print(qq)
    B.hat[,qq]<-matrix(coef(cvfit,s=cvfit$lambda.min)[2:(p+1),1])
    #Yhat.test.rr[,qq]<-predict(cvfit, newx = X.test, s = "lambda.min")
    Yhat.train.rr[,qq]<-predict(cvfit, newx = X.train, s = "lambda.min")
  }
  ####
  Yhat.train<-Yhat.train.rr
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=1)
  return(fitted)
}


ridge.cw.gcv.fit<-function(X.train,Y.train){
  ## Ridge ##
  q = dim(Y.train)[2]
  p= dim(X.train)[2]
  n=dim(X.train)[1]
  Yhat.train.rr<-matrix(0,n,q)
  B.hat<-matrix(0,dim(X.train)[2],q)
  for (qq in seq(dim(Y.train)[2])){
    cvfit = cv.glmnet(X.train, Y.train[,qq],nfolds=5,nlambda=20)
    B.hat[,qq]<-matrix(coef(cvfit,s=cvfit$lambda.min)[2:(p+1),1])
    Yhat.train.rr[,qq]<-predict(cvfit, newx = X.train, s = "lambda.min")
  }
  ####
  Yhat.train<-Yhat.train.rr
  
  #####
  r = p/n ### ADJUST FOR RIDGE ##### 
  # not clear if they do a seperate r for each q?

  #### Canonical Correlations on XB_ridge, rather than X
  tryCatch({
      Q<-solve(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%Yhat.train)%*%solve(t(Yhat.train)%*%Yhat.train)%*%(t(Yhat.train)%*%Y.train)
      e<-eigen(Q)
      That=e$vectors
      ## shrink the canonical correlations based on generalized cross validation
      C2=diag(q)*e$values # the eigen values are the squared canonical correlations.
      #if (p>q){
      #  diag(C2)[(p+1):q]=0 # oterwise it repeats ccors 
      #}
      d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
      Dhat = diag(q)*diag(d)
      Dhat[Dhat<0]=0 # take positive part as B-F recomment
      S = That%*%Dhat%*%solve(That) # project into canonical directions, and project back.. 
      succeed=1
      B.hat<-B.hat%*%S
      fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
      return(fitted)
  },error= function(c){

      Q=ginv(t(Y.train)%*%Y.train)%*%(t(Y.train)%*%Yhat.train)%*%ginv(t(Yhat.train)%*%Yhat.train)%*%(t(Yhat.train)%*%Y.train)
      e<-eigen(Q)
      That=Re(e$vectors) # I have no idea if this is ok
      C2=diag(q)*Re(e$values)
        #if (p>q){
        #  diag(C2)[(p+1):q]=0 # oterwise it repeats ccors 
        #}
      d=(1-r)*(C2-r) / (((1-r)^2*C2) + r^2*(1-C2)) 
      Dhat = diag(q)*diag(d)
      Dhat[Dhat<0]=0
      S = That%*%Dhat%*%ginv(That)
      succeed=0
      B.hat<-B.hat%*%S
      fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=succeed)
      return(fitted)
  })
  
}







###########################################
##### B-H like Simulations #########
############################################

simulation<-function(params,nsims){
  
  #set up storage
  est.mse<-c(0) #probably want per q too..!
  pred.mse<-c(0)
  pred.cor<-c(0)
  est.mse2<-c(0)
  succeed<-c(0)
  params$modellist
  model=modellist[1]
  d<-1
  df<-data.frame(model,params$n,params$p,params$q,params$bsig,params$SNR,est.mse,pred.mse,pred.cor,est.mse2,succeed,d)
  levels(df$model)=modellist
  
  ### Start Simulation
  simulations = seq(nsims)
  
  # mse=matrix(0,nsims,3) # to compare
  
  for (s in simulations){

    
    #generate data
    data<-generate_data(params)
    
    # loop through models
    mm=1
    for (mname in modellist){
      #fit 
      funname = paste(mname,".fit",sep='')
      fitted<-do.call(funname,list(X.train=data$X.train,Y.train=data$Y.train))
      d = 1
      if (mname == "ols.cw.cv"){
        d<-fitted$d
      }
      
      # estimation error
      results.e<-estimation.error(data$B,fitted$B.hat,data$V)
      
      #predict test set
      Yhat.test<-predict.me(fitted,data$X.test)
      
      #prediction error
      results.p<-prediction.error(data$Y.test,Yhat.test)
      #print(results.p$pred.mse)
      #mse[s,mm]<-mean((Yhat.test-data$Y.test)**2)
      
      #stack the dataframe 
      newrow = c(1:dim(df)[2])
      newrow[1]=mname
      newrow[2]=params$n
      newrow[3]=params$p
      newrow[4]=params$q
      newrow[5]=params$bsig
      newrow[6]=params$SNR
      newrow[7]=results.e$est.mse
      newrow[8]=results.p$pred.mse
      newrow[9]=results.p$pred.cor
      newrow[10]=results.e$est.mse2
      newrow[11]=fitted$succeed
      newrow[12]=d
      df = rbind(df,newrow)
      mm=mm+1
    }
  }

  df<-df[2:dim(df)[1],]
  for (dd in seq(2,dim(df)[2])){
    df[,dd]<-as.numeric(df[,dd])
  }
  
  return(df)
  
}


estimation.error<-function(B,B.hat,V=NULL){
  p=dim(B)[1]
  q=dim(B)[2]
  se<-matrix(0,q)
  if (is.null(V)){
    V=diag(p)
  }
  for (qq in seq(q)){
    se[qq]<-t(B[,qq])%*%V%*%B.hat[,qq]
  }
  # criteria 5.12 basically. 
  return(list(est.mse=sum(se)/(p*q),est.mse2=mean((B-B.hat)**2)))
}

prediction.error<-function(Y.test,Yhat.test){
  n=dim(Y.test)[1]
  q=dim(Y.test)[2]
  pred.mse<-mean((Yhat.test-Y.test)**2)
  #pred.cor<-mean(diag(cor(Y.test,Yhat.test)))#correlation of each yi with corresponding yhati, across n
  cors<-matrix(0,q)
  for (qq in seq(q)){
    cors[qq]<-cor(Y.test[,qq],Yhat.test[,qq])
  }
  
  mse<-matrix(0,q)
  for (qq in seq(q)){
    mse[qq]<-mean((Y.test[,qq]-Yhat.test[,qq])**2)
  }
  
  results<-list(pred.mse=mean(mse),pred.cor=mean(cors,na.rm=TRUE),map.cor=cors)
  return(results)
}



predict.me<-function(fitted,X.test){
  Yhat.test<-X.test%*%fitted$B.hat
  return(Yhat.test)
}







library(gridExtra)
library(ggplot2)

modellist<-c("ols","ols.cw.gcv","ridge.cw.gcv","ols.cw.cv4","ridge.cw.cv4","ridge.cw.cv5")
params=list(n=200,p=50,q=25,bsig=.1,SNR=1)
nsims<-5
df<-simulation(params,nsims)

df1 <- subset(df, params.p < params.q )
p1<-ggplot(data=df,aes(x=model,y=est.mse))+geom_boxplot()
p2<-ggplot(data=df,aes(x=model,y=pred.mse))+geom_boxplot()
p3<-ggplot(data=df,aes(x=model,y=pred.cor))+geom_boxplot()
grid.arrange(p1,p2,p3, nrow = 2, ncol = 2)


