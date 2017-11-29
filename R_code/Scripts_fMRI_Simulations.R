


library(neuRosim)
source('Scripts_LinearModels.R')
library('fields')
library('neuRosim')

###########################
## Generate Design ##
###########################


## Event Related Design ##
generate_event_design<-function(design){
  X.binary<-matrix(0,design$TRS,length(design$conditions))
  TR=5 # wait first 5 TR's
  while(TR < design$TRS){
    X.binary[TR,sample(1:length(design$conditions),1)]=1 # pick a random condition
    TR = TR+1 # count that TR
    TR = TR + sample(design$iti[1]:design$iti[2],1) # add some TRs
  }
  X.conv<-X.binary
  for (cond in design$conditions){
    #X.conv[,cond] <- specifydesign(totaltime = design$TRS, onsets = which(X.binary[,cond]==1),durations = 1, effectsize = 1, TR = design$TRsec, conv = "Balloon")
    X.conv[,cond]<-specifydesign(totaltime = design$TRS, onsets =which(X.binary[,cond]==1), durations = 1, effectsize = 1, TR = design$TRsec, conv = "double-gamma")
  }
  design$X.binary<-X.binary
  design$X.conv<-X.conv
  
  return(design)
}


###########################
## Generate fMRI Data ##
###########################


## Data Generation ##
generate_regions_signal<-function(design,map){
  
 
  p<-dim(design$X.conv)[2]
  n<-dim(design$X.conv)[1]
  q<-map$dim1*map$dim2

  ### B  (Generate Spatial Beta Maps) (magnitudes as gaussians) 
  map$map.B<-array(0,c(q,p)) # B spatial map per condition # voxel x p
  map$map.B.region<-array(0,c(q,p,length(map$coords))) # B spatial map per condition per region # voxel x p x region
  c=1
  for (coo in map$coords){
    b <-map$which.b[[c]]
    temp <- specifyregion(dim = c(map$dim1, map$dim2), coord = coo, radius = map$radius[c],form = "sphere", fading = 0.5)
    map$map.B.region[,b,c]<-as.vector(temp)
    map$map.B[,b]<-temp+map$map.B[,b]
    c=c+1
  }
   
  # get region centers
  rc<-matrix(0,map$dim1,map$dim2)
  rr=1
  map$region.centers = c()
  for (coor in map$coords){
    rc[coor[1],coor[2]]=rr
    map$region.centers[rr]=which(rc==rr)
    rr=rr+1
  }

  ### XB (Add Temporal Signal) ### 
  map$map.ts.XB <-design$X.conv%*%t(map$map.B)
 

  # calculate sigma for SNR
  mean_signal=mean(map$map.ts.XB[map$map.ts.XB!=0])
  sigma <- mean_signal/map$SNR
 
  # E (Add in Noise)
  n.white <- systemnoise(dim = c(map$dim1,map$dim2), nscan = design$TRS, sigma = sigma, type = "rician")
  n.low <- lowfreqdrift(dim = c(map$dim1,map$dim2), nscan = design$TRS, TR = design$TRsec, freq = 120)
  n.phys <- physnoise(dim = c(map$dim1,map$dim2), nscan = design$TRS, sigma = sigma, TR = design$TRsec)
  n.temp <- temporalnoise(dim = c(map$dim1,map$dim2), sigma =sigma, nscan = design$TRS, rho = c(0.4, -0.2))
  n.spat <- spatialnoise(dim = c(map$dim1,map$dim2), sigma = sigma, nscan = design$TRS, method = "gaussRF", FWHM = 4)
  weights<-map$weights
  w <- weights
  noise <- (w[1] * n.white + w[2] * n.temp + w[3] * n.low + 
          w[4] * n.phys + w[5]  * n.spat)/sqrt(sum(w^2))
  noise<-apply(noise,3,rbind) # make 2D
  map$map.ts.E=t(noise) 
  mean_noise<-mean(map$map.ts.E[map$map.ts.E!=0])
 
  # Combine Signal + Noise
  map$map.ts.Y <- map$map.ts.XB + map$map.ts.E - mean_noise
  
 
  # MASK (Masks/Underlay (0's not Nans))
  map$underlay<-as.vector(data.matrix(read.table(map$underlayname)))
  map$mask<-as.vector(data.matrix(read.table(map$maskname)))
  ## mask should be 0's 
  
  
  # extract relevant matrices for analysis. 
  map$e.idx<-which(map$mask!=0) # that way can put back into map. 
  map$Y = map$map.ts.Y[,map$e.idx]
  map$X = design$X.conv
  map$B = map$map.B[map$e.idx,]
  map$E = map$map.ts.E[,map$e.idx]
 
  # test new map
  # Bhat <-map$B-30 
  # newmap<-matrix(0,98*98)
  # newmap[map$e.idx]<-matrix(Bhat)
  # newmap<-matrix(newmap,98,98)
  # image.plot(newmap)
  
  #### Add Nan's for Plotting.. 
  map$masknan = map$mask
  map$masknan[map$mask==0]=NA
  for (t in seq(n)){
      map$map.ts.Y[t,]<-map$map.ts.Y[t,]*map$masknan
      map$map.ts.XB[t,]<-map$map.ts.XB[t,]*map$masknan
      map$map.ts.E[t,]<-map$map.ts.E[t,]*map$masknan
  }
  for (pp in seq(p)){
      map$map.B[,pp]<-map$map.B[,pp]*map$masknan
  }
  
  return(map)
}


#v=3000
#map$fmri.data[v,]<-map$fmri.data[v,]*1000
#meandata<-matrix(apply(map$fmri.data, c(1), mean),dim1,dim2)
#image(meandata)

# count up.. 
#v=20*98-50 # should be about halfway
#m<-matrix(0,98,98)
#mv<-as.vector(m)
#m[v]=1
#image(matrix(m,98,98))

#### 3D arrays
##### library(abind)
#### test <- replicate( 5,  matrix(runif(25),ncol=5), simplify=FALSE )
#### testArray <- do.call(abind, c(test, along = 3))
#a <- specifyregion(dim = c(20, 20), coord = c(10, 10), radius = 10,form = "sphere", fading = 0.5)
#av<-as.vector(a) # turn region into vector 
#avt<-av%*%t(seq(10))
#avt<-av%*%t(ts)
#at<-matrix(avt[,2]+avt[,3],20,20)
#image(at)
#plot(avt[30,],type='l')

#map$map.ts.signal[(dim1*40)+20,]=1000
#v.idx<-matrix(seq(4*4),4,4)




###########################################
##### fMRI Simulations #####################
############################################


fmri.modelcomparisons<-function(nsims,map,design){
  

  
  ##set up storage
  est.mse<-c(0) #probably want per q too..!
  pred.mse<-c(0)
  pred.cor<-c(0)
  modellist=map$modellist
  model=map$modellist[1]
  df<-data.frame(model,est.mse,pred.mse,pred.cor,SNR=map$SNR,p=length(design$conditions),n=design$TRS)
  levels(df$model)=modellist
  
  
  maps.Bhat<-array(0,c(map$dim1,map$dim2,length(design$conditions),length(modellist),nsims))
  maps.pred.cor<-array(0,c(map$dim1,map$dim2,length(modellist),nsims))
  #maps.B<-array(0,c(map$dim1,map$dim2,dim(B)[1],length(modellist)))
  peaks.ts.Y=array(0,c(design$TRS,length(map$which.b),length(modellist)))
  peaks.ts.Yhat=array(0,c(design$TRS,length(map$which.b),length(modellist)))
  
  
  for (subject in seq(nsims)){
    
  
  
  
  ##### Collect Training Data
  design.train<-generate_event_design(design)
  map.train<-generate_regions_signal(design.train,map)
  
  #### Collect Testing Data
  #### maybe change up study parameters.. ### 
  design.test<-generate_event_design(design)
  map.test<-generate_regions_signal(design.test,map)
  
  # get distance matrix for search light
  distrank<-data.matrix(read.table("sim_distrank_masked.txt"))
  #distrank<-distrank[map.train$e.idx,]
  # change the distrank, so the indexes reflect mask

  
  # extract data
  X.train<-map.train$X # nxq
  Y.train<-map.train$Y # needs to be nxq
  B<-t(map.train$B)
  X.test<-map.test$X
  Y.test<-map.test$Y # needs to be nxq
  
  p<-length(design$conditions)
    
  
  # loop through models
  mm=1
  for (mname in modellist){
    
    funname = paste(mname,".fit",sep='')
    
    # Fit 
    if (grepl('cw',mname)){
      # Search Light
      fitted<-searchlight.fit(funname,X.train,Y.train,distrank,map$search.size)
    }else{
      # fit on all voxels
      fitted<-do.call(funname,list(X.train=X.train,Y.train=Y.train))
    }
    

    # estimation error
    results.e<-estimation.error(B,fitted$B.hat)
    
    # grab Beta maps.. for each model
    for (bb in seq(p)){
      emptymap<-matrix(NA,map$dim1*map$dim2)
      emptymap[map.train$e.idx]<-fitted$B.hat[bb,]
      maps.Bhat[,,bb,mm,subject]=matrix(emptymap,map$dim1,map$dim2)
    }
    #predict test set
    Yhat.test<-predict.me(fitted,X.test)
    
    #prediction error
    results.p<-prediction.error(Y.test,Yhat.test)

    emptymap<-matrix(NA,map$dim1*map$dim2)
    emptymap[map.train$e.idx]<-results.p$map.cor
    maps.pred.cor[,,mm,subject]=emptymap
    
   # for (rr in map$which.b){
  #    peaks.ts.Y[,rr,mm]=Y.test[,map$region.centers[rr]]
  #    peaks.ts.Yhat[,rr,mm]=Yhat.test[,map$region.centers[rr]]
  #  } 
    
    ### Store performance ####
    newrow = c(1:dim(df)[2])
    newrow[1]=mname
    newrow[2]=results.e$est.mse2
    newrow[3]=results.p$pred.mse
    newrow[4]=results.p$pred.cor
    newrow[5]=map$SNR
    newrow[6]=p
    newrow[7]=design$TRS
    df = rbind(df,newrow)
    mm=mm+1
    
  }

  
  }
  
  ## make numeric
  df<-df[2:dim(df)[1],]
  for (dd in seq(2,dim(df)[2])){
    df[,dd]<-as.numeric(df[,dd])
  }
  results=list(df=df,maps.Bhat=maps.Bhat,peaks.ts.Y=peaks.ts.Y,peaks.ts.Yhat=peaks.ts.Yhat)
  return(results)
}






searchlight.fit<-function(funname,X.train,Y.train,distrank,size){
  # if don't have voxelidx.. just use seq(1:end)
  Yhat.train<-matrix(0,dim(Y.train)[1],dim(Y.train)[2]) # nxq
  B.hat<-matrix(0,dim(X.train)[2],dim(Y.train)[2]) # pxq

  for (voxel in seq(dim(Y.train)[2])){
    Y.train.search = Y.train[,distrank[voxel,1:size]] # n x searchlight size ## 0's in the index get skipped cool.. THat's a more useful feature
    if(is.null(dim(Y.train.search))){
      # only one voxel in search
      Y.train.search=matrix(Y.train.search)
    }
    fitted.search<-do.call(funname,list(X.train=X.train,Y.train=Y.train.search)) 
    B.hat[,voxel]<-fitted.search$B.hat[,1] # should be the first voxel, but take the index in distrank == voxel number 
    Yhat.train[,voxel]<-fitted.search$Yhat.train[,1] 
  }
  fitted<-list(B.hat=B.hat,Yhat.train=Yhat.train,succeed=1)
  return(fitted)
}


calc_dist<-function(dim1,dim2,topn=50){
  #for 2D.. 
  #
  distrank=matrix(0,dim1*dim2,topn) # each voxel give its closest neighbors (in flattened coordinates)
  v = 1
  for (x in seq(dim1)){
    for (y in seq(dim2)){ # go down first
      one.voxel.dist<-matrix(0,dim1*dim2)
      ov=1
      for (xx in seq(dim1)){
        for (yy in seq(dim2)){ 
          one.voxel.dist[ov]=sqrt((x-xx)**2+(y-yy)**2)
          ov=ov+1
        }
      }
      distrank[v,] = sort(one.voxel.dist,index.return=TRUE)[1:topn][[2]][1:topn]
      v=v+1 # this is fine and in order.. 
    }
    print(x)
  }
  
  ## saved 
  
  
  ### Take the full 98x98 distance rankings and transform their indexes into this smaller s
  ### there's defintely a formula for this
  distrank.new<-distrank[map.test$e.idx,]

  for (dd in seq(dim(distrank.new)[1])){
    for (ddd in seq(dim(distrank.new)[2])){
      #print('here')
      em = matrix(0,map.train$dim1,map.train$dim2)
      em<-as.vector(em)
      em[distrank.new[dd,ddd]]=1
      em<-em[map.train$e.idx]
      if((sum(em)==1)==1){
        distrank.new[dd,ddd]<-which(em==1) # its position in the masked matrix
      }else{
        distrank.new[dd,ddd]<-0
      }
    }

  }
  write.table(distrank.new,file="sim_distrank_masked.txt")
  
  return(distrank)
}



