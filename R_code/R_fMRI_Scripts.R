
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}



########################################################################


layout(matrix(c(1,2,3,0,4,0), nrow=2, ncol=3), widths=c(4,4,1), heights=c(4,1))
layout.show(4)

pal.1=colorRampPalette(c("blue","black", "red", "yellow"), space="rgb")
breaks <- seq(min(n.gaus.image), max(n.gaus.image),length.out=100)
par(mar=c(1,1,1,1))
image(n.gaus.image, col=pal.1(length(breaks)-1), breaks=breaks, xaxt="n", yaxt="n", ylab="", xlab="")
points(4,4, pch=2, lwd=2, cex=2,col="blue")
par(mar=c(3,1,1,1))
image.scale(n.gaus.image, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)
box()



########### HIGHER LEVEL FUNCTIONS ###########


TR <- 2
nscan <- 100
total <- TR * nscan
os1 <- seq(1, total, 40)
os2 <- seq(15, total, 40)
dur <- list(20, 7)
os <- list(os1, os2)
effect <- list(3, 10)

design <- simprepTemporal(totaltime = total, onsets = os, durations = dur,
                          effectsize = effect, TR = TR, hrf = "double-gamma")

par(mfrow=c(1,1))
w <- c(0.3, 0.3, 0.01, 0.09, 0.3)
ts <- simTSfmri(design = design, base = 10,
                SNR = 2, noise = "mixture", type = "rician", weights = w, verbose = FALSE)
plot(ts,type='l')
plot(ts_signal,type='l')

# mixture contains Rician system noise, temporal noise of order 1, low-frequency
# drift, physiological noise and task-related noise and has a baseline value of 10



#### generating volumes of data ###

regions <- simprepSpatial(regions = 2, coord = list(c(10, 5, 24),
                          c(53, 29, 24)), radius = c(10, 5), form = "sphere")



onset <- list(os, os)
duration <- list(dur, dur)
effect1 <- list(2, 9)
effect2 <- list(14, 8)
design2 <- simprepTemporal(regions = 2, onsets = onset,
                          durations = duration, TR = TR, hrf = "double-gamma",
                          effectsize = list(effect1, effect2), totaltime = total)

w <- c(0.3, 0.3, 0.01, 0.09, 0.1, 0.2)
data <- simVOLfmri(dim = c(64, 64, 64), base = 100, design = design2,
                      image = regions, SNR = 10, noise = "mixture", type = "rician",
                      weights = w, verbose = FALSE)

#Note that with simTSfmri and simVOLfmri it is also
#possible to simulate data that contain only activation or only noise




library("oro.nifti")
epi <- readNIfTI("EPI.nii")
baseline <- epi
baseline.bin <- ifelse(baseline > .4, 1, 0)
ix <- which(baseline == 1)
baseline[-ix] <- 0
image(baseline[,,30])
image(baseline.bin[,,30])



n.gaus <- spatialnoise(dim = c(91,109,91), sigma = 15, nscan = 100, method = "gaussRF", FWHM = 4)
n.gaus.image <-apply(n.gaus, c(1,2,3), mean)

n.gaus.image.masked = n.gaus.image*baseline.bin
#random <- nifti(array(runif(91*109*91), c(91,109,91)))
baseline@.Data = n.gaus.image.masked
writeNIfTI(baseline, "randomfield")


T1 <- readNIfTI("T1.nii")
writeNIfTI(baseline, "T11")



#
URL <- "http://imaging.mrc-cbu.cam.ac.uk/downloads/Colin/colin_1mm.tgz"
urlfile <- "colin_1mm.tgz"
if (!file.exists(urlfile)){
  download.file(URL, dest=urlfile, quiet=TRUE)
  untar(urlfile)
}

img <- readNIfTI("colin_1mm")
img = cal_img(img)
writeNIfTI(img, "T11")
tfile = tempfile()
writeNIfTI(img, filename = tfile, verbose=TRUE)



Z = matrix(.9,2,2)
diag(Z)=1
out<-mvrnorm(100,matrix(0,2),Z)
y1<-out[,1]
x<-out[,2]
#plot(x,y1)
y2 <-rnorm(100,0,1)
print(cor(x,y1))

Y = cbind(y1,y2)
cc<-cancor(x,Y)
print(cc)

plot3d(y1,y2,x)

Yn = t(cc$ycoef%*%t(Y))
y1n=Yn[,1]
y2n=Yn[,2]
#print(Yn)
#scatterplot3d(y1n,y2n,x)
plot3d(y1n,y2n,x)




#### Their Example ### 

TR <- 2
nscan <- 100
total <- TR * nscan
os1 <- seq(1, total, 40)
os2 <- seq(15, total, 40)
dur <- list(20, 7)
os <- list(os1, os2)
effect <- list(3, 10)

regions <- simprepSpatial(regions = 2, coord = list(c(10, 5, 24),
                                                    c(53, 29, 24)), radius = c(10, 5), form = "sphere")

onset <- list(os, os)
duration <- list(dur, dur)
effect1 <- list(2, 9)
effect2 <- list(14, 8)
design2 <- simprepTemporal(regions = 2, onsets = onset,
                              durations = duration, TR = TR, hrf = "double-gamma",
                              effectsize = list(effect1, effect2), totaltime = total)

w <- c(0.3, 0.3, 0.01, 0.09, 0.1, 0.2)
data <- simVOLfmri(dim = c(64, 64, 64), base = 100, design = design2,
                      image = regions, SNR = 10, noise = "mixture", type = "rician",
                      weights = w, verbose = FALSE)

