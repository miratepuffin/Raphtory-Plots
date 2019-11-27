# time series analysis of the "proportion" data
# from Ben and Naomi

library('pracma')
library('psd')


dat <- read.table('timeProportionNov15.dat')

#old method to remove the trend of the data, it works but try the pracma library
meanDat <- mean(dat$V1)
#shiftDat <- dat$V1-meanDat

# plot original data
pdf("periodicTS.pdf",height=5,width=8)
plot.ts(dat$V1);
points(rep(meanDat,length(dat$V1)),type='l',col='red',lw=2)
dev.off();

#using package pracma
shiftDat <- detrend(dat$V1,t='linear')
#pdf("peridicShift.pdf",height=5,width=8)
#plot.ts(shiftDat);
#def.off();

#using packet psd two parts am intersted $spec $freq
datTS <- ts(shiftDat)
datPrw <- pspectrum(datTS,plot=T)
resSortedFull <- sort(datPrw$spec,index.return=T,decreasing=T)


#fourier transform of the detrended data
fftDat <- fft(shiftDat)

sizeFFT <- length(fftDat)/2

# plot the modulous of the  first 4000 frequencies
tikz("fourierProportion.tex",width=3,height=3)
plot(Mod(fftDat)[0:4000],xlim=c(0,4000),type='l')
dev.off()

tmp <- datPrw$spec
tikz("PowerSPProportion.tex",width=3,height=3)
plot(tmp,xlim=c(0,4000),type='l')
dev.off()


# the auto correlation function to check if the data is stationary or not
tikz("acfProportion.tex",width=3,height=3)
acf(shiftDat,lag.max=100)
dev.off()

resSortedHalf <- sort(Mod(fftDat)[0:sizeFFT],index.return=T,decreasing=T)
resSortedFull <- sort(Mod(fftDat),index.return=T,decreasing=T)

# the hum 

ind <- seq(1,20)
tmp <- numeric(length(fftDat))
harmo <- complex(real=tmp,imaginary=tmp);
harmo[ind] <- fftDat[ind]
res1 <- fft(harmo,inverse=T)


pdf("theHum.pdf",width=8,height=5)
plot(shiftDat,type='l')
points(Re(res1)/length(res1),type='l',lw=2,col='red')
dev.off()

#first ten harmonics

ind <- resSortedFull$ix[1:10]
tmp <- numeric(length(fftDat))
harmo <- complex(real=tmp,imaginary=tmp);
harmo[ind] <- fftDat[ind]
res1 <- fft(harmo,inverse=T)


tikz("TenHarmonics.tex",width=3,height=3)
plot(shiftDat[1000:1168],type='l')
points(Re(res1)[1000:1168]/length(res1),type='l',lw=2,col='red')
dev.off()

ind <- resSortedFull$ix[1:240]
tmp <- numeric(length(fftDat))
harmo <- complex(real=tmp,imaginary=tmp);
harmo[ind] <- fftDat[ind]
res <- fft(harmo,inverse=T)
tikz("TwoHundredAndFiftyHarmonics.tex",width=3,height=3)
plot(shiftDat[1000:1168],type='l')
points(Re(res)[1000:1168]/length(res),type='l',lw=2,col='red')
dev.off()

datRemoveHar <- shiftDat - Re(res)/length(res)
#tikz("decorrelatedACF.tex",width=3,height=3)
acf(datRemoveHar,lag.max=100)
#dev.off()

#
# removing cyclcal dependency the data using multiples of 24
datToAv <- shiftDat[1:(626*24)]
mDat <- t(matrix(data=datToAv,nrow=24))


seasonalDat <- colMeans(mDat,na.rm=T)
tikz("onePeriodSeasonal.tex",width=3,height=3)
plot(seasonalDat,type='o')
dev.off();

ave <- NULL
dev <- NULL
for(i in 1:length(mDat[1,])){
  ave <- rbind(ave,mean(mDat[,i]));
  dev <- rbind(dev,sd(mDat[,i]));
}

tikz("onePeriodWithSD.tex",width=2,height=3)
plot(ave,type='l',ylim=c(-0.5,0.5))
points(ave+(dev),type='l',col='red')
points(ave-(dev),type='l',col='red')
dev.off()
#library(matrixStats)
#transform(X, SD=rowSds(X, na.rm=TRUE))


# repeat the seasonal dat 626 cylcles
trep <- rep(seasonalDat,626)

tikz("compareSeasonal.tex",width=3,height=3)
plot(shiftDat[1000:1168],type='l')
points(trep[1000:1168],type='l',lw=2,col='red')
#points(Re(res)[1000:1168]/length(res),type='l',lw=2,col='blue')
points(Re(res1)[1000:1168]/length(res1),type='l',lw=2,col='green')
dev.off()


#try to remove the seasonal of 24 hrs
deseasonDat <- datToAv-trep
tikz("acfDeseasonal.tex",width=3,height=3)
acf(deseasonDat,lag.max=100)
dev.off()


lagpad <- function(x, k) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}


# these two tex files are too big for LaTeX
#tikz("errorLag.tex",width=3,height=3)
resLag <- lagpad(datRemoveHar,4)
resCombo <- cbind(datRemoveHar,resLag)
plot(resCombo,type='p')
#dev.off();

# histogram of the error

tikz("histogramError.tex", width=3,height=3)
hist(resCombo)
dev.off();

tikz("seasonalLag4.tex",width=3,height=3)
resLag <- lagpad(Re(res)/length(res),4)
resCombo <- cbind(Re(res)/length(res),resLag)
plot(resCombo[1000:1168,],type='l')
dev.off()

#q("no")


