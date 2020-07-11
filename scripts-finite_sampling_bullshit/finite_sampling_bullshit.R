## Calculations of long-run and sample mutual info ##

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('LaplacesDemon')
#library('RNetCDF')
#library('Rmpfr')
options(bitmapType='cairo')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()

set.seed(181225)



## freqs[,S] = response freqs for stimulus S
## assumes all stimuli equally probable
mutualinfo <- function(freqs,base=2){##in bits
    stimfreqs <- 1/ncol(freqs)
    jointfreqs <- freqs*stimfreqs
    respfreqs <- rowSums(jointfreqs)

    sum(jointfreqs*log2(jointfreqs/outer(respfreqs,rep(stimfreqs,ncol(freqs)))),na.rm=TRUE)/log2(base)
}

freqsampling <- function(freqs,nsamples=20){
    nresps <- nrow(freqs)
    foreach(s=1:ncol(freqs),.combine=cbind)%do%{
        onesample <- sample(1:nresps,nsamples,replace=TRUE,prob=freqs[,s])
        tabulate(onesample,nbins=nresps)
    }    
}


plotsingle <- function(rfreqs0,filename,tit0,nmcsamples=1000,base=2){
    nstim <- ncol(rfreqs0)
    milongrun0 <- mutualinfo(rfreqs0)

    misamples0 <- foreach(i=1:nmcsamples,.combine=c)%do%{
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            rdirichlet(n=1,alpha=rfreqs0[,s])
        }
        mutualinfo(fsample)
    }

    hires <- hist(misamples0,breaks='FD',plot=FALSE)
    maxy <- max(hires$counts)
    maxx <- max(hires$breaks,milongrun0*1.1)
    pdff(paste0(filename))
    matplot(x=misamples0,y=rep(-maxy/20,nmcsamples),type='p',lty=1,lwd=3,pch='.',cex=2,col=myblue,xlim=c(0,maxx),ylim=c(-maxy/20,maxy),xlab='sampled I/bit',ylab='',main=tit0)
    hist(misamples0,breaks='FD',xlim=c(0,maxx),ylim=c(-maxy/20,maxy),xlab='I',ylab='',add=TRUE)
##legend('top',paste0('\nblack: mean, blue: median\nyellow: 16%q'),bty='n')
    matlines(x=rep(milongrun0,2),y=c(0,maxy),type='l',lty=1,lwd=3,pch='.',col=myred)
dev.off()
}

freqs0 <- matrix(1/10,10,2)
plotsingle(freqs0,'hist_caseA', 'case A')


f1 <- c(rep(1/7,7),rep(0,3))
freqs0 <- cbind(f1,rev(f1))
plotsingle(freqs0,'hist_caseB', 'case B')

