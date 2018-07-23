## Application of the two Ising parameter priors to neurons from Nicola's data
## cross-validation-like analysis

library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
# library('plot3D')
# library('doParallel')
library('R.matlab')
library('Matrix')

## colour-blind-friendly palette
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
mmtoin <- 0.0393701

workdir <- 'ising-results1/'
dir.create(workdir)

## data files and bin width
dataf <- c('data/BEN_T7C1.mat', 'data/BEN_T8C1.mat')
binwidth <- 0.003

calcevidence <- function(dir='./',dataf,par1,par2){
	slabel <- substring(dir,nchar(dir))
    if(slabel!='/'){dir <- paste0(dir,'/')}
    ## read the spike times and find corresponding bins
    nn <- 2
    spikes <- sapply(1:nn,function(i){readMat(paste0(dir,dataf[i]))$cellTS %/% binwidth})

    ## set first bin at 1 and find number of bins
    tempmin <- min(c(unlist(spikes)))
    n <- max(c(unlist(spikes)))-tempmin+1

    ## sparse matrix with activity sequences (logical); one column per neuron
    train <- sparseMatrix(i=(unlist(spikes)-tempmin+1), j=unlist(sapply(1:nn,function(i){rep(i,length(spikes[[i]]))})))

    ## sparse vector with sequence of states
    states <- as(c(Inf,train[,1]+2*train[,2]),"sparseVector")

    ## frequency of state s after obs observations
    freq <- function(s,obs){substates <- states[1:(obs+1)]
        length(substates[c(substates)==s])}

    ## log-evidence after obs observations
    logevidence <- function(obs,l,f){lgamma(l)-lgamma(l+obs)+
                                  sum(sapply(0:3,function(i){lgamma(l*f+freq(i,obs))-lgamma(l*f)}))}

    lev1 <- logevidence(n,par1[1],par1[2])
    lev2 <- logevidence(n,par2[1],par2[2])
    c(lev1, lev2, lev1-lev2)
}

collectevidences <- function(datadir,par1,par2,nsamples,savedir='./',label='',seed=999){
    set.seed(seed)
	slabel <- substring(datadir,nchar(datadir))
    if(slabel!='/'){datadir <- paste0(datadir,'/')}
	slabel <- substring(savedir,nchar(savedir))
    if(slabel!='/'){savedir <- paste0(savedir,'/')}

    filelist <- list.files(path=datadir,pattern='T.*.mat')
    nfiles <- length(filelist)

    mat <- matrix(NA,nsamples,2+3)
    for(i in 1:nsamples){
	cat('\rsample ',i)
        filenumbers <- sample(1:nfiles,2)
        dataf <- filelist[filenumbers]
        mat[i,] <- c(filenumbers, calcevidence(datadir,dataf,par1,par2))
    }
	cat('\n')
    
    write.table(mat,paste0(savedir,label,'_evidences.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
	mat
}

stop('end of script')


## read the spike times and find corresponding bins
nn <- length(dataf)
spikes <- sapply(1:nn,function(i){readMat(dataf[i])$cellTS %/% binwidth})

## set first bin at 1 and find number of bins
tempmin <- min(c(unlist(spikes)))
n <- max(c(unlist(spikes)))-tempmin+1

## sparse matrix with activity sequences (logical); one column per neuron
train <- sparseMatrix(i=(unlist(spikes)-tempmin+1), j=unlist(sapply(1:nn,function(i){rep(i,length(spikes[[i]]))})))

## sparse vector with sequence of states
states <- as(c(Inf,train[,1]+2*train[,2]),"sparseVector")

freq <- function(s,obs){substates <- states[1:(obs+1)]
    length(substates[c(substates)==s])}

# freqs <- function(obs){sapply(0:3,function(i){freq(i,obs)})}

# prob <- function(s,obs,l,f){(l*f+freq(s,obs))/(l+obs)}

## a <- 0; la <- 1e-12; nu <- rep(1,4)/4
## for(i in 1:n){cat('\r',i)
##     a <- a -log(nu[states[i+1]+1])
##     nu <- (la*nu + (states[i+1]==(0:3)))/(la+1)
##     la <- la + 1}
## a/n
## saveRDS(a/n,"suprise_prior_mu.rds")

logevidence <- function(obs,l,f){lgamma(l)-lgamma(l+obs)+
                                  sum(sapply(0:3,function(i){lgamma(l*f+freq(i,obs))-lgamma(l*f)}))}

