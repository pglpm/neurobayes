## Simple study of overlearning

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('RColorBrewer')
library('png')
library('pracma')
#library('plot3D')
#library('ggplot2')
#library('cowplot')
#library('doParallel')
#library('GA')

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


workdir <- 'overlearning/'
suppressWarnings(dir.create(workdir))


## normalization (not used)
z <- function(t){1+exp(2*t)+exp(t)}

allp <- function(p){c(p,1-sum(p))}
norma <- function(p){p/sum(p)}

## two of the three probabilities for the likelihood
pr <- function(t){c(1/(1+exp(2*t)+exp(t)), 1/(exp(-2*t)+exp(-t)+1), 1/(exp(-t)+1+exp(t)))}

prior <- function(t){dnorm(t,mean=0,sd=100)}

prfromdata <- function(data,priorf,pp=rep(1/2,2)){
    ldata <- length(data[1,])
    ## current evidence: each row = class, each col = first two probs
    integ <- matrix(NA,2,2)
    ## likelihood sequence: each row = class, each col = first two probs
    likelihood <- array(NA,c(2,2,ldata))
    ## probs for classes
    prob <- matrix(NA,2,ldata)
    evidence <- rep(1,2)
    logevidences <- rep(NA,ldata)
    ## frequencies: each row = class, each col = frequencies
    fr <- matrix(0,2,3)
    ## utility scores
    score <- rep(NA,ldata)
   ## print(integ);print(fr);print(evidence);print('')
    for(d in 1:ldata){
        integrand <- function(t,i,h){pr(t)[i] * pr(t)[1]^fr[h,1] * pr(t)[2]^fr[h,2] * pr(t)[3]^fr[h,3] * priorf(t)}
        vintegrand <- Vectorize(integrand,'t')
##        integrand <- function(t,i,h){pr(t)[i] * prod(allp(pr(t))^(fr[h,])) * priorf(t)}
       invisible(capture.output({ integ<- sapply(1:2,function(i){
            sapply(1:2,
                   function(h){((integrate(vintegrand,-Inf, Inf, abs.tol=1e-52,i=i,h=h)$value))})}) }))
        
        likelihood[,,d] <- integ/evidence
        class <- data[1,d]
        outcome <- data[2,d]
        ## probabilities for the classes
        prob[,d] <- norma(apply(likelihood[,,d],1,allp)[outcome,] * pp)

        ## utility
        score[d] <- (sign(prob[class,d]-0.5)+1)/2

        evidence[class] <- c(integ[class,],
                             evidence[class]-sum(integ[class,]))[outcome]
        logevidences[d] <- sum(log(evidence))
        fr[class,outcome] <- fr[class,outcome]+1
        ##print(integ);print(likelihood[,,d]);print(fr);print(evidence);print('')
    }
    list(likelihoods=likelihood,probs=prob,scores=score,logevidences=logevidences,finfreq=fr)
}

generatedata <- function(nsamples,pfreqs,pp=rep(1,2)/2){
    data <- matrix(NA,2,nsamples)
    data[1,] <- sample(1:2,nsamples,replace=T,prob=pp)
    le <- sapply(1:2,function(i){sum(data[1,]==i)})
for(i in 1:2){
    data[2,data[1,]==i] <- sample(1:3,le[i],replace=T,prob=pfreqs[,i])}
    data}

averagefromdata <- function(pfreqs,priorf,nsamples=100,nshuffles=100,label='',pp=rep(1/2,2),seed=999){
    if(label==''){label=format(Sys.time(),'%y%m%dT%H%M')}
    set.seed(seed)
    pb <- txtProgressBar(1,nshuffles,1,style=3)

    data <- generatedata(nsamples,pfreqs,pp)

    alllikelihood <- array(NA,c(2,2,nsamples,nshuffles))
    allscores <- matrix(NA,nsamples,nshuffles)
    alllogevidences <- matrix(NA,nsamples,nshuffles)
    res <- list()
    
    for(s in 1:nshuffles){
        sdata <- data[,sample(1:nsamples)]
        res[[s]] <- prfromdata(sdata,priorf,pp)

        alllikelihood[,,,s] <- res[[s]]$likelihoods
        allscores[,s] <- res[[s]]$scores
        alllogevidences[,s] <- res[[s]]$logevidences
        setTxtProgressBar(pb,s)
    }
    close(pb)
    avglikelihood1 <- apply(alllikelihood[1,,,],c(1,2),mean)
    avglikelihood2 <- apply(alllikelihood[2,,,],c(1,2),mean)
    avgscore <- apply(allscores,1,mean)
    avglogevidence <- apply(alllogevidences,1,mean)

    saveRDS(res,paste0('results_',label,'_',nsamples,'_',nshuffles,'.rds'))
    write.table(avglikelihood1,paste0('lh1_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    write.table(avglikelihood2,paste0('lh2_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    write.table(avgscore,paste0('scores_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    write.table(avglogevidence,paste0('logev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    write.table(res[[1]]$finfreq,paste0('finalfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    
    res
}

pfreqs <- matrix(c(1,1,8,4,4,2),3,2)/10

## nshuffles = 100 * 5e3
totals <- averagefromdata(pfreqs,prior,nsamples=100,nshuffles=500,label='std100')
## 3 -> abs.tol=1e-52
## std100 -> prior with std 100
