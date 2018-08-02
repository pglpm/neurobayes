## Simple study of overlearning - no classification problem
## Using parallel processing

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('RColorBrewer')
library('png')
library('pracma')
library('foreach')
library('doParallel')
#library('plot3D')
#library('ggplot2')
#library('cowplot')
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
pr1 <- function(t){c(1/(1+exp(2*t)+exp(t)), 1/(exp(-2*t)+exp(-t)+1), 1/(exp(-t)+1+exp(t)))}

pr2 <- function(t){c(1/(1+exp(t))^2, 1/(1+exp(-t))^2, 1/(1+cosh(t)))}

pr3 <- function(t){c(1/(1+5*exp(t)+exp(2*t)), 1/(1+5*exp(-t)+exp(-2*t)), 5/(5+exp(-t)+exp(t)))}

prior2 <- function(t){dnorm(t,mean=0,sd=100)}


prfromdatafull <- function(data,priorf,pr){
    ldata <- length(data)
    ## current evidence
    integ <- rep(NA,3)
    ## probs for outcomes
    prob <- matrix(NA,3,ldata)
    
    evidence <- 1
    logevidence <- rep(NA,ldata)
    ## frequencies
    fr <- rep(0,3)
    ## utility scores
    score <- rep(NA,ldata)
    ## surprises
    surprise <- rep(NA,ldata)
    
    for(d in 1:ldata){
        integrand <- function(t,i){pr(t)[i] * pr(t)[1]^fr[1] * pr(t)[2]^fr[2] * pr(t)[3]^fr[3] * priorf(t)}
        vintegrand <- Vectorize(integrand,'t')
##        integrand <- function(t,i,h){pr(t)[i] * prod(allp(pr(t))^(fr[h,])) * priorf(t)}
       invisible(capture.output({ integ <- sapply(1:3,function(i){integrate(vintegrand,-Inf, Inf, abs.tol=0,i=i)$value})}))
        
        probt <- integ/evidence
        prob[,d] <- probt
        
        outcome <- data[d]

        ## surprise
        surprise[d] <- -log(probt[outcome])

        ## utility
        ma <- max(probt)
        score[d] <- (probt[outcome]==ma)/length(probt[probt==ma])

        evidence <- integ[outcome]
        logevidence[d] <- log(evidence)

        fr[outcome] <- fr[outcome]+1
        ##print(integ);print(likelihood[,,d]);print(fr);print(evidence);print('')
    }
    list(probs=prob,scores=score,surprises=surprise,logevidences=logevidence,finfreq=fr)
}

generatedata <- function(nsamples,pp){sample(1:3,nsamples,replace=T,prob=pp)}

averagenewdata <- function(pp,priorf,pr,nsamples,nshuffles=2,label='',seed=999,cores=1){
    if(label==''){label=format(Sys.time(),'%y%m%dT%H%M')}
#    pb <- txtProgressBar(1,nshuffles,1,style=3)

    write.table(pp,paste0('targetfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    message('starting calculations...')
    if(cores>1){
        message('in parallel')
        cl <- makeForkCluster(cores)
        registerDoParallel(cl)
    }
    allres <- foreach(s=1:nshuffles) %dopar% {
        set.seed(seed+s)
        sdata <- generatedata(nsamples,pp)
        res <- prfromdatafull(sdata,priorf,pr)
        res
    }
    if(is.integer(cores) && cores>1){
    stopCluster(cl)
    }
    message('...done')
    
    lallres <- do.call(rbind,allres)
    
    allprobs <- unlist(lallres[,1])
    dim(allprobs) <- c(3,nsamples,nshuffles)
    avgprobs <- apply(allprobs,c(1,2),mean,na.rm=T)
    sdprobs <- apply(allprobs,c(1,2),sd,na.rm=T)

    allscores <- unlist(lallres[,2])
    dim(allscores) <- c(nsamples,nshuffles)
    avgscore <- apply(allscores,1,mean,na.rm=T)
    sdscore <- apply(allscores,1,sd,na.rm=T)

    allsurprises <- unlist(lallres[,3])
    dim(allsurprises) <- c(nsamples,nshuffles)
    avgsurprise <- apply(allsurprises,1,mean,na.rm=T)
    sdsurprise <- apply(allsurprises,1,sd,na.rm=T)

    avghit <- apply(exp(-allsurprises),1,mean,na.rm=T)
    sdhit <- apply(exp(-allsurprises),1,sd,na.rm=T)

    alllogevidences <- unlist(lallres[,4])
    dim(alllogevidences) <- c(nsamples,nshuffles)
    avglogevidence <- apply(alllogevidences,1,mean,na.rm=T)
    sdlogevidence <- apply(alllogevidences,1,sd,na.rm=T)

    allfreqs <- unlist(lallres[,5])
    dim(allfreqs) <- c(3,nshuffles)
    avgfreqs <- apply(allfreqs,1,mean,na.rm=T)
    sdfreqs <- apply(allfreqs,1,sd,na.rm=T)

    message('saving data...')
    
    saveRDS(lallres,paste0('_results_',label,'_',nsamples,'_',nshuffles,'.rds'))

    write.table(avgscore,paste0('avgscores_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdscore,paste0('sdscores_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    
    write.table(avglogevidence,paste0('avglogev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdlogevidence,paste0('sdlogev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avgprobs,paste0('avgprobs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdprobs,paste0('sdprobs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avgsurprise,paste0('avgsurprise_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdsurprise,paste0('sdsurprise_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avghit,paste0('avghit_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdhit,paste0('sdhit_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avgfreqs,paste0('avgfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdfreqs,paste0('sdfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    message('Finished.')
    list(avgscore=avgscore,avglogevidence=avglogevidence,avgprobs=avgprobs,avgsurprise=avgsurprise,avgfreqs=avgfreqs)
}





### below: not used anymore





averagefromdata <- function(pfreqs,priorf,pr,nsamples=100,nsubsamples=NULL,nshuffles=100,label='',pp=rep(1/2,2),seed=999,cores=20){
    if(label==''){label=format(Sys.time(),'%y%m%dT%H%M')}
    if(is.null(nsubsamples)){nsubsamples <- nsamples}
    set.seed(seed)
#    pb <- txtProgressBar(1,nshuffles,1,style=3)

    data <- generatedata(nsamples,pp)

    message('starting parallel calculations...')
    cl <- makeForkCluster(cores)
    registerDoParallel(cl)

    allres <- foreach(s=1:nshuffles) %dopar% {
        set.seed(seed+s)
        sdata <- data[,sample(1:nsamples)[1:nsubsamples]]
        res <- prfromdatafull(sdata,priorf,pr,pp)
        if(s==1){
            write.table(res$finfreq,paste0('finalfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
            write.table(pfreqs,paste0('targetfreqs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
        }
        list(res$likelihoods, res$scores, res$logevidences,res$probs,res$surprises)
    }
    stopCluster(cl)
    message('...done')
    
    lallres <- do.call(rbind,allres)
    
    alllikelihood <- unlist(lallres[,1])
    dim(alllikelihood) <- c(2,3,nsamples,nshuffles)
    avglikelihood1 <- apply(alllikelihood[1,,,],c(1,2),mean,na.rm=T)
    avglikelihood2 <- apply(alllikelihood[2,,,],c(1,2),mean,na.rm=T)
    sdlikelihood1 <- apply(alllikelihood[1,,,],c(1,2),sd,na.rm=T)
    sdlikelihood2 <- apply(alllikelihood[2,,,],c(1,2),sd,na.rm=T)

    allscores <- unlist(lallres[,2])
    dim(allscores) <- c(nsamples,nshuffles)
    avgscore <- apply(allscores,1,mean,na.rm=T)
    sdscore <- apply(allscores,1,sd,na.rm=T)

    alllogevidences <- unlist(lallres[,3])
    dim(alllogevidences) <- c(nsamples,nshuffles)
    avglogevidence <- apply(alllogevidences,1,mean,na.rm=T)
    sdlogevidence <- apply(alllogevidences,1,sd,na.rm=T)

    allprobs <- unlist(lallres[,4])
    dim(allprobs) <- c(2,nsamples,nshuffles)
    avgprobs <- apply(allprobs,c(1,2),mean,na.rm=T)
    sdprobs <- apply(allprobs,c(1,2),sd,na.rm=T)

    allsurprises <- unlist(lallres[,5])
    dim(allsurprises) <- c(nsamples,nshuffles)
    avgsurprise <- apply(allsurprises,1,mean,na.rm=T)
    sdsurprise <- apply(allsurprises,1,sd,na.rm=T)

    message('saving data...')
    
    saveRDS(lallres,paste0('_results_',label,'_',nsamples,'_',nshuffles,'.rds'))

    write.table(avglikelihood1,paste0('avglh1_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdlikelihood1,paste0('sdlh1_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    
    write.table(avglikelihood2,paste0('avglh2_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdlikelihood2,paste0('sdlh2_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    
    write.table(avgscore,paste0('avgscores_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdscore,paste0('sdscores_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    
    write.table(avglogevidence,paste0('avglogev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdlogevidence,paste0('sdlogev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avgprobs,paste0('avgprobs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdprobs,paste0('sdprobs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avgsurprise,paste0('avgsurprise_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(sdsurprise,paste0('sdsurprise_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    message('Finished.')
    lallres
}

# pfreqs <- matrix(c(1,1,8,4,4,2),3,2)/10
#totals <- averagefromdata(pfreqs,prior,nsamples=100,nshuffles=50000,label='testpar')

## nshuffles = 100 * 5e3
## 3 -> abs.tol=1e-52
## std100 -> prior with std 100



recalculate <- function(label,nsamples,nshuffles){
    message('reading data...')
    lallres <- do.call(rbind,readRDS(paste0('_results_',label,'_',nsamples,'_',nshuffles,'.rds')))
    message('...done')
    
    alllikelihood <- unlist(lallres[,1])
    dim(alllikelihood) <- c(2,2,nsamples,nshuffles)
    avglikelihood1 <- apply(alllikelihood[1,,,],c(1,2),mean,na.rm=T)
    avglikelihood2 <- apply(alllikelihood[2,,,],c(1,2),mean,na.rm=T)

    alllogevidences <- unlist(lallres[,4])
    dim(alllogevidences) <- c(nsamples,nshuffles)
    avglogevidence <- apply(alllogevidences,1,mean,na.rm=T)

    allprobs <- unlist(lallres[,2])
    dim(allprobs) <- c(2,nsamples,nshuffles)
    avgprobs <- apply(allprobs,c(1,2),mean,na.rm=T)

    write.table(avglikelihood1,paste0('lh1_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    write.table(avglogevidence,paste0('nalogev_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    write.table(avgprobs,paste0('probs_',label,'_',nsamples,'_',nshuffles,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    message('Finished.')
    NA
}

debugprfromdata <- function(data,priorf,pr,pp=rep(1/2,2)){
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
    ## surprises
    surprise <- rep(NA,ldata)
    ## print(integ);print(fr);print(evidence);print('')
    #if(verb=1){ fileConn<-file("intlog.txt")}
    
    for(d in 1:ldata){
        integrand <- function(t,i,h){pr(t)[i] * pr(t)[1]^fr[h,1] * pr(t)[2]^fr[h,2] * pr(t)[3]^fr[h,3] * priorf(t)}
        vintegrand <- Vectorize(integrand,'t')
##        integrand <- function(t,i,h){pr(t)[i] * prod(allp(pr(t))^(fr[h,])) * priorf(t)}
       invisible(capture.output({ integ<- sapply(1:2,function(i){
            sapply(1:2,
                   function(h){((integrate(vintegrand,-Inf, Inf, abs.tol=1e-82,i=i,h=h)$value))})}) }))

               invisible(capture.output({ integ2<- sapply(1:3,function(i){
            sapply(1:2,
                   function(h){((integrate(vintegrand,-Inf, Inf, abs.tol=1e-82,i=i,h=h)$value))})}) }))

        
        likelihood[,,d] <- integ/evidence
        print(d)
        print('this:')
        print(cbind(integ,evidence-apply(integ,1,sum)))
        print(integ2)
        print('diff:')
        print(integ2-cbind(integ,evidence-apply(integ,1,sum)))
        print('')
        print(evidence)
        print(integ/evidence)
        class <- data[1,d]
        outcome <- data[2,d]
        ## probabilities for the classes
        prob[,d] <- norma(apply(likelihood[,,d],1,allp)[outcome,] * pp)

        # surprise
        surprise[d] <- -log(prob[class,d])

        ## utility
        score[d] <- (sign(prob[class,d]-0.5)+1)/2

        evidence[class] <- c(integ[class,],
                             evidence[class]-sum(integ[class,]))[outcome]
        logevidences[d] <- sum(log(evidence))

        fr[class,outcome] <- fr[class,outcome]+1
        ##print(integ);print(likelihood[,,d]);print(fr);print(evidence);print('')
    }
    list(likelihoods=likelihood,probs=prob,scores=score,surprises=surprise,logevidences=logevidences,finfreq=fr)
}
