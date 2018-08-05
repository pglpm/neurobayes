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

## states
states <- 0:2

## normalization function of conjugate prior
## this is slower
ibeta <- function(a,b){integrate(function(t){t^(a-1)*(1-t)^(b-1)},0,-1)$value}
qq2 <- function(k,n){(-1)^k*(ibeta(k, 1 - 2*n) + ibeta(-k + 2*n, 1 - 2*n))}

## this is faster
binte <- function(t,a,b){t^(a-1)*(1-t)^(b-1)}
qq <- function(k,n){
    (-1)^(k+1)*(integrate(binte,-1,0,a=k,b=1-2*n)$value
    +integrate(binte,-1,0,a=2*n-k,b=1-2*n)$value)
}


## likelihood, with a prior distribution g of (1,1,2)
g <- c(1,2,1)
likelihood <- c(function(t){1/(1+exp(t))^2},
                function(t){1/(1+cosh(t))},
                function(t){1/(1+exp(-t))^2}
                )
flikelihood <- function(t){c(1/(1+exp(t))^2,
                             1/(1+cosh(t)),
                             1/(1+exp(-t))^2)}

## argmax of conjugate prior
tmax <- function(k,n){-log((2*n-k)/k)}

maxlikelihood <- c(function(k,n){((2*n-k)/(2*n))^2},
                   function(k,n){k*(2*n-k)/(2*n^2)},
                   function(k,n){(k/(2*n))^2}
                   )

fmaxlikelihood <- function(k,n){c(((2*n-k)/(2*n))^2,
                                   k*(2*n-k)/(2*n^2),
                                   (k/(2*n))^2)}


probsfromdata <- function(data,k0,n0){
    ndata <- length(data)
    data <- data+1
    ## predictive probabiliites
    predprobs <- matrix(NA,3,ndata)

    ## predictive probabilities from max likelihood approx
    mlpredprobs <- matrix(NA,3,ndata)

    ## probability of outcome
    outcomeprobs <- rep(NA,ndata)

    ## probability of outcome, max-likelihood approx
    outcomemlprobs <- rep(NA,ndata)

    ## utility scores (unit-diagonal utility matrix)
    scores <- rep(NA,ndata)

    ## log-evidences
    logevidences <- rep(NA,ndata)

    ## frequencies
    freqs <- rep(0,3)
    
    qq0 <- qq(k0,n0)
    oldqq <- qq0    
    
    for(d in 1:ndata){
        kk <- k0 + sum(states * freqs)
        nn <- n0+d-1

        newqq <- g*sapply(states,function(i){prod(g^freqs)*qq(kk+i,nn+1)})

        ## predictive probabilities & their max-likelihood approximations
        prob <- newqq/oldqq
        predprobs[,d] <- prob
        mlpredprobs[,d] <- fmaxlikelihood(kk,nn)

        outcome <- data[d]
        outcomeprobs[d] <- prob[outcome]
        outcomemlprobs[d] <- mlpredprobs[outcome,d]
        
        ## utility
        ma <- max(prob)
        scores[d] <- (prob[outcome]==ma)/length(prob[prob==ma])

        oldqq <- newqq[[outcome]]
        logevidences[d] <- log(oldqq/qq0)
        freqs[outcome] <- freqs[outcome]+1
}
    list(predprobs=predprobs,mlpredprobs=mlpredprobs,outcomeprobs=outcomeprobs,outcomemlprobs=outcomemlprobs,scores=scores,logevidences=logevidences,finfreqs=freqs)
}

generatedata <- function(ndata,relfreqs){sample(states,ndata,replace=T,prob=relfreqs)}

averagenewdata <- function(relfreqs,ndata,nsamples,k0=1,n0=1,probs2save=NULL,samples2save=NULL,label='',seed=666,cores=1,jobname=NULL){

    if(label==''){label=format(Sys.time(),'%y%m%dT%H%M')}
    label <- paste0(label,'_',ndata,'_',nsamples)

    if(is.null(probs2save)){probs2save <- 1:min(15,ndata)}
    if(is.null(samples2save)){samples2save <- nsamples}

    if(is.character(jobname)){file.copy(paste0(jobname,'.R'),paste0('defs_',label,'.R.txt'))}

    write.table(relfreqs,paste0('relfreqs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    message('starting calculations...')
    if(cores>1){
        message('in parallel')
        cl <- makeForkCluster(cores)
        registerDoParallel(cl)
    }
    allresults <- foreach(s=1:nsamples) %dopar% {
        set.seed(seed+s)
        
        data <- generatedata(ndata,relfreqs)
        probsfromdata(data,k0,n0)
    }
    if(is.integer(cores) && cores>1){
    stopCluster(cl)
    }
    message('...done. Saving data...')
    
    allresults <- do.call(rbind,allresults)
    saveRDS(allresults,paste0('_results_',label,'.rds'))
    
    allpredprobs <- unlist(allresults[,1])
    dim(allpredprobs) <- c(3,ndata,nsamples)
    avgpredprobs <- apply(allpredprobs,c(1,2),mean,na.rm=T)
    write.table(avgpredprobs,paste0('avgpredprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgpredprobs)
    stdpredprobs <- apply(allpredprobs,c(1,2),sd,na.rm=T)
    write.table(stdpredprobs,paste0('stdpredprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdpredprobs)
    rm(allpredprobs)
    
    allmlpredprobs <- unlist(allresults[,2])
    dim(allmlpredprobs) <- c(3,ndata,nsamples)
    avgmlpredprobs <- apply(allmlpredprobs,c(1,2),mean,na.rm=T)
    write.table(avgmlpredprobs,paste0('avgmlpredprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgmlpredprobs)
    stdmlpredprobs <- apply(allmlpredprobs,c(1,2),sd,na.rm=T)
    write.table(stdmlpredprobs,paste0('stdmlpredprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdmlpredprobs)
    rm(allmlpredprobs)
    
    alloutcomeprobs <- unlist(allresults[,3])
    dim(alloutcomeprobs) <- c(ndata,nsamples)
    avgoutcomeprobs <- apply(alloutcomeprobs,1,mean,na.rm=T)
    write.table(avgoutcomeprobs,paste0('avgoutcomeprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgoutcomeprobs)
    stdoutcomeprobs <- apply(alloutcomeprobs,1,sd,na.rm=T)
    write.table(stdoutcomeprobs,paste0('stdoutcomeprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdoutcomeprobs)
    avgsurprises <- apply(-log(alloutcomeprobs),1,mean,na.rm=T)
    write.table(avgsurprises,paste0('avgsurprises_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgsurprises)
    stdsurprises <- apply(-log(alloutcomeprobs),1,sd,na.rm=T)
    write.table(stdsurprises,paste0('stdsurprises_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdsurprises)
    write.table(alloutcomeprobs[probs2save,samples2save],paste0('probsequence_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')

    rm(alloutcomeprobs)

    alloutcomemlprobs <- unlist(allresults[,4])
    dim(alloutcomemlprobs) <- c(ndata,nsamples)
    avgoutcomemlprobs <- apply(alloutcomemlprobs,1,mean,na.rm=T)
    write.table(avgoutcomemlprobs,paste0('avgoutcomemlprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgoutcomemlprobs)
    stdoutcomemlprobs <- apply(alloutcomemlprobs,1,sd,na.rm=T)
    write.table(stdoutcomemlprobs,paste0('stdoutcomemlprobs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdoutcomemlprobs)
    rm(alloutcomemlprobs)
    
    allscores <- unlist(allresults[,5])
    dim(allscores) <- c(ndata,nsamples)
    avgscores <- apply(allscores,1,mean,na.rm=T)
    write.table(avgscores,paste0('avgscores_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgscores)
    stdscores <- apply(allscores,1,sd,na.rm=T)
    write.table(stdscores,paste0('stdscores_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdscores)
    rm(allscores)
    
    alllogevidences <- unlist(allresults[,6])
    dim(alllogevidences) <- c(ndata,nsamples)
    avglogevidences <- apply(alllogevidences,1,mean,na.rm=T)
    write.table(avglogevidences,paste0('avglogevidences_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avglogevidences)
    stdlogevidences <- apply(alllogevidences,1,sd,na.rm=T)
    write.table(stdlogevidences,paste0('stdlogevidences_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdlogevidences)
    rm(alllogevidences)

    allfinfreqs <- unlist(allresults[,7])
    dim(allfinfreqs) <- c(3,nsamples)
    avgfinfreqs <- apply(allfinfreqs,1,mean,na.rm=T)
    write.table(avgfinfreqs,paste0('avgfinfreqs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(avgfinfreqs)
    stdfinfreqs <- apply(allfinfreqs,1,sd,na.rm=T)
    write.table(stdfinfreqs,paste0('stdfinfreqs_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Null')
    rm(stdfinfreqs)
    rm(allfinfreqs)
    
    message('Finished.')
    NULL
}

if(FALSE){
stop('OK: end of script')


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
}
