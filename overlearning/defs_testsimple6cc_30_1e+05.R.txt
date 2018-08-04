source('overlearning_simple1.R')

pp0 <- c(c(1/2,1/2)*4/10, 6/10)

##prior0 <- function(t){(dnorm(t,mean=6,sd=0.5)+dnorm(t,mean=-6,sd=0.5))/2}

prior0 <- function(t){dnorm(t,mean=0,sd=4)}

ext <- 11

pr0  <- function(t){c(1/(1+ext*exp(t)+exp(2*t)), 1/(1+ext*exp(-t)+exp(-2*t)),
ext/(ext+exp(-t)+exp(t)))}

test <- averagenewdata(pp=pp0,priorf=prior0,pr=pr0,nsamples=30,nshuffles=1e5,label='testsimple6cc',cores=20,jobname='j_simple2')
