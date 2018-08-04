source('overlearning_simple1.R')

pp0 <- c(c(1/2,1/2)*3/4, 1/4)

prior0 <- function(t){(dnorm(t,mean=6,sd=0.5)+dnorm(t,mean=-6,sd=0.5))/2}

ext <- 11

pr0  <- function(t){c(1/(1+ext*exp(t)+exp(2*t)), 1/(1+ext*exp(-t)+exp(-2*t)),
ext/(ext+exp(-t)+exp(t)))}

test <- averagenewdata(pp=pp0,priorf=prior0,pr=pr0,nsamples=50,nshuffles=1e5,label='testsimple5np4',cores=20,jobname='j_simple2')
