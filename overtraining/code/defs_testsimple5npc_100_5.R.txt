source('overlearning_simple1.R')

pp0 <- c(c(1/2,1/2)*2/3, 1/3)

prior0 <- function(t){(dnorm(t,mean=5,sd=0.5)+dnorm(t,mean=-5,sd=0.5))/2}

ext <- 11

pr0  <- function(t){c(1/(1+ext*exp(t)+exp(2*t)), 1/(1+ext*exp(-t)+exp(-2*t)),
ext/(ext+exp(-t)+exp(t)))}

test <- averagenewdata(pp=pp0,priorf=prior0,pr=pr0,nsamples=100,nshuffles=5,label='testsimple5npc',cores=20,jobname='j_simple2')
