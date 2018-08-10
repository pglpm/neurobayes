source('overlearning_simple1.R')

pp1 <- c(c(1/2,1/2)*4/4, 0/4)

prior4 <- function(t){(dnorm(t,mean=5,sd=0.5)+dnorm(t,mean=-5,sd=0.5))/2}

test <- averagenewdata(pp=pp1,priorf=prior4,pr=pr2,nsamples=20,nshuffles=32768*5,label='testsimple3np',cores=20)
