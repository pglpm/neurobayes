source('overlearning_simple1.R')

pp1 <- c(c(1/2,1/2)*4/4, 0/4)

prior1 <- function(t){dnorm(t,mean=0,sd=10)}


test <- averagenewdata(pp=pp1,priorf=prior1,pr=pr2,nsamples=50,nshuffles=100000,label='testsimple2',cores=20)
