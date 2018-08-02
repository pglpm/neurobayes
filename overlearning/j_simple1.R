source('overlearning_simple1.R')

pp1 <- c(c(1/2,1/2)*3/4,1/4)

prior1 <- function(t){dnorm(t,mean=0,sd=100)}


test <- averagenewdata(pp=pp1,priorf=prior1,pr=pr2,nsamples=100,nshuffles=1000,label='testsimple1',cores=20)
