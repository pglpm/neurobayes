source('overlearning4.R')

pfreqs <- matrix(c(
    c(1/2,1/2)*1/4, 3/4,
    c(1/2,1/2)*4/5, 1/5
),3,2)

prior3 <- function(t){dnorm(t,mean=0,sd=1)}

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs=pfreqs,priorf=prior3,pr=pr2,nsamples=20,nshuffles=5000,label='sd1o5',cores=20)
