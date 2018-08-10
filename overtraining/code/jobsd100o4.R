source('overlearning4.R')

pfreqs <- matrix(c(
    c(1/2,1/2)*0, 1,
    c(1/2,1/2)*1, 0
),3,2)

prior3 <- function(t){dnorm(t,mean=0,sd=1.686834785358595)}

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs=pfreqs,priorf=prior3,pr=pr2,nsamples=30,nshuffles=100000,label='sd1centr0',cores=20)
