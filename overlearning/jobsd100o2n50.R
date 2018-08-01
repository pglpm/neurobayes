source('overlearning4.R')

pfreqs <- matrix(c(1/2,1/2,10,4.5,4.5,2),3,2)/11

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs,prior,nsamples=20,nshuffles=1e6,label='sd100onew20',cores=20)
