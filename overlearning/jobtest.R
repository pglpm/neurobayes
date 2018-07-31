source('overlearning4.R')

pfreqs <- matrix(c(1,1,9,5,5,1),3,2)/11

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs,prior,nsamples=100,nshuffles=400,label='testnos',cores=20)
