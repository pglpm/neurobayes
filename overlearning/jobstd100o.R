source('overlearning4.R')

pfreqs <- matrix(c(1,1,8,4,4,2),3,2)/10

## nshuffles = 100 * 5e3
totals <- averagefromdata(pfreqs,prior,nsamples=100,nshuffles=100000,label='sd100o')
