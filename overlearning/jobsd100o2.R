source('overlearning4.R')

pfreqs <- matrix(c(1/2,1/2,10,4.5,4.5,2),3,2)/11

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs=pfreqs,priorf=prior2,pr=pr2,nsamples=75,nshuffles=5000,label='sd100o2',cores=20)
