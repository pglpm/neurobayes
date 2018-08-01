source('overlearning4.R')

pfreqs <- matrix(c(1/2,1/2,10,4.5,4.5,2),3,2)/11

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs=pfreqs,priorf=prior2,pr=pr2,nsamples=10,nshuffles=1e5,label='sd100o3',cores=20)
