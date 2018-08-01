source('overlearning4.R')

pfreqs <- matrix(c(
    c(1/2,1/2)*9/11, 2/11,
    c(1/2,1/2)*2/11, 9/11
),3,2)

## nshuffles = 100 * 5e3
totals <- averagenewdata(pfreqs=pfreqs,priorf=prior2,pr=pr3,nsamples=100,nshuffles=500,label='sd100o4near',cores=20)
