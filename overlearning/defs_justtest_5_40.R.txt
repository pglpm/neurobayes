source('overlearning_conj1.R')

relfreqs0 <- c(1/2*9, 1, 1/2*9)/10


test <- averagenewdata(relfreqs=relfreqs0,ndata=5,nsamples=40,k0=1,n0=1,label='justtest',probs2save=1:5,samples2save=NULL,cores=20,jobname='j_conj1')
