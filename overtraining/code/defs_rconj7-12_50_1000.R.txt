source('overlearning_conj1.R')

relfreqs0 <- c(7*1/2, 5, 7*1/2)/12


test <- averagenewdata(relfreqs=relfreqs0,ndata=50,nsamples=1e3,k0=1,n0=1,label='rconj7-12',probs2save=1:15,samples2save=NULL,cores=20,jobname='j_conj1')
