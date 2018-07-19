## Application of the two Ising parameter priors to neurons from Nicola's data
## cross-validation-like analysis

library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
# library('plot3D')
# library('doParallel')
library('R.matlab')
library('Matrix')

## colour-blind-friendly palette
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()
mmtoin <- 0.0393701

workdir <- 'ising-results1/'
dir.create(workdir)

## data files and bin width
dataf <- c('data/BEN_T7C1.mat', 'data/BEN_T8C1.mat')
binwidth <- 0.003

## read the spike times and find corresponding bins
nn <- length(dataf)
spikes <- sapply(1:nn,function(i){readMat(dataf[i])$cellTS %/% binwidth})

## set first bin at 1 and find number of bins
tempmin <- min(c(unlist(spikes)))
n <- max(c(unlist(spikes)))-tempmin+1

## sparse matrix with activity sequences (logical); one column per neuron
train <- sparseMatrix(i=(unlist(spikes)-tempmin+1), j=unlist(sapply(1:nn,function(i){rep(i,length(spikes[[i]]))})))

## sparse vector with sequence of states
states <- as(train[,1]+2*train[,2],"sparseVector")




