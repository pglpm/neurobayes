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
