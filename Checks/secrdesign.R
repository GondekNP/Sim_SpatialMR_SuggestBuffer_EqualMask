#'From previous simulations with my own code
library(secrdesign)
sig<-sqrt((1500*1500)/pi)
traplocs<-make.grid(nx=6, ny=6, spacing = 800)

#'Combinations of scenarios to test
scens<-make.scenarios(trapsindex = 1, noccasions = c(2,6), nrepeats = 1, D = c(.012,12), sigma=sig, g0=c(.1,.5))

#'Note: allowing secr to generate its own mask using make.mask()
ran<-run.scenarios(nrepl=100, scens, traplocs, fit=TRUE, ncores=3, seed=123124)