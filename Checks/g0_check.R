
source('C:/Users/Gonde/OneDrive/Documents/R Projects/Sim_SpatialMR/SimFunctions/Simulation.R')
##create trap locations
library(secr)
traplocs<-make.grid(nx=6, ny=6, spacing = 800)
Detector<-rownames(traplocs)
traplocs<-data.frame(Detector, X=traplocs$x, Y=traplocs$y)
trapPath<-tempfile(fileext = ".csv")
write.table(x = traplocs, file = trapPath, sep=",", col.names = FALSE, row.names = FALSE) 
patht0<-"temp.csv"

masks<-list()

for (j in 1:10)
{
samps<-sim.bear(known = known, sig = sig, int.g0 =.5 , traplocs=traplocs, behav=0, IH=0, sessions=6, redun=0, stratDensity=0)
write.table(samps, file=patht0, sep = ",")
t0caphist<-read.capthist(captfile = patht0, trapfile = trapPath, detector = 'proximity')
fitted<-secr.fit(t0caphist, model = list(g0~1), buffer = 3000, trace = FALSE, CL=TRUE, detectfn = 0, start = list(.4, 0, 650))
masks[1][j]<-suggest.buffer(fitted)

samps<-sim.bear(known = known, sig = sig, int.g0 =.5 , traplocs=traplocs, behav=0, IH=0, sessions=6, redun=0, stratDensity=0)
write.table(samps, file=patht0, sep = ",")
t0caphist<-read.capthist(captfile = patht0, trapfile = trapPath, detector = 'proximity')
fitted<-secr.fit(t0caphist, model = list(g0~bk), buffer = 3000, trace = FALSE, CL=TRUE, detectfn = 0, start = list(.4, 0, 650))
masks[2][j]<-suggest.buffer(fitted)

samps<-sim.bear(known = known, sig = sig, int.g0 =.5 , traplocs=traplocs, behav=0, IH=0, sessions=6, redun=0, stratDensity=0)
write.table(samps, file=patht0, sep = ",")
t0caphist<-read.capthist(captfile = patht0, trapfile = trapPath, detector = 'proximity')
fitted<-secr.fit(t0caphist, model = list(g0~t), buffer = 3000, trace = FALSE, CL=TRUE, detectfn = 0, start = list(.4, 0, 650))
masks[3][j]<-suggest.buffer(fitted)

samps<-sim.bear(known = known, sig = sig, int.g0 =.5 , traplocs=traplocs, behav=0, IH=0, sessions=6, redun=0, stratDensity=0)
write.table(samps, file=patht0, sep = ",")
t0caphist<-read.capthist(captfile = patht0, trapfile = trapPath, detector = 'proximity')
fitted<-secr.fit(t0caphist, model = list(g0~bk+t), buffer = 3000, trace = FALSE, CL=TRUE, detectfn = 0, start = list(.4, 0, 650))
masks[4][j]<-suggest.buffer(fitted)

}

favstats(masks[[1]])
favstats(masks[[2]])
favstats(masks[[3]])
favstats(masks[[4]])