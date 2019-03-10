#'creating comparative AC's
library(latex2exp)
egAC <- function(stratDensity, inhib, known){
library(secr)
traplocs<-make.grid(nx=6, ny=6, spacing = 800)
Detector<-rownames(traplocs)
traplocs<-data.frame(Detector, X=traplocs$x, Y=traplocs$y)
trapPath<-tempfile(fileext = ".csv")
write.table(x = traplocs, file = trapPath, sep=",", col.names = FALSE, row.names = FALSE) 

library(sp)
library(mosaic)

#'Defining 'observation window' in spatstat for rSSI to work properly
library(spatstat)
traprange<-owin(xrange=c(min(as.numeric(traplocs$X))-500, max(as.numeric(traplocs$X))+500),
                yrange=c(min(as.numeric(traplocs$Y))-500, max(as.numeric(traplocs$Y))+500))
traprange$units$singular<-"meter"
traprange$units$plural<-"meters"

#Simulating AC's for 15 bears, with inhibition range 'r' defined by homerange radius.
#If strat density is nonzero, create disproportionate grid 
if(stratDensity!=0){
  #Get some AC's for the whole grid first
  AC <- rSSI(r = inhib, n=round(length(known)*(1-stratDensity)), win = traprange, giveup = 10000)
  #Half the size of the grid
  traprange$xrange<-c(0, max(traprange$xrange)/2)
  #Generate the rest of the points on the half grid
  AC2 <- rSSI(r = inhib, n=round(length(known)*stratDensity), win = traprange, giveup = 10000)
  AC <- rbind(data.frame(AC), data.frame(AC2))
}

else {AC <- rSSI(r = inhib, n=length(known), win = traprange, giveup = 10000)}

#AC$x <- AC$x+max(traprange$xrange) ##rSSI creates origin at 0,0, so translating to 

ACs<-data.frame(AC, ID="BEAR")

return (ACs)
}

strat1<-egAC(stratDensity = .75, inhib = .2, known=c(letters, LETTERS, c("Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg","Hh"))[1:30])
strat2<-egAC(stratDensity = 0, inhib = .2, known=c(letters, LETTERS, c("Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg","Hh"))[1:30])
traplocs<-cbind(make.grid(nx=6, ny=6, spacing = 800), ID="TRAP")
strat1<-rbind(strat1,traplocs)
strat2<-rbind(strat2,traplocs)


