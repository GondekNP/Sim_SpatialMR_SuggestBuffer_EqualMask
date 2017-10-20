source('C:/Users/Gonde/OneDrive/Documents/R Projects/Sim_SpatialMR_1.1/SimFunctions/Simulation.R')

##Genetically identified individuals (instead of 16-34-299 for eg, letters for simplicity)
known<-c(letters, LETTERS, c("Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg","Hh"))[1:30] ##30 known bears
int.g0<- .1
behav<-0
redun<-0
sessions<-2
IH<-0

##Path for trapping grid object and csv
library(secr)
traplocs<-make.grid(nx=6, ny=6, spacing = 800)
Detector<-rownames(traplocs)
traplocs<-data.frame(Detector, X=traplocs$x, Y=traplocs$y)
trapPath<-tempfile(fileext = ".csv")
write.table(x = traplocs, file = trapPath, sep=",", col.names = FALSE, row.names = FALSE) 
pathSanity<-tempfile(fileext = ".csv")

##Sigma value denoting AC's generated
sig<-sqrt((1500*1500)/pi)

runningDen<-NULL
runningModels<-list()

for (o in 1:500){

    ##Generating the dataset and capture history
    SanityData<-sim.bear(known = known, sig = sig, int.g0, traplocs=traplocs, behav, IH, sessions, redun)
    write.table(SanityData, file=pathSanity, sep = ",") 
    SanityCH<-read.capthist(captfile = pathSanity, trapfile = trapPath, detector = 'proximity')
    
    ##Fitting null model g0 ~ 1
    fitted<-secr.fit(capthist = SanityCH, model = as.formula("g0 ~ 1"), buffer = 10000, trace = FALSE, CL=TRUE, detectfn = 0)
    densityEstimate<-derived(fitted)[2,1]
    
    runningDen<- c(runningDen, densityEstimate)
    runningModels[[o]]<-fitted
    print(o)
}
##Seems like autoini is failing at .5... using -2.19, which creates a g0 of .1 (exactly where autoini starts) works just fine, and the density estimate is pretty close to what it should be

