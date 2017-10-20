#setwd("C:/Users/Gonde/OneDrive/Documents/R Projects/Sim_SpatialMR")

setwd("D:/Sim_SpatialMR_SuggestBuffer_EqualBuffer")


sim.bear <- function (known, sig, traplocs, int.g0=1, behav= -.7, IH=0, sessions=2, redun = 0, inhib=.2, stratDensity=0){
  library(sp)
  library(mosaic)
  
  #'Defining 'observation window' in spatstat for rSSI to work properly
  library(spatstat)
  #'range is 500 m buffered on the min/max of x and y locations
  traprange<-owin(xrange=c(min(as.numeric(traplocs$X))-500, max(as.numeric(traplocs$X))+500),
                  yrange=c(min(as.numeric(traplocs$Y))-500, max(as.numeric(traplocs$Y))+500))
  traprange$units$singular<-"meter"
  traprange$units$plural<-"meters"
  
  maskBase<-readRDS("MaskObjectForRSSI.rds")
  
  #Simulating AC's for 15 bears, with inhibition range 'r' defined by homerange radius.
  #If strat density is nonzero, create disproportionate grid 
  if(stratDensity!=0){
    #Get some AC's for the whole grid first
    AC <- rSSI(r = inhib, n=round(length(known)*(1-stratDensity)), win = maskBase, giveup = 10000)
    #Half the size of the grid
    traprange$xrange<-c(0, max(traprange$xrange)/2)
    #Generate the rest of the points on the half grid
    AC2 <- rSSI(r = inhib, n=round(length(known)*stratDensity), win = maskBase, giveup = 10000)
    AC <- rbind(data.frame(AC), data.frame(AC2))
  }
  
  else {AC <- rSSI(r = inhib, n=length(known), win = maskBase, giveup = 10000)}
  
  ACs<-data.frame(AC, ID=known, captured=rep(FALSE,length(known)), IHconstant = rnorm(n = length(known), mean=0, sd = IH))
  
  BearSamps<-data.frame()
  
  
  for (s in 1:sessions){ #Captures for each session
    
    for (b in known){ #Captures for each bear in each session
      bAC<-as.numeric(filter(ACs, ID==b)[,c("x","y")])
      #Euclidean distance between AC for this bear and each trap location, and subsequent half-normal capture prob
      trapMatrix<-data.matrix(traplocs[,c("X", "Y")])
      dists<-data.frame(dist=spDistsN1(pts=trapMatrix, pt = bAC), trapID=traplocs$Detector)
      
      for (h in 1:nrow(dists)){ #For each individual trap
        # intercept capture prob + effect from individual heterogeneity
        logit.g0<- int.g0 + filter(ACs, ID==b)[,"IHconstant"]
        g0<-plogis(logit.g0)
        g <- g0 * exp((-(dists$dist^2))/(2*sig^2)) #halfnormal curve as defined in secr documentation
        dists$g<-g
        
        #Now, a distribution for after a bear is captured
        logit.g0bk<- int.g0 + behav + filter(ACs, ID==b)[,"IHconstant"]
        g0bk<-plogis(logit.g0bk)
        gbk <- g0bk * exp((-(dists$dist^2))/(2*sig^2)) #halfnormal curve as defined in secr documentation
        dists$gbk<-gbk
        
        #To figure out which to use, determine whether bear has been captured at each detector and assign accordingly
        if (nrow(filter(BearSamps, ID==b, site==dists$trapID[h]))>0)
        {capProb <- dists$gbk[h]} #if captured already, give gbk capprob
        else 
        {capProb <- dists$g[h]} #if not, give naive g capprob
        
        IHc<-filter(ACs, ID==b)[,"IHconstant"]
        
        if ( rbinom(n=1, size=1, prob=capProb) == 1 ){ ##Coin flip - if captured (evals to 1), add a row to the samps, mark bear as captured
          newSamp<-data.frame(type="BearMR", ID = b, Period=s, site=dists$trapID[h]) #first (non-redundant) sample
          BearSamps<-rbind(BearSamps, newSamp)
          if (redun!=0){
            for (v in (1:(rpois(1, (redun + exp(IHc)))))) {BearSamps<-rbind(BearSamps, newSamp)} #if redun is 0, evals to 1, only one samp
          }
          ACs$captured[which(ACs$ID==b)] <- TRUE ##Bear is captured, next time the cap prob will change depending on 'behav'
        }
        
      }
      
    }
    
  }
  BearSamps$Period<-as.numeric(BearSamps$Period)
  return(BearSamps)
}

secr.from.samples<-function (full, samps, trapcsv, subtype, modEval, trial, number, sizechar, fullN, source, sizeSuccess)  { 
  sampsAnalyzed<-samps
  samps<-samps[!duplicated(samps$INDuniqID),]
  library(secr)
  fitted<-NULL
  if (sizeSuccess==TRUE){
    strt<-Sys.time()
    patht0<-tempfile(fileext = ".csv")
    
    if(source=="EmpData"){
      #Gender: Sex = 204.25 -> Male, Sex= 250.25 ->Female
      samps$Group<-rep("M",nrow(samps))
      samps$Group[samps$Sex==250.25]<-"F" 
      
      #Create a data set that contains a count of the number of times each bear
      #was seen for each unique site x period combination. 
      caphist<- samps%>%group_by(ID, site, Period)%>%
        summarize(Count= n())
      sexid<-unique(select(samps, ID, Group))
      caphist2<-merge(caphist, sexid, all=FALSE)
      samps<-data.frame(Session="BearMR", ID=caphist$ID, Occassion=caphist2$Period, Detector=caphist2$site, Sex=caphist2$Group)
    }
    
    write.table(samps, file=patht0, sep = ",") 
    
    try({
      t0caphist<-read.capthist(captfile = patht0, trapfile = trapcsv, detector = 'proximity')
      buff<-suggest.buffer(t0caphist, detectfn = 0)
      fitted<-secr.fit(t0caphist, model = modEval, buffer = buff, trace = FALSE, CL=TRUE, detectfn = 0)})
    
    if(!is.null(fitted)){outcome<-TRUE} else{outcome<-FALSE; print("Model fit failed.")}
    fitted$timeElapsed<-Sys.time() - strt
    fitted$sampsAnalyzed <- sampsAnalyzed
    fitted$subsamps<-samps
    fitted$outcome<-outcome
    
  } 
  fitted$fullsamps<-full
  fitted$fullN<-fullN
  fitted$sizeSuccess<-sizeSuccess
  
  #now save the object in some logical way as RDS
  modEval<-Reduce(paste, deparse(modEval))
  modelPathName<-gsub(pattern = " ~ ", replacement = " tilde " , x = modEval)
  pathRDS<-paste(source,"/", trial, "/" , modelPathName , "/", sizechar, "/", subtype, number, ".rds", collapse="", sep="")
  fitted$derived<-derived(fitted)
  print(fitted$derived)
  print(pathRDS)
  saveRDS(fitted, file = pathRDS)
}

get.rds.starts <- function(){
  starts<-data.frame()
  for(j in c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")){
    for (k in c("g0 tilde bk", "g0 tilde 1", "g0 tilde bk + t", "g0 tilde t")){
      for (l in c("250", "550", "850")) {
        for (h in c("SimData", "EmpData")){
          path<-paste(h, "/", j, "/", k, "/", l, sep="", collapse="")
          files<-list.files(path)
          
          if (length(which(files == "SimpleRandom10001.rds"))==0){
            newLine<-data.frame(trial=j,model=k, startno=10001, source=h, size=l)
          }
          else{
            maxFile<-max(files)
            maxFile<-strsplit(maxFile, "")[[1]]
            last<-maxFile[(length(maxFile)-8):(length(maxFile)-4)] ##This will only pick spread.one (or whatever the max char subsample is in the future, but that shouldn't matter because the trials are paired, and their max numbers are the same 
            last<-as.numeric(paste(last, sep="", collapse=""))
            newLine<-data.frame(trial=j,model=k, startno=last+1, source=h, size=l)
          }
          starts<-rbind(starts, newLine)
        }
      }
    }
  }
  colnames(starts)<-c("trial", "model", "startno", "source", "size")
  return(starts)
}

bear.init <- function(sim, known, trial, model, sig, size, traplocs, trapPath, behav=0, IH=0, sessions=12, redun=0, int.g0=.16, stratDensity=0){
  #' Now, fit the models in parallel like in original project
  library(doParallel)
  library(foreach)
  
  #Single simulation for both fittings 
  attempts<-0 #5 attempts to reach size desired
  sizeSuccess<-TRUE #if false, pass and save just the nrow of full samps
  fullsamps<-sim.bear(known = known, sig = sig, int.g0, traplocs=traplocs, behav, IH, sessions, redun, stratDensity)
  if (sim == TRUE) {
    while (nrow(fullsamps)<size && attempts<6){ #Make sure at least desired size
      fullsamps <- sim.bear(known = known, sig = sig, int.g0, traplocs=traplocs, behav, IH, sessions, redun, stratDensity)
      attempts <- attempts+1
      }
    
    if (nrow(fullsamps)<size){
      sizeSuccess<-FALSE
      } 
    
    }else {fullsamps <- read.csv("samps1_uniq.csv"); trapPath<-"detectorfileScaled.csv"; trial<-"t1"} #note that trial is not relevant for empirical data, but calling it t1 helps it save in the same format as simulated stuff. 
  
  
  #setup parallel backend to use all processors - note that this is not generally recommended for
  #computers that are actually in use, but since this was run primarily on a server and/or broken
  #laptop, I opted to use all cores, because I don't need any cores to run other tasks. 
  cl<-makeCluster(detectCores())
  registerDoParallel(cl)
  foreach (l=c("SimpleRandom", "Spread.one", "Full")) %dopar% {
    #for (l in c("Full", "SimpleRandom", "Spread.one")) { #just for debugging - otherwise impossible to know where errors occur
    source('SimFunctions/Simulation.R') ##source this script on each core so that all functions exist...
    source('SimFunctions/BearSubsample.R') ##subsampling functions
    j<-trial
    k<-model
    library(spatstat)
    library(secr)
    library(mosaic)
    library(sp)
    library(foreach)
    fullN<-nrow(fullsamps)
    starts<-get.rds.starts()
    notilde<-gsub(pattern = " ~ ", replacement = " tilde " , x = k)
    if(sim==TRUE){datatype<-"SimData"}else{datatype<-"EmpData"}
    #weird bug... using as.character inside of filter doesnt work?
    sizechar<-as.character(size)
    if (size==10000){ size <- nrow(fullsamps); sizechar<-"10000" }
    startno<-filter(starts, trial==j, model==notilde, source==datatype, size==sizechar)[,"startno"]
    if (sizeSuccess==TRUE & l!="Full"){
      
      hairsamps<-BearSubsample(data = fullsamps, type = l, n = size)
    
    }else{
      
      for (i in 1:nrow(fullsamps)){
        
      fullsamps$uniqueID[i]<-paste(fullsamps$Period[i], fullsamps$site[i], sep="", collapse="") 
      fullsamps$INDuniqID[i]<-paste(fullsamps$Period[i], fullsamps$site[i], fullsamps$ID[i],sep="", collapse="") 
        }
      hairsamps<-fullsamps
      
      } ##If size is large enough, proceed to subsample, if not, assign the full samps (it won't matter because secr.fit will never be called)
    secr.from.samples(full=fullsamps, samps = hairsamps, trapcsv = trapPath, modEval = as.formula(k), trial = j, number = startno, subtype = l, sizechar=sizechar, fullN=fullN, source=datatype, sizeSuccess)
    
  }
  stopCluster(cl)
}

bear.setup <- function(startstage=1){
  for(p in 1:500){
    print(paste("Completed ", p-1 , " iterations. System time: ", Sys.time(), sep="", collapse = ""))
    oneDone<-FALSE
    if (oneDone==TRUE) {rm(list=ls())}
    
    ##create trap locations
    library(secr)
    traplocs<-make.grid(nx=6, ny=6, spacing = 800)
    Detector<-rownames(traplocs)
    traplocs<-data.frame(Detector, X=traplocs$x, Y=traplocs$y)
    trapPath<-tempfile(fileext = ".csv")
    write.table(x = traplocs, file = trapPath, sep=",", col.names = FALSE, row.names = FALSE)

    ##Genetically identified individuals (instead of 16-34-299 for eg, letters for simplicity)
    known<-c(letters, LETTERS, c("Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg","Hh"))[1:30] ##30 known bears
    
    #sig<-sqrt(3/pi)/1000
    sig<-sqrt((1500*1500)/pi)
    
    if (startstage==1 || oneDone==TRUE){
      #empirical, all models and subtypes
      bear.init(sim = FALSE, known, sig, size=250, model="g0 ~ 1", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=550, model="g0 ~ 1", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=850, model="g0 ~ 1", trial = "t1", traplocs = traplocs, trapPath = trapPath)

      bear.init(sim = FALSE, known, sig, size=250, model="g0 ~ bk", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=550, model="g0 ~ bk", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=850, model="g0 ~ bk", trial = "t1", traplocs = traplocs, trapPath = trapPath)

      bear.init(sim = FALSE, known, sig, size=250, model="g0 ~ t", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=550, model="g0 ~ t", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=850, model="g0 ~ t", trial = "t1", traplocs = traplocs, trapPath = trapPath)

      
      bear.init(sim = FALSE, known, sig, size=250, model="g0 ~ bk + t", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=550, model="g0 ~ bk + t", trial = "t1", traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = FALSE, known, sig, size=850, model="g0 ~ bk + t", trial = "t1", traplocs = traplocs, trapPath = trapPath)

      oneDone<-TRUE
    }
    
    if (startstage==2 || oneDone==TRUE){
      #sim g0 ~ 1, n = 250
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ 1", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==3 || oneDone==TRUE){
      #sim g0 ~ 1, n = 550
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ 1", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==4 || oneDone==TRUE){
      #sim g0 ~ 1, n = 850
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ 1", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
  
    
    if (startstage==6 || oneDone==TRUE){
      #sim g0 ~ bk, n = 250
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==7 || oneDone==TRUE){
      #sim g0 ~ bk, n = 550
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==8 || oneDone==TRUE){
      #sim g0 ~ bk, n = 850
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
 
    if (startstage==9 || oneDone==TRUE){
      #sim g0 ~ t, n = 250
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==10 || oneDone==TRUE){
      #sim g0 ~ t, n = 550
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==11 || oneDone==TRUE){
      #sim g0 ~ t, n = 850
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    

    
    if (startstage==13 || oneDone==TRUE){
      #sim g0 ~ bk + t, n = 250
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=250, model="g0 ~ bk + t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==14 || oneDone==TRUE){
      #sim g0 ~ bk + t, n = 550
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=550, model="g0 ~ bk + t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }
    
    if (startstage==15 || oneDone==TRUE){
      #sim g0 ~ bk + t, n = 850
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t1", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t2", int.g0 = .5, behav = 1, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t3", int.g0 = .5, behav = 0, IH=1.25, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t4", int.g0 = .5, behav = 0, IH=0, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t5", int.g0 = .5, behav = 0, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t6", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath)
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t7", int.g0 = .5, behav = 1, IH=1.25, redun=1, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75) 
      bear.init(sim = TRUE, known, sig, size=850, model="g0 ~ bk + t", trial = "t8", int.g0 = .5, behav = 0, IH=0, redun=0, sessions=6, traplocs = traplocs, trapPath = trapPath, stratDensity = .75)
      oneDone<-TRUE
    }

  }
}

