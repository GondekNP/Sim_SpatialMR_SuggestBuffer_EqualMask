BearSubsample<-function(data, type, n){
  require(dplyr)
  if (type=="SimpleRandom"){ ##Simple random (all n taken at random with no regard to period or site)
    for (i in 1:nrow(data)){
      data$uniqueID[i]<-paste(data$Period[i], data$site[i], sep="", collapse="") 
      data$INDuniqID[i]<-paste(data$Period[i], data$site[i], data$ID[i],sep="", collapse="") 
    } 
    
    selected<-mosaic::sample(data, size=n)
  }
  
  

  if (type=="Stratified"){##For each of the six periods 
    for (i in 1:nrow(data)){
      data$uniqueID[i]<-paste(data$Period[i], data$site[i], sep="", collapse="")
      data$INDuniqID[i]<-paste(data$Period[i], data$site[i], data$ID[i],sep="", collapse="") 
    }
    
    selected<-NULL
    for (i in (1:6)){         
      new<-mosaic::sample(data[data$Period==i,], size = n/6)
      selected<-rbind(selected, new)
      selected<-filter(selected, !(is.na(Period)))
    } 
  }
  
  if (type=="Spread.one"){##First take 1 from each site X period, then random
    selected<-NULL
    droprow<-NULL
    for (i in 1:nrow(data)){
      data$uniqueID[i]<-paste(data$Period[i], data$site[i], sep="", collapse="")
      data$INDuniqID[i]<-paste(data$Period[i], data$site[i], data$ID[i],sep="", collapse="") 
    }
    #data<-mutate(data, uniqueID=c(Period, site))
    data<-sample(data)##Mix them up, then take the first one from each uniqueID
    
    selected<-data[!duplicated(data$uniqueID),]
    droprow<-(!duplicated(data$uniqueID)) * c(1:nrow(data))
    data<- data[-droprow,]
    
    if((n-nrow(selected))>0)
      {(selected<-rbind(selected, mosaic::sample(data, size=(n-nrow(selected)))))}
    else{selected<-mosaic::sample(selected, size=n)}
    # selected<-selected[,1:15]
  }
  
  if (type=="Spread.two"){##First take 2 from each site X period, then random
    selected<-NULL
    droprow<-NULL
    for (i in 1:nrow(data)){
      data$uniqueID[i]<-paste(data$Period[i], data$site[i], sep="", collapse="")
      data$INDuniqID[i]<-paste(data$Period[i], data$site[i], data$ID[i],sep="", collapse="") 
    }
    data<-sample(data)##Mix them up, then take the first one from each uniqueID
    
    for (s in (1:2)){ ##Do it twice this time
    selected<-rbind(data[!duplicated(data$uniqueID),], selected)
    droprow<-(!duplicated(data$uniqueID)) * c(1:nrow(data))
    data<- data[-droprow,]
      }
    
    
    if((n-nrow(selected))>0)
    {(selected<-rbind(selected, mosaic::sample(data, size=(n-nrow(selected)))))}
    else{selected<-mosaic::sample(selected, size=n)}  }
  
  return(selected)
  }