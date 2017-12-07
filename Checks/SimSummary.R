compile.secr.results <- function (whichSource="EmpData"){
library(secr)
  library(stringr)
compiled<-data.frame()
trials<-c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")
if (whichSource=="EmpData"){trials<-"t1"}
for(j in trials){
  for (k in c("g0 tilde bk", "g0 tilde 1", "g0 tilde bk + t", "g0 tilde t")){
    for (l in c("250", "550", "850")){
      path<-paste(whichSource,"/", j, "/", k,"/",l, sep="", collapse="")
      files<-list.files(path)
      whichDesktop<-which(files=="desktop.ini")
      whichIcon<-which(files=="Icon\r")
      if (length(whichDesktop)>0){files<-files[-whichDesktop]}
      if (length(whichIcon)>0){files<-files[-whichIcon]}
      files<-paste(path, files, sep="/")
      
      for (m in 1:length(files)){
        newSECR<-readRDS(files[m])
        try({colnames(newSECR$subsamps)[which(colnames(newSECR$subsamps)=="Occassion")]<-"Period"})
        try({colnames(newSECR$fullsamps)[which(colnames(newSECR$fullsamps)=="Occassion")]<-"Period"})
        
        try({colnames(newSECR$subsamps)[which(colnames(newSECR$subsamps)=="Detector")]<-"site"})
        try({colnames(newSECR$fullsamps)[which(colnames(newSECR$fullsamps)=="Detector")]<-"site"})
        
        
       for (i in 1:nrow(newSECR$fullsamps)){
          newSECR$fullsamps$uniqueID[i]<-paste(newSECR$fullsamps$Period[i], newSECR$fullsamps$site[i], sep="", collapse="")
          newSECR$fullsamps$INDuniqID[i]<-paste(newSECR$fullsamps$Period[i], newSECR$fullsamps$site[i], newSECR$fullsamps$ID[i],sep="", collapse="") 
       }
        if (length(newSECR$subsamps)>0){
        for (i in 1:nrow(newSECR$subsamps)){
          newSECR$subsamps$uniqueID[i]<-paste(newSECR$subsamps$Period[i], newSECR$subsamps$site[i], sep="", collapse="")
          newSECR$subsamps$INDuniqID[i]<-paste(newSECR$subsamps$Period[i], newSECR$subsamps$site[i], newSECR$subsamps$ID[i],sep="", collapse="") 
          }
        }
        
        
        #Redundancy info
        FullN.notRedun=sum(!duplicated(newSECR$fullsamps$INDuniqID))
        SubN.notRedun=sum(!duplicated(newSECR$subsamps$INDuniqID))
        SubProp.notRedun=SubN.notRedun/as.numeric(paste(l))
        
        #Density info
        DenhatBoth<-NA
        Denhat<-NA
        Denhat.SE<-NA
        DenhatBoth<-newSECR$derived
        
       
       

        trial<-j
        model<-k
        size<-l
        
        
        if (!is.null(DenhatBoth[[1]][1])&&!is.na(DenhatBoth[[1]][1])){
        Denhat = DenhatBoth[2,1]
        Denhat.SE = DenhatBoth[2,2]
        fullN<-newSECR$fullN
        
          if (length(newSECR$fit$par)==2){
            g0<-newSECR$fit$par[1]
            g0.bTRUE<-NA
            sigma<-newSECR$fit$par[2]
          } else {
            g0<-newSECR$fit$par[1]
            g0.bTRUE<-newSECR$fit$par[2]
            sigma<-newSECR$fit$par[3]
          } 
        } else{g0<-NA;g0.bTRUE<-NA;sigma<-NA;fullN<-NA}
        
        #Compiling
        newLine<-data.frame(FullN.notRedun, SubProp.notRedun, SubN.notRedun, Denhat, Denhat.SE, g0, g0.bTRUE, sigma, trial, model, size, fullN, whichSource, sim=files[m])
        compiled<-rbind(compiled, newLine)
        print(newLine)
      }
    }
  }
}
return(compiled)

}

compiledEmp<-compile.secr.results()
compiledSim<-compile.secr.results("SimData")
compiled<-rbind(compiledEmp, compiledSim)
library(mosaic)

write.csv(compiledEmp, "Checks/CompiledEmp.csv")
write.csv(compiledSim, "Checks/CompiledSim.csv")

#compiledEmp<-read.csv("Checks/CompiledEmp.csv")
#compiledSim<-read.csv("Checks/CompiledSim.csv")

#need to split up sim into simnumber and subtype 
split<-strsplit(as.character(compiled$sim), split = "")
simNo<-NULL
subtype<-NULL
for (j in 1:nrow(compiled)){
  new<-split[[j]]
  subNew <- paste(new[(length(new)-20):(length(new)-9)], sep="", collapse="")
  simNoNew <- as.numeric(paste(new[(length(new)-8):(length(new)-4)], sep="", collapse=""))
  subtype<-c(subtype, subNew)
  simNo<- c(simNo, simNoNew)
}
compiled$subtype<-subtype
compiled$simNo<-simNo

#Now we need the discrepancy for each individual simulation
#Todo: fix for two diff models
discrep<-NULL
for (k in c("SimData", "EmpData")){
  for (h in c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")){
    for (i in c(250,550,850)){
      for (t in  c("g0 tilde bk", "g0 tilde 1", "g0 tilde bk + t", "g0 tilde t")){
          j<-10001
          for (u in 1:(nrow(compiled))){
              browser()
              nextSim<-filter(compiled, trial==h, simNo==j, whichSource==k, model==t, size==i)
              if(nrow(nextSim)==0){break}
              newLine<-nextSim
              if(nrow(filter(nextSim, subtype=="Full"))>0){newLine$DenhatBiasFull<-nextSim[,4]-filter(nextSim, subtype=="Full")[,"Denhat"]} #Universal secr density is animals per hectare
              #newLine$DenhatBiasTrue<-nextSim[,4]-.012 #Universal secr density is animals per hectare
              
              print(newLine$DenhatBias)

              discrep<-bind_rows(discrep, newLine)
              j<-j+1
          }
      }
    }
  }
}


##cleanup
colnames(discrep)[15]<-"SubsamplingType"
discrep$SubsamplingType[which(discrep$SubsamplingType=="0/Spread.one")]<-"SSP"
compiled$subtype[str_detect(compiled$subtype, "Spread")]<-"SPR"

discrep$SubsamplingType[which(discrep$SubsamplingType=="SimpleRandom")]<-"SRS"
compiled$subtype[which(compiled$subtype=="SimpleRandom")]<-"SRS"

compiled$subtype[str_detect(compiled$subtype, "Ful")]<-"Full"

levels(compiled$trial)<-c(levels(compiled$trial),"Empirical")
compiled$trial[which(compiled$whichSource=="EmpData")]<-"Empirical"



#Money plots

ggplot(data=filter(compiled, trial==c("t4","t5","t6","t7", "Empirical")))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=SubProp.notRedun, fill=subtype)) +
  ylab("Proportion of \nnon-redundant Samples") +
  xlab("Trial") +
  ggtitle("Sample Redundancy \nvs Subsampling Type") +
  theme_grey(14) +
  facet_wrap(~size)

ggplot(data=filter(compiled, whichSource=="SimData", trial==c("t4","t5","t6","t7")))+
  scale_fill_grey() +
  geom_boxplot(aes(x=trial, y=SubProp.notRedun, fill=subtype)) +
  ylab("Proportion of \nnon-redundant Samples") +
  xlab("Trial") +
  ggtitle("Sample Redundancy \nvs Subsampling Type (Simulated)") +
  theme_light() +
  facet_wrap(~size)

ggplot(data=filter(compiled, whichSource=="EmpData", size!="999"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Emulated)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t1"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t1)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t2"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t2)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)


ggplot(data=filter(compiled, whichSource=="SimData", trial=="t3"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t3)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t4"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t4)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t5"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t5)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t6"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t6)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t7"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t7)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)

ggplot(data=filter(compiled, whichSource=="SimData", trial=="t8"))+
  scale_fill_hue() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=subtype)) +
  ylab("Derived density") +
  xlab("Trial") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated, t8)") +
  theme_light() +
  facet_wrap(model~size,ncol=3)


ggplot(data=filter(discrep, whichSource=="SimData"))+
  scale_fill_grey() +
  geom_boxplot(aes(x=trial, y=Denhat, fill=SubsamplingType)) +
  ylab("d-hat / d") +
  ggtitle("d-hat \nvs Subsampling Type (Simulated)") +
  theme_grey(14) +
  facet_wrap(~size+model, ncol=4)

ggplot(data=filter(discrep, whichSource=="SimData"))+
  scale_fill_grey() +
  geom_boxplot(aes(x=trial, y=DenhatBias, fill=SubsamplingType)) +
  ylab("d-hat bias") +
  ggtitle("d-hat bias \nvs Subsampling Type (Simulated)") +
  theme_grey(14) +
  facet_wrap(~size+model, ncol=4)

ggplot(data=filter(discrep, whichSource=="SimData", size!=850))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=plogis(g0), fill=SubsamplingType)) +
  ylab("g0") +
  ggtitle("g0 vs Subsampling Type (Simulated)") +
  theme_grey(14) +
  facet_wrap(size~model, ncol=4) +
  geom_hline(aes(linetype=5), yintercept=.5)

#Looking into weird results at 850... seemingly random trials have wider boxes, or is they all just polymodal?
densityplot(filter(discrep, whichSource=="SimData", trial=="t6", model=="g0 tilde bk + t", size==850, SubsamplingType=="SimpleRandom")$nhatBias)

densityplot(filter(discrep, whichSource=="SimData", trial=="t6", model=="g0 tilde bk + t", size==850, SubsamplingType=="Spread.one")$nhatBias)

ggplot(data=filter(discrep, whichSource=="SimData"))+
  scale_fill_grey() +
  geom_boxplot(aes(x=trial, y=nhat.SE, fill=SubsamplingType)) +
  ylab("n-hat - Known Individuals (30)") +
  ggtitle("n-hat SE \nvs Subsampling Type (Simulated)") +
  theme_light() +
  facet_wrap(~size+model, ncol=4)


  


