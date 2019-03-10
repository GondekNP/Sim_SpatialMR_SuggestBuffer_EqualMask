# setwd('C:/Users/ConservationMetrics/Documents/NPG/Sim_SpatialMR_SuggestBuffer_EqualMask')
setwd('C:/Users/Nick/Documents/GitHub/spatialMR-master/Sim_SpatialMR_SuggestBuffer_EqualMask')
library(mosaic)
library(latex2exp)
library(readr)
library(ggplot2)
library(tidyverse)
library(secr)
library(spatstat)

maskBase<-readRDS("MaskObjectForRSSI.rds")
trueDensity = 30 / (area.owin(as.owin(maskBase)) * .0001) #1 sq meter is .0001 hectare

compile.secr.results <- function (whichSource="EmpData"){
library(secr)
  library(stringr)
compiled<-data.frame()
trials<-c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")
if (whichSource=="EmpData"){trials<-"t1"}
for(j in trials){
    for (k in c("g0 tilde bk", "g0 tilde 1", "g0 tilde bk + t", "g0 tilde t")){
    for (l in c("850", "550", "250")){
      path<-paste(whichSource,"/", j, "/", k,"/",l, sep="", collapse="")
      files<-list.files(path)
      whichDesktop<-which(files=="desktop.ini")
      whichIcon<-which(files=="Icon\r")
      if (length(whichDesktop)>0){files<-files[-whichDesktop]}
      if (length(whichIcon)>0){files<-files[-whichIcon]}
      files<-paste(path, files, sep="/")
      
      for (m in 1:length(files)){
        newSECR<-readRDS(files[m])
        subtype = str_extract(files[m], "(?<=\\d{3}/).*(?=\\d{5})")
        if(str_detect(files[m], "Full")){size = nrow(newSECR$fullsamps)}else{size = str_extract(files[m], "(?<=_).*(?=.rds)")}
        simNo = paste0(str_extract(files[m], "(?<=/)\\d{3}(?=/)"), "_", str_extract(files[m], "\\d{5}(?=.*.rds)"))
        trial<-j
        model<-k
        OGpath = paste0(whichSource, "_", trial, "_", model, "_", simNo)
        
        if (is.null(newSECR)){
          newLine<-data.frame(FullN.notRedun = NA, SubProp.notRedun = NA, SubN.notRedun = NA,
                              Denhat = NA, Denhat.SE = NA, g0 = NA, g0.bTRUE = NA, sigma = NA, trial, model,
                              size, simNo, subtype, fullN = NA, whichSource, sim=files[m], FullProp.notRedun = NA, OGpath)
          compiled<-rbind(compiled, newLine)
          print(newLine)
          next()
        }
        if (!newSECR$sizeSuccess){
          newLine<-data.frame(FullN.notRedun = NA, SubProp.notRedun = NA, SubN.notRedun = NA,
                                                      Denhat = NA, Denhat.SE = NA, g0 = NA, g0.bTRUE = NA, sigma = NA, trial, model,
                                                      size, simNo, subtype, fullN = NA, whichSource, sim=files[m], FullProp.notRedun = NA, OGpath)
        compiled<-rbind(compiled, newLine)
        print(newLine)
        next()}
        
        try({colnames(newSECR$sampsAnalyzed)[which(colnames(newSECR$sampsAnalyzed)=="Occassion")]<-"Period"})
        try({colnames(newSECR$fullsamps)[which(colnames(newSECR$fullsamps)=="Occassion")]<-"Period"})
        
        try({colnames(newSECR$sampsAnalyzed)[which(colnames(newSECR$sampsAnalyzed)=="Detector")]<-"site"})
        try({colnames(newSECR$fullsamps)[which(colnames(newSECR$fullsamps)=="Detector")]<-"site"})
        
        
       for (i in 1:nrow(newSECR$fullsamps)){
          newSECR$fullsamps$uniqueID[i]<-paste(newSECR$fullsamps$Period[i], newSECR$fullsamps$site[i], sep="", collapse="")
          newSECR$fullsamps$INDuniqID[i]<-paste(newSECR$fullsamps$Period[i], newSECR$fullsamps$site[i], newSECR$fullsamps$ID[i],sep="", collapse="") 
       }
        if (length(newSECR$sampsAnalyzed)>0){
        for (i in 1:nrow(newSECR$sampsAnalyzed)){
          newSECR$sampsAnalyzed$uniqueID[i]<-paste(newSECR$sampsAnalyzed$Period[i], newSECR$sampsAnalyzed$site[i], sep="", collapse="")
          newSECR$sampsAnalyzed$INDuniqID[i]<-paste(newSECR$sampsAnalyzed$Period[i], newSECR$sampsAnalyzed$site[i], newSECR$sampsAnalyzed$ID[i],sep="", collapse="") 
          }
        }
        
        #Redundancy info
        FullN.notRedun=sum(!duplicated(newSECR$fullsamps$INDuniqID))
        SubN.notRedun=sum(!duplicated(newSECR$sampsAnalyzed$INDuniqID))
        SubProp.notRedun=SubN.notRedun/nrow(newSECR$sampsAnalyzed)
        
        
        #Density info
        DenhatBoth<-NA
        Denhat<-NA
        Denhat.SE<-NA
        DenhatBoth<-newSECR$derived

        if (!is.null(DenhatBoth[[1]][1])&&!is.na(DenhatBoth[[1]][1])){
        Denhat = DenhatBoth[2,1]
        Denhat.SE = DenhatBoth[2,2]
        fullN<-newSECR$fullN
        FullProp.notRedun = FullN.notRedun/ fullN
        
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
     
        newLine<-data.frame(FullN.notRedun, SubProp.notRedun, SubN.notRedun,
                                  Denhat, Denhat.SE, g0, g0.bTRUE, sigma, trial, model,
                                  size, simNo, subtype, fullN, whichSource, sim=files[m], FullProp.notRedun, OGpath)
        compiled<-rbind(compiled, newLine)
        print(newLine)
      }
    }
  }
}
return(compiled)

}

# compiledEmp<-compile.secr.results()
# compiledSim<-compile.secr.results("SimData")
# compiled<-rbind(compiledEmp, compiledSim)

# 
# write.csv(compiledEmp, "Checks/CompiledEmp2019.csv")
# write.csv(compiledSim, "Checks/CompiledSim2019.csv")


discrep<-NULL # need a data frame describing the discrepancy between full sims and true density
avgEmpFull <-  filter(compiled, whichSource=='EmpData', subtype == 'Full') %>% summarise()
skipped<-0
noFull <- data.frame()
OGskipped<-character()
for (j in unique(compiled$OGpath)){
  subSims <- filter(compiled, OGpath == j, subtype != 'Full', !is.na(g0))
  if (nrow(subSims)<1) {
    skipped<-skipped+1
    print(paste0(j, " skipped: n skipped = ", skipped))
    OGskipped <- c(OGskipped, j)
    next()}
  fullSim <- filter(compiled, OGpath == j, subtype == 'Full')
  if (nrow(fullSim)<1){noFull <- bind_rows(subSims, noFull); next()}
  
  for (i in 1:(nrow(subSims) + 1)){
    if(i == (nrow(subSims) + 1)){
      newLine <- fullSim
    } else{
      newLine <- subSims[i,]
    }
    newLine$DenhatBiasFullPct <- newLine$Denhat / fullSim$Denhat
    newLine$DenhatBiasFull <- newLine$Denhat - fullSim$Denhat
    newLine$DenhatBiasTrue <- newLine$Denhat - trueDensity #Universal secr density is animals per hectare
    newLine$DenhatBiasTruePct <- newLine$Denhat / trueDensity #Universal secr density is animals per hectare

    discrep <- bind_rows(discrep, newLine)
    print(newLine)
  }
}
#Skip we need the discrepancy for each individual simulation
# #Todo: fix for two diff models
# # discrep<-NULL
# # for (k in c("SimData", "EmpData")){
# #   for (h in c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8")){
# #     for (i in c(250,550,850)){
# #       for (t in  c("g0 tilde bk", "g0 tilde 1", "g0 tilde bk + t", "g0 tilde t")){
# #           j<-10001
# #           for (u in 1:(nrow(compiled))){
# #               nextSim<-filter(compiled, trial==h, simNo==j, whichSource==k, model==t, size==i)
# #               if(nrow(nextSim)==0){break}
# #               newLine<-nextSim
# #               if(nrow(filter(nextSim, subtype=="Full"))>0){newLine$DenhatBiasFull<-nextSim[,4]-filter(nextSim, subtype=="Full")[,"Denhat"]} #Universal secr density is animals per hectare
# #               #newLine$DenhatBiasTrue<-nextSim[,4]-trueDensity #Universal secr density is animals per hectare
# #               
# # 
# #               discrep<-bind_rows(discrep, newLine)
# #               j<-j+1
# #           }
# #       }
# #     }
# #   }
# # }

# write.csv(discrep, 'Checks/discrep_newSims.csv')

## cleanup
# colnames(discrep)[13]<-"SubsamplingType"

compiledEmp <- read.csv("Checks/CompiledEmp2019.csv")
compiledSim <- read.csv("Checks/CompiledSim2019.csv")
compiled <- rbind(compiledSim, compiledEmp)
discrep <- read.csv('Checks/discrep_newSims.csv')

compiledA <- filter(compiled, subtype == "Full")
compiledB <- filter(compiled, subtype != "Full")
 compiledA$sizeRep <- 'Full'
 compiledB$sizeRep <- compiledB$size
 
compiledA$SubProp.notRedun<-compiledA$FullProp.notRedun
compiled<-rbind(compiledA, compiledB)
compiled$subtype <- str_replace(compiled$subtype, "Spread.one", "SPR")
compiled$subtype <- str_replace(compiled$subtype, "SimpleRandom", "SRS")
compiled$SubsamplingType <- compiled$subtype

compiled <- mutate(compiled, OGpath = paste0(whichSource, "_", trial, "_", model, "_", simNo))
compiled <- mutate(compiled, SubProp.notRedun = SubN.notRedun / size)
# 
discrep$SubsamplingType <- str_replace(discrep$SubsamplingType, "Spread.one", "SPR")
discrep$SubsamplingType <- str_replace(discrep$SubsamplingType, "SimpleRandom", "SRS")

discrepA <- discrep %>% filter(SubsamplingType == 'Full')
discrepB <- discrep %>% filter(SubsamplingType != 'Full')
discrepA$sizeRep <- 'Full'
discrepB$sizeRep <- as.character(discrepB$sizeRep)
discrep<- bind_rows(discrepA, discrepB)

discrepA <- discrep %>% filter(whichSource == 'EmpData')
discrepB <- discrep %>% filter(whichSource != 'EmpData')
discrepA$trial <- 'Empirical'
discrep<- bind_rows(discrepA, discrepB)

discrep <- discrep %>% mutate(model = str_replace(string = model, "tilde", "\\~"))


# #Money plots


failedSize <- compiled %>% filter(is.na(g0))
compiled <- compiled %>% filter(!is.na(g0))

# ggplot(data=filter(compiled, trial %in% c("t4","t5","t6","t7"), subtype!="Full"))+
#   scale_fill_grey(name="Subsampling\nType") +
#   geom_boxplot(aes(x=trial, y=SubProp.notRedun, fill=subtype)) +
#   ylab("Proportion of \nnon-redundant Samples") +
#   xlab("Scenario") +
#   ggtitle("Sample Redundancy vs \nSubsampling Type (Simulated)") +
#   theme_grey(14) +
#   facet_wrap(~size) + ggsave("Figures/PropNonRedun_Sim.png", height=5, width = 7.5)

PropNotRedun_Sim <- ggplot(data=filter(compiled, trial %in% c("t4","t5","t6","t7"), SubProp.notRedun>0))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=SubProp.notRedun, fill=subtype)) +
  ylab("Proportion of \nnon-redundant Samples") +
  xlab("Scenario") +
  ggtitle("Sample Redundancy vs \nSubsampling Type (Simulated)") +
  scale_fill_manual(values= c( 'gray', 'dodgerblue1', 'dodgerblue4')) +
  theme_grey(14) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~sizeRep)+ ggsave("Figures/Color_PropNonRedun_Sim_withFull.png", height=5, width = 7.5)

ggplot(data=filter(compiled, trial %in% c("t4","t5","t6","t7"), SubProp.notRedun>0), 
       aes(x=trial, y=SubProp.notRedun, fill=subtype))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_violin(scale = 'width') +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9))+
  ylab("Proportion of \nnon-redundant Samples") +
  xlab("Scenario") +
  ggtitle("Sample Redundancy vs \nSubsampling Type (Simulated)") +
  scale_fill_manual(values= c( 'gray', 'dodgerblue1', 'dodgerblue4')) +
  theme_grey(14) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~sizeRep, ncol = 4)+ ggsave("Figures/violin_Color_PropNonRedun_Sim_withFull.png", height=5, width = 7.5)

d1 = filter(discrep, whichSource == "EmpData")
d2 = filter(discrep, whichSource == "SimData")
d1$trial = "Empirical"
discrep = rbind(d1,d2)

ggplot(data=filter(discrep, whichSource=="SimData", size=="250", model == "g0 ~ 1" | model == "g0 ~ bk", SubsamplingType!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType)) +
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n (Simulated, n=250)')) +
  # theme_grey(14) +  
  scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  xlab("Trial") +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  facet_wrap(~model, ncol=4) + ggsave("Figures/Dsub_Dfull_simmed.png", height=5, width = 7.5)

meanTrials <- discrep %>% group_by(trial, SubsamplingType, model) %>%
  summarise(meanT = mean(DenhatBiasFullPct)) %>%
  filter(model == "g0 tilde 1" | model == "g0 tilde bk")

ggplot(data=filter(discrep, size=="250", model == "g0 ~ 1" | model == "g0 ~ bk", SubsamplingType!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4),
       aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_violin() +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9))+
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n (Simulated, n=250)')) +
  # theme_grey(14) +
  scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  xlab("Trial") +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~str_replace(string = model, "tilde", "\\~"), ncol=4) + ggsave("Figures/Violin_Dsub_Dfull_all.png", height=5, width = 7.5)
# 
# ggplot(data=filter(discrep, size=="250", model == "g0 tilde 1" | model == "g0 tilde bk", SubsamplingType!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4))+
#   scale_fill_grey(name="Subsampling\nType") +
#   geom_violin(aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType), addMean=TRUE) +
#   ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
#   ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n (Simulated, n=250)')) +
#   # theme_grey(14) +
#   scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
#   xlab("Trial") +
#   geom_hline(aes(yintercept=1), linetype='dashed')+
#   theme(legend.position="bottom", legend.title=element_blank()) +
#   facet_wrap(~str_replace(string = model, "tilde", "\\~"), ncol=4) + ggsave("Figures/Violin_Dsub_Dfull_both.png", height=5, width = 7.5)

ggplot(data=filter(discrep, size=="250", model == "g0 ~ 1" | model == "g0 ~ bk", SubsamplingType!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType)) +
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n  (Simulated, n=250)')) +
  theme_grey(14) +  
  xlab("Scenario") +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  facet_wrap(~model, ncol=4) + ggsave("Figures/Dsub_Dfull_simmed.png", height=5, width = 7.5)

# 
# ggplot(data=filter(discrep, whichSource!='EmpData', model == "g0 ~ 1" | model == "g0 ~ bk", subtype!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4))+
#   geom_boxplot(aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType)) +
#   ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
#   ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n  (Simulated)')) +
#   theme_grey(14) +  
#   xlab("Scenario") +
#   geom_hline(aes(yintercept=1), linetype='dashed')+
#   scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
#   theme(legend.position="bottom", legend.title=element_blank()) +
#   facet_grid(sizeRep~model) + ggsave("Figures/color_Dsub_Dfull_simmed_allSize.png", height=7.5, width = 7.5)

ggplot(data=filter(discrep, whichSource!='EmpData', model == "g0 ~ 1" | model == "g0 ~ bk", SubsamplingType!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4),
       aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType))+
  geom_violin() +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9))+
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n  (Simulated)')) +
  theme_grey(14) +  
  xlab("Scenario") +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_grid(size~model) + ggsave("Figures/violin_color_Dsub_Dfull_simmed_allSize.png", height=7.5, width = 7.5)

ggplot(data=filter(discrep, whichSource=="EmpData", 
                   SubsamplingType!="Full"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x = as.character(size), y=DenhatBiasFullPct, fill=SubsamplingType)) +
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n (Empirical)')) +
  theme_grey(14) +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~model,ncol=4) +
  xlab("Size of Subsample") + ggsave("Figures/Color_Dsub_Dfull_emp.png", height=5, width = 7.5)

ggplot(data=filter(discrep, whichSource=="EmpData", 
                   SubsamplingType!="Full"), aes(x = as.character(size), y=DenhatBiasFullPct, fill=SubsamplingType))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_violin()+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9)) +
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n (Empirical)')) +
  theme_grey(14) +
  geom_hline(aes(yintercept=1), linetype='dashed')+
  scale_fill_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~model,ncol=4) +
  xlab("Size of Subsample") + ggsave("Figures/violin_Color_Dsub_Dfull_emp.png", height=5, width = 7.5)




ggplot(data=filter(discrep, whichSource=="SimData", model == "g0 ~ 1" | model == "g0 ~ bk", trial %in% c("t1", "t3","t4", "t5", "t8"),
                   DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4, g0 <2.5), aes(x=trial, y=g0, fill=SubsamplingType))+
  geom_boxplot()+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9))+
  ylab(TeX('$\\hat{g}_{0}$')) +
  xlab("Scenario") +
  ggtitle(TeX("$\\g_0$ vs Subsampling Type (Simulated)")) +
  theme_grey(14) +
  facet_grid(size~model) +
  scale_fill_manual(values= c('gray', 'dodgerblue1', 'dodgerblue4')) +
  theme_grey(14) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(aes(linetype=5), yintercept=.5) + ggsave("Figures/color_g0_simmed.png", height=7.55, width = 7.5)


ggplot(data=filter(discrep, whichSource=="SimData", model == "g0 ~ 1" | model == "g0 ~ bk", trial %in% c("t1", "t3","t4", "t5", "t8"),
                   DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4, g0 <2.5), aes(x=trial, y=g0, fill=SubsamplingType))+
  geom_violin()+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9)) +
  ylab(TeX('$\\hat{g}_{0}$')) +
  xlab("Scenario") +
  ggtitle(TeX("$\\g_0$ vs Subsampling Type (Simulated)")) +
  theme_grey(14) +
  facet_grid(sizeRep~model) +
  scale_fill_manual(values= c('gray', 'dodgerblue1', 'dodgerblue4')) +
  theme_grey(14) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(aes(linetype=5), yintercept=.5) + ggsave("Figures/violin_color_g0_simmed.png", height=7.55, width = 7.5)


#from john - Modify Fig 5 so that it displays the mean (D^/D true) + a confidence interval for D^/Dtrue - rather than all of the individual results.  This will allow us to tell where there is strong evidence for bias; by contrast, if CIs overlap 1, then differences might just be explained by sampling error/small numbers of sims.  To get the CI, use  mean (D^/Dtrue)  +/- 1.96*sd(D^/Dtrue)/sqrt(nsims).  
# Except the below is wrong because it is comparing to full D hat rather than d true
discrep2 <- discrep %>% 
  filter(whichSource!='EmpData', !is.na(DenhatBiasTruePct))%>%
  group_by(SubsamplingType, trial, model, sizeRep) %>%
  summarise(meanDenhatFullPct = mean(DenhatBiasFullPct), CI_denhatBias = 1.96*sd(DenhatBiasFullPct)/sqrt(n()),
            meanDenhatTruePct = mean(DenhatBiasTruePct), CI_denhatTrueBias = 1.96*sd(DenhatBiasTruePct)/sqrt(n()), N = n())

ggplot(data=filter(discrep2, SubsamplingType!='Full', trial!="Empirical", size==250, model %in% c("g0 ~ 1", "g0 ~ bk")), aes(color = SubsamplingType,x=trial, group = SubsamplingType))+
  geom_point(aes(y=meanDenhatFullPct, shape=SubsamplingType), position=position_dodge(.4), size=3, fill="white", stat = 'identity') +
  geom_errorbar(aes(ymin=meanDenhatFullPct - CI_denhatBias, ymax=meanDenhatFullPct + CI_denhatBias), position=position_dodge(.4), colour="black", width = .4) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ (Simulated, n=250)')) +
  ylab("Derived density of subsample / \nDerived density of full sample (95% CI)") +
  xlab("Scenario") +
  theme_grey(14) +
  scale_color_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  facet_wrap(~model, ncol=2)+
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=1, linetype =2)  + ggsave("Figures\\pointsCI_Dsub_Dfull_simmed.png", height=5, width = 7.5)

ggplot(data=filter(discrep2, SubsamplingType!='Full',trial!="Empirical", model %in% c("g0 ~ 1", "g0 ~ bk")), aes(color = SubsamplingType,x=trial, group = SubsamplingType))+
  geom_point(aes(y=meanDenhatFullPct, shape=SubsamplingType), position=position_dodge(.4), size=3, fill="white", stat = 'identity') +
  geom_errorbar(aes(ymin=meanDenhatFullPct - CI_denhatBias, ymax=meanDenhatFullPct + CI_denhatBias), position=position_dodge(.4), colour="black", width = .4) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ (Simulated, n=250)')) +
  ylab("Derived density of subsample / \nDerived density of full sample (95% CI)") +
  xlab("Scenario") +
  theme_grey(14) +
  scale_color_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  facet_grid(size~model)+
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=1, linetype =2)  + ggsave("Figures\\pointsCI_ALL_Dsub_DFull_simmed.png", height=5, width = 7.5)

# now for what he ACTUALLY wanted.... dsub over true D rather than full Dhat estimate
ggplot(data=filter(discrep2, trial!="Empirical", sizeRep %in% c("250", "Full"), model %in% c("g0 ~ 1", "g0 ~ bk")), aes(color = SubsamplingType,x=trial, group = SubsamplingType))+
  geom_point(aes(y=meanDenhatTruePct, shape=SubsamplingType), position=position_dodge(.4), size=3, fill="white", stat = 'identity') +
  geom_errorbar(aes(ymin=meanDenhatTruePct - CI_denhatTrueBias, ymax=meanDenhatTruePct + CI_denhatTrueBias), position=position_dodge(.4), colour="black", width = .4) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $D$ (Simulated, n=250)')) +
  ylab("Derived density / True density (95% CI)") +
  xlab("Scenario") +
  theme_grey(14) +
  scale_color_manual(values= c('seagreen','dodgerblue1', 'dodgerblue4')) +
  facet_wrap(~model, ncol=2)+
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=1, linetype =2)  + ggsave("Figures\\pointsCI_Dsub_D_simmed.png", height=5, width = 7.5)

# 
ggplot(data=filter(discrep2, trial!="Empirical", model %in% c("g0 ~ 1", "g0 ~ bk")), aes(color = SubsamplingType,x=trial, group = SubsamplingType))+
  geom_point(aes(y=meanDenhatTruePct, shape=SubsamplingType), position=position_dodge(.4), size=3, fill="white", stat = 'identity') +
  geom_errorbar(aes(ymin=meanDenhatTruePct - CI_denhatTrueBias, ymax=meanDenhatTruePct + CI_denhatTrueBias), position=position_dodge(.4), colour="black", width = .4) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $D$ (Simulated, n=250)')) +
  ylab("Derived density / True density (95% CI)") +
  xlab("Scenario") +
  theme_grey(14) +
  scale_color_manual(values= c('seagreen','dodgerblue1', 'dodgerblue4')) +
  facet_wrap(sizeRep~model, ncol=2)+
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=1, linetype =2)  + ggsave("Figures\\pointsCI_ALL_Dsub_D_simmed.png", height=5, width = 7.5)

ggplot(data=filter(discrep2, trial!="Empirical", SubsamplingType!='Full', size==250, model %in% c("g0 ~ 1", "g0 ~ bk")), aes(color = SubsamplingType,x=trial, group = SubsamplingType))+
  geom_point(aes(y=meanDenhatTruePct, shape=SubsamplingType), position=position_dodge(.4), size=3, fill="white", stat = 'identity') +
  geom_errorbar(aes(ymin=meanDenhatTruePct - CI_denhatTrueBias, ymax=meanDenhatTruePct + CI_denhatTrueBias), position=position_dodge(.4), colour="black", width = .4) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $D$ (Simulated, n=250)')) +
  ylab("Derived density / True density (95% CI)") +
  xlab("Scenario") +
  theme_grey(14) +
  scale_color_manual(values= c('dodgerblue1', 'dodgerblue4')) +
  facet_wrap(~model, ncol=2)+
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=1, linetype =2)  + ggsave("Figures\\pointsCI_Dsub_D_noFull_simmed.png", height=5, width = 7.5)


#Stale money plots

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

ggplot(data=filter(discrep, whichSource=="SimData"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=plogis(g0), fill=SubsamplingType)) +
  ylab("g0") +
  ggtitle("g0 vs Subsampling Type (Simulated)") +
  theme_grey(14) +
  facet_wrap(size~model, ncol=4) +
  geom_hline(aes(linetype=5), yintercept=plogis(.5))

  discrep = mutate(discrep, model = str_replace(model, "tilde", "\\~"))

ggplot(data=filter(discrep, whichSource=="SimData",trial=="t2" | trial=="t6" | trial =="t7", model == "g0 ~ bk"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=plogis(g0.bTRUE + g0), fill=SubsamplingType)) +
  ylab(TeX("g_{bk}")) +
  xlab("Trial") + 
  ggtitle(TeX("g_{bk} vs Subsampling Type (Simulated, Behavioral Effect Present)")) +
  theme_grey(14) +
  geom_hline(aes(linetype=5), yintercept=plogis(1.5))

ggplot(data=filter(discrep, whichSource=="SimData",trial=="t2" | trial=="t6" | trial =="t7", model == "g0 ~ bk"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=plogis(g0), fill=SubsamplingType)) +
  ylab(TeX("g_0")) +
  xlab("Trial") + 
  ggtitle(TeX("g_0 vs Subsampling Type (Simulated, Behavioral Effect Present)")) +
  theme_grey(14) +
  geom_hline(aes(linetype=5), yintercept=plogis(.5))

tmp = filter(discrep, whichSource=="SimData",trial=="t2" | trial=="t6" | trial =="t7", model == "g0 ~ bk")
tmp$disp = 'g[0]'
tmp$yval = plogis(tmp$g0)

tmp2 = tmp
tmp2$disp = 'g[bk]'
tmp2$yval = plogis(tmp2$g0.bTRUE + tmp2$g0)
tmp = rbind(tmp, tmp2)
Zline = data.frame(disp = c("g[0]", "g[bk]"), Z = c(plogis(.5), plogis(1.5)))

ggplot(data=tmp)+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=yval, fill=SubsamplingType)) +
  ylab(TeX("g_0 or g_{bk}")) +
  xlab("Scenario") + 
  ggtitle(TeX("g_0 and g_{bk} vs Subsampling Type (Simulated, Behavioral Effect Absent)")) +
  theme_grey(14) +
  facet_grid(size~disp, labeller = label_parsed) +
  geom_hline(data = Zline, aes(yintercept=Z), linetype=5)

tmp = filter(discrep, whichSource=="SimData",trial=="t2" | trial=="t6" | trial =="t7", model == "g0 ~ bk")
tmp$disp = 'g[0]'
tmp$yval = tmp$g0

tmp2 = tmp
tmp2$disp = 'g[bk]'
tmp2$yval = tmp2$g0.bTRUE
tmp = rbind(tmp, tmp2)
Zline = data.frame(disp = c("g[0]", "g[bk]"), Z = c(.5, 1))

ggplot(data=filter(tmp, yval<4))+
  geom_boxplot(aes(x=trial, y=yval, fill=SubsamplingType)) +
  ylab(TeX("g_0 or g_{bk}")) +
  xlab("Scenario") + 
  ggtitle(TeX("g_0 and g_{bk} vs Subsampling Type")) +
  theme_grey(14) +
  facet_grid(sizeRep~disp, labeller = label_parsed) +
  scale_fill_manual(values= c( 'gray', 'dodgerblue1', 'dodgerblue4')) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(data = Zline, aes(yintercept=Z), linetype=5) + ggsave("Figures/color_g0_gbk_notPlogis.png", height=7.5, width = 7.5)


ggplot(data=filter(tmp, yval<4), aes(x=trial, y=yval, fill=SubsamplingType))+
  geom_violin()+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point",
               size = .75,
               shape = 23,
               stroke = 1.5,
               color = 'white',
               position = position_dodge(.9)) +
  ylab(TeX("g_0 or g_{bk}")) +
  xlab("Scenario") + 
  ggtitle(TeX("g_0 and g_{bk} vs Subsampling Type")) +
  theme_grey(14) +
  facet_grid(sizeRep~disp, labeller = label_parsed) +
  scale_fill_manual(values= c( 'gray', 'dodgerblue1', 'dodgerblue4')) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(data = Zline, aes(yintercept=Z), linetype=5) + ggsave("Figures/violin_color_g0_gbk_notPlogis.png", height=7.5, width = 7.5)

ggplot(data=filter(discrep, whichSource=="SimData"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=plogis(g0), fill=SubsamplingType)) +
  ylab(TeX("g_0")) +
  xlab("Trial") + 
  ggtitle(TeX("g_0 vs Subsampling Type (Simulated, Behavioral Effect Absent)")) +
  theme_grey(14) +
  facet_wrap(size~model, ncol=4) +
  geom_hline(aes(linetype=5), yintercept=plogis(.5))

ggplot(data=filter(discrep, whichSource=="SimData", subtype!="Full", DenhatBiasFullPct > .5 & DenhatBiasFullPct < 1.4, model == "g0 ~ bk" | model == "g0 ~ 1"))+
  scale_fill_grey(name="Subsampling\nType") +
  geom_boxplot(aes(x=trial, y=DenhatBiasFullPct, fill=SubsamplingType)) +
  ylab(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$')) +
  ggtitle(TeX('$\\hat{D}_{Sub}$ / $\\hat{D}_{Full}$ \n vs Subsampling Type and Simulation Scenario (Simulated)')) +
  theme_grey(14) +
  xlab("Scenario") + 
  geom_hline(aes(yintercept=1), linetype='dashed')+
  facet_wrap(size~str_replace(string = model, "tilde", "\\~"), ncol=2)

#Looking at the actual bk values from the original data for text
bkt <- readRDS('C:/Users/Nick/Documents/GitHub/spatialMR-master/Sim_SpatialMR_SuggestBuffer_EqualMask/EmpData/t1/g0 tilde bk + t/850/Full10015.rds')
bk  <- readRDS('C:/Users/Nick/Documents/GitHub/spatialMR-master/Sim_SpatialMR_SuggestBuffer_EqualMask/EmpData/t1/g0 tilde bk/850/Full10015.rds')

bkt$fit$par
bkt$betanames

bk$fit$par
bk$betanames



