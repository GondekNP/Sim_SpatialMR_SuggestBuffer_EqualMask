#This function will take on either an individual secr model or a single trial x subsampling combo and overlay the detection curves.
rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)

empBetas <- readRDS("Checks\\AllBetasEmpData.rds")
simBetas <- readRDS("Checks\\AllBetasSimData.rds")
betaFrame <- data.frame()
for (j in 1:length(empBetas)){
  
  if (!is.null(empBetas[[j]][[2]])){
    
      betas = unlist(empBetas[[j]][2])
      names(betas) = unlist(empBetas[[j]][3])
      new = data.frame(t(as.matrix(betas)))
  
      new$path = unlist(empBetas[[j]][1])[1]
      new$trial = unlist(empBetas[[j]][1])[2]
      new$modEval = unlist(empBetas[[j]][1])[3]
      new$size = unlist(empBetas[[j]][1])[4]
      new$source = "Emp"
      
      betaFrame = bind_rows(betaFrame, new)
    }
}

for (j in 1:length(simBetas)){
  
  if (!is.null(simBetas[[j]][[2]])){
    
    betas = unlist(simBetas[[j]][2])
    names(betas) = unlist(simBetas[[j]][3])
    new = data.frame(t(as.matrix(betas)))
    
    new$path = unlist(simBetas[[j]][1])[1]
    new$trial = unlist(simBetas[[j]][1])[2]
    new$modEval = unlist(simBetas[[j]][1])[3]
    new$size = unlist(simBetas[[j]][1])[4]
    new$source = "Sim"
    
    betaFrame = bind_rows(betaFrame, new)
  }
}

for (k in 1:nrow(betaFrame)){
    if (str_detect(betaFrame$path[k], "Full")) {betaFrame$subtype[k] <- "Full"}
    if (str_detect(betaFrame$path[k], "Simple")) {betaFrame$subtype[k] <- "SRS"}
    if (str_detect(betaFrame$path[k], "Spread")) {betaFrame$subtype[k] <- "SPR"}
  }

plotDetectionCurves <- function(Betas, source1, model, trial1, sub){
   if (sub != "Full"){Betas <- filter(Betas, modEval == model, trial == trial1, source==source1, subtype == sub, size == size1)
   titleForm = paste0("Input Detection Curve vs Simulated Detection Curves\n (Trial ",trial1,", Model ", str_replace(model, "tilde", "~"), ", Type ", sub,", Size ", size1, ")")
   titleForm2 = paste0("Input Post-capture Detection Curve vs Simulated Post-capture Detection Curves\n (Trial ",trial1,", Model ", str_replace(model, "tilde", "~"), ", Type ", sub,", Size ", size1, ")")}
  
  if (sub == "Full"){Betas <- filter(Betas, modEval == model, trial == trial1, source==source1, subtype == sub)
  titleForm = paste0("Input Detection Curve vs Simulated Detection Curves\n (Trial ",trial1,", Model ", str_replace(model, "tilde", "~"), ", Type ", sub, ")")
  titleForm2 = paste0("Input Post-capture Detection Curve vs Simulated Post-capture Detection Curves\n (Trial ",trial1,", Model ", str_replace(model, "tilde", "~"), ", Type ", sub,")")}

   dists <- seq(0,3000, length=100)
   sig_hat <- sqrt((1500*1500)/pi)
   
   g0_hat <- plogis(.5)
   DP_hat <- g0_hat * exp((-(dists^2))/(2*sig_hat^2)) #true detection curve for t1
   
   g0bk_hat<- plogis(.5 + 1)
   DPbk_hat <- g0bk_hat * exp((-(dists^2))/(2*sig_hat^2)) #true detection curve for t1
   
   combined<-data.frame()
       for (j in 1:nrow(Betas)){
         g0_new <- plogis(Betas[j,"g0"])
         sig_new <- exp(Betas[j, "sigma"])
         DP_new <- g0_new * exp((-(dists^2))/(2*sig_new^2))
         newLine <- data.frame(Type = "Simulated", dummy = j, Sim = Betas[j,"path"], X=dists, Y = DP_new, bk = FALSE)
         combined<-bind_rows(combined, newLine)
         
         g0_new <- plogis(Betas[j,"g0"] + Betas[j,"g0.bkTRUE"])
         DP_new <- g0_new * exp((-(dists^2))/(2*sig_new^2))
         newLine <- data.frame(Type = "Simulated", dummy = j, Sim = Betas[j,"path"], X=dists, Y = DP_new, bk = TRUE)
         combined<-bind_rows(combined, newLine)
       }
   
    truLine = data.frame(Type = "Input", Sim = "HAT", X=dists, Y = DP_hat, dummy = 10000000, bk = FALSE)
    truLine_bk = data.frame(Type = "Input", Sim = "HAT", X=dists, Y = DPbk_hat, dummy = 10000001, bk = TRUE)
    
    combined<-bind_rows(truLine, truLine_bk, combined)

        int_curve <- ggplot(data=filter(combined, bk == FALSE), aes(x=X,y=Y, color = Type, group=dummy))+
                geom_line(size=1) +
                #geom_line(data=truLine, aes(x=X, y=Y), colour = "red", size=2) + 
                labs(title=titleForm,
                     x="Meters from AC", y="Detection Probability")+
                scale_color_manual(values = c('Simulated' = 'gray50', 'Input' = 'red'))
        
        bk_curve <- ggplot(data=filter(combined, bk == TRUE), aes(x=X,y=Y, color = Type, group=dummy))+
          geom_line(size=1) +
          #geom_line(data=truLine, aes(x=X, y=Y), colour = "red", size=2) + 
          labs(title=titleForm2,
               x="Meters from AC", y="Detection Probability after Capture")+
          scale_color_manual(values = c('Simulated' = 'gray50', 'Input' = 'red'))
        
        return(list(int_curve, bk_curve))
}

