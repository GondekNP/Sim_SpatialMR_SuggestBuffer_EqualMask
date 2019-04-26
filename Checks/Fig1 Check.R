##Figure 1 check - based on empirical data, how many traps did bears vist per period? and how many per trap?

library(tidyverse)
# setwd("C:/Users/ConservationMetrics/Documents/NPG/Sim_SpatialMR_SuggestBuffer_EqualMask")
setwd('C:/Users/Nick/Documents/GitHub/spatialMR-master/Sim_SpatialMR_SuggestBuffer_EqualMask')

fullsamps <- read.csv("samps1_uniq.csv")
head(fullsamps)

byPeriod <- fullsamps %>% group_by(Period, ID) %>% summarise(sites = length(unique(site))) # grouped by period and bear, how many unique sites?

byPeriodHist <- ggplot(byPeriod, aes(sites)) +
  geom_histogram(bins =15, color = 'black', fill = 'dodgerblue1') + 
  ylab('N sites visited in a period') +
  xlab('Frequency') +
  ggsave("Figures/NumSitesPerPeriod.png", height=5, width = 7.5)

byTraps <- fullsamps %>% group_by(Period, ID, site) %>% summarise(N = n()) # grouped by period and bear, how many unique sites?

byTrapsHist <- ggplot(byTraps, aes(N)) +
  geom_histogram(bins =11, color = 'black', fill = 'dodgerblue1') + 
  ylab('Number of bears') +
  xlab('Number of genotyped samples per individual at each site-session') +
  ggsave("Figures/NumClustersperSite.png", height=5, width = 7.5)

bySites <- fullsamps %>% group_by(Period, site) %>% summarise(N = n()) 

byTrapsHist <- ggplot(bySites, aes(N)) +
  geom_histogram(bins =11,color = 'black', fill = 'dodgerblue1') + 
  ylab('Number of site-sessions') +
  xlab('Number of samples genotyped') +
  ggsave("Figures/NumClustersperSite_forreal.png", height=5, width = 7.5)
