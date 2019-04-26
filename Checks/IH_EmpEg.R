#Read in emp data and create a figure showing the heterogeneity in samples deposited
setwd("C:/Users/Nick/Documents/GitHub/Sim_SpatialMR_SuggestBuffer_EqualMask")
samps<-read.csv("samps1_uniq.csv")

library(dplyr)
library(mosaic)

unique(samps$ID) #43 Animals
tally_bearMR <- as.data.frame(tally(data=samps, ~INDuniqID))

ggplot(data = tally_bearMR, aes(x=Freq)) +
  geom_histogram(bins = 11, color = "black") +
  scale_x_continuous(breaks = 1:11)+
  xlab("Frequency") +
  ylab("Count")