#Load libraries
library(gridExtra)
library(ggplot2)
library(viridis)
library(mgcv)
library(reshape2)
library(ggrepel)
library(cowplot)
library(plotly)
library(dplyr)
library(Hmisc)

###############################
# ggplots (bar and dot) of yield and recruitment timescales
# meanage from influence function (YPR) 
###############################

# read in datasets; results from matlab convolution code
yield_ts_recruitment_ts=read.csv("species_halfr_meanage_af.csv")

# derived quantities
yield_ts_recruitment_ts$diff=yield_ts_recruitment_ts$Y95-yield_ts_recruitment_ts$B95
yield_ts_recruitment_ts$one_M=1/yield_ts_recruitment_ts$M
yield_ts_recruitment_ts$Sp=as.factor(yield_ts_recruitment_ts$Species)
yield_ts_recruitment_ts$FM=yield_ts_recruitment_ts$F/yield_ts_recruitment_ts$M

# mold dataset for ggplot
ts=subset(yield_ts_recruitment_ts,select=-c(F, M, Af, MeanAge, diff, FM, Sp, one_M))
molten.ts=melt(ts,id=c("Species"))
molt.ts=inner_join(molten.ts, yield_ts_recruitment_ts, by = "Species")

#----------------------------------------------------------------------------------
# Plot of ratio of convolution 95% timescale and biomass timescale in a bar plot
quartz()
Fig.5<-ggplot(molt.c,aes(x = reorder(Species, M), y=value, fill=as.factor(M), position=variable),color="black") +
  geom_bar(stat="identity", position = position_dodge(width = 0.5),color="black",show.legend = FALSE)+
  #geom_text(aes(x=Species,y=value+1.5,label=value))+
  scale_fill_viridis(discrete = T,option="C")+
  theme_bw()+
  labs(x= "", y = "Years")+
  ylim(0, 65) +
  coord_flip()
Fig.5

#----------------------------------------------------------------------------------
# Plot of difference between yield and recruitment timescales vs Meanage from YPR
quartz()
Fig.6<-ggplot(b, aes(x=MeanAge,y=diff)) +
  geom_point(aes(col=Species,shape=Species), size=4)  +
  scale_shape_manual(values=1:nlevels(b$Species))+
  theme_bw()+
  xlab('Mean Age')+
  ylab('Difference between 95% yield time & 95% biomass time')
Fig.6
