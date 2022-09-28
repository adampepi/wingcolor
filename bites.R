#Analysis of plasticine clay moth predation
#Written in Vancouver 8-2-19 (only portland data atm)
#Cleaned sept 2022 for publication
rm(list=ls())

library(lme4)
library(ggplot2)
library(reshape)
library(glmmTMB)
library(emmeans)
library(car)

setwd('~/Documents/Research/dissertation/coloration/wingcolor')
bites<-read.csv('bites.csv')
str(bites)
bitespois<-read.csv("bitespois.csv")
str(bitespois)
bites2<-melt(bites[,1:6], id=c("id","morph","site")) 
str(bites2)

bites2p<-melt(bitespois[,1:6], id=c("id","morph","site")) 
str(bites2p)

bites$sure.bites

qplot(morph,value,facets=.~site,geom="violin",data=bites2)+geom_count()
qplot()

bites2POB<-bites2[bites2$site=="POB",]
bites2BMR<-bites2[bites2$site=="BMR",]


bitesPOB<-bites[bites$site=="POB",]
bitesBMR<-bites[bites$site=="BMR",]


glm2<-glm(sure.bites~morph*site,data=bites, family=binomial)
summary(glm2)
glm2r<-glm(sure.bites~morph+site,data=bites, family=binomial)
summary(glm2r)

anova(glm2r, glm2, test = "LRT")
anova(glm2r, glm2)

glm2bmr<-glm(sure.bites~morph,data=bitesBMR, family=binomial)
summary(glm2bmr)

glm2pob<-glm(sure.bites~morph,data=bitesPOB, family=binomial)
summary(glm2pob)

glm2.1<-glm(sure.bites~morph*site,data=bitespois, family=poisson)
summary(glm2.1)

lsmeans(glm2,~site*morph,transform='response')




ggplot(bites, aes(x = site, y = sure.bites,fill= morph)) + 
  stat_summary(fun.y="mean", geom= "bar", position="dodge")+
  stat_summary(fun.data = "mean_se", colour = "black", geom="errorbar", position="dodge")+
  xlab("Site") +
  ylab("Attack rate") + 
  scale_fill_manual(name="Morph", values=c("D"="grey30", "L"="orange"), 
                    labels=c("Dark","Light")) +
  labs(fill = "Morph")+theme_classic()+
  theme(text = element_text(size=30))

