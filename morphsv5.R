#Analysis of categorical moth coloration with PRISM climate normals, 1981-2010
#Written in in NYC 8-16-19
#Adjusted November 2020 for writing final paper=
#Cleaned for archive with finished paper in Sept 2022.

rm(list=ls())
library(lme4)
library(ggplot2)
library(reshape)
library(glmmTMB)
library(emmeans)
library(gridExtra)
library(betareg)
library(gstat)
library(tidyverse)
library(raster)
library(maptools)
library(sp)
library(AICcmodavg)
library(ggeffects)
library(ggnewscale)
library(scatterpie)
library(rgdal)
library(ggrepel)
setwd('~/Documents/Research/dissertation/coloration/wingcolor')

morphs<-read.csv('morphs3.csv')
str(morphs)

morphs$proplight<-morphs$light/(morphs$light+morphs$dark)
morphs$proplight[1:10]<-morphs$proplight[1:10]-0.0000001
morphs$proplight[21]<-morphs$proplight[21]-0.0000001
morphs$proplight[11]<-morphs$proplight[11]+0.0000001
morphs$proplight
morphs$weight<-(morphs$light+morphs$dark)
morphs$Site<-morphs$site

morphsp<-read.csv('mothcolorprop.csv')
str(morphsp)

morphs2<-left_join(morphsp,morphs,by='Site')
str(morphs2)

morphs$pos <- numFactor(morphs$long,morphs$lat)
morphs$group<-1
morphs$group<-as.factor(morphs$group)
morphs$x <- 110.*(morphs$lat-42)
morphs$y <- 110.*(morphs$long-(-125))*cos(morphs$long/(360/(2*pi)))
morphs2$pos <- numFactor(morphs2$long,morphs2$lat)
morphs2$group<-1
morphs2$group<-as.factor(morphs2$group)
morphs2$x <- 110.*(morphs2$lat-42)
morphs2$y <- 110.*(morphs2$long-(-125))*cos(morphs2$long/(360/(2*pi)))
morphs2<-morphs2[!is.na(morphs2$melanin.prop),]
morphs2$vpdmean<-(morphs2$vpdmax_jja+morphs2$vpdmin_jja)/2

ggplot(morphs2,aes(x=melanin.prop))+geom_density()+theme_classic()+xlim(.3,1)+xlab("Melanin Proportion")+ylab('Density')
str(morphs2)

ghi <- shapefile(
  "~/Documents/Research/dissertation/coloration/wingcolor/maps/l48_ghi_10km.shp")
str(ghi)


pts <- SpatialPoints(cbind(morphs2$lat,morphs2$long), 
                     proj4string = CRS(proj4string(ghi)))

ghi1<-over(pts, ghi)
morphs2$ghijja<-rowMeans(ghi1[,10:12])
morphs2$ghijj<-rowMeans(ghi1[,10:11])
morphs2$ghija<-rowMeans(ghi1[,11:12])
morphs2$ghi<-morphs2$ghijj

morphs2$ghi[morphs2$site=='LAM']<-morphs2$ghija[morphs2$site=='LAM']
morphs2$ghi[morphs2$site=='WAC']<-morphs2$ghija[morphs2$site=='WAC']
morphs2$ghi[morphs2$site=='WIL']<-morphs2$ghija[morphs2$site=='WIL']

morphs2$stmean_jja<-scale(morphs2$tmean_jja)
morphs2$stmin_jja<-scale(morphs2$tmin_jja)
morphs2$stmax_jja<-scale(morphs2$tmax_jja)
morphs2$sppt_jja<-scale(morphs2$ppt_jja)
morphs2$sghi_jja<-scale(morphs2$ghijja)
morphs2$svpdmean<-scale(morphs2$vpdmean)


#Models
bglmm1sp<- glmmTMB(melanin.prop ~ tmean_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm1sp)

bglmm2sp<- glmmTMB(melanin.prop ~ tmin_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm2sp)

bglmm3sp<- glmmTMB(melanin.prop ~ tmax_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm3sp)


bglmm4sp<- glmmTMB(melanin.prop ~ ppt_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm4sp)

bglmm5sp<- glmmTMB(melanin.prop ~ log(ghijja) +gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm5sp)


bglmm6sp<- glmmTMB(melanin.prop ~ ghijja +gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm6sp)

bglmm7sp<- glmmTMB(melanin.prop ~ vpdmax_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm7sp)

bglmm8sp<- glmmTMB(melanin.prop ~ vpdmin_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm8sp)

bglmm9sp<- glmmTMB(melanin.prop ~ tdmean_jja + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm9sp)

bglmm10sp<- glmmTMB(melanin.prop ~ vpdmean + gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm10sp)

bglmm11sp<- glmmTMB(melanin.prop ~ vpdmean +tmean_jja+ gau(pos + 0|group)+(1|site),family=beta_family(), data=morphs2)
summary(bglmm11sp)


aictab(list(bglmm1sp,bglmm2sp,bglmm4sp,bglmm5sp,bglmm7sp,bglmm10sp,bglmm11sp))


mydf <- ggpredict(bglmm5sp, terms = c("ghijja"))
str(mydf)

mydf$ghijja<-mydf$x
mydf$melanin.prop<-mydf$predicted
mydfl<-mydf
mydfl$melanin.prop<-mydf$conf.low
mydfh<-mydf
mydfh$melanin.prop<-mydfh$conf.high
str(morphs2)
ghiplot<-ggplot(morphs2,aes(x=ghijja,y=melanin.prop))+scale_x_log10()+
  geom_density_2d_filled(alpha=0.7)+geom_point()+geom_line(data=mydf, size=1.5)+geom_line(data=mydfl, lty=2)+geom_line(data=mydfh, lty=2) + xlab("GHI") + ylab("Proportion melanized")+ theme_classic()+ggtitle('c')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20),legend.position = 'none')

ghiplot

mydf <- ggpredict(bglmm1sp, terms = c("tmean_jja"))
str(mydf)

mydf$tmean_jja<-mydf$x
mydf$melanin.prop<-mydf$predicted
mydfl<-mydf
mydfl$melanin.prop<-mydf$conf.low
mydfh<-mydf
mydfh$melanin.prop<-mydfh$conf.high

tmeanplot<-ggplot(morphs2,aes(x=tmean_jja,y=melanin.prop))+
  geom_point()+geom_line(data=mydf, size=1.5)+geom_line(data=mydfl, lty=2)+geom_line(data=mydfh, lty=2) + xlab("Mean Temperature") + ylab("Proportion melanized")+ theme_classic()+ggtitle('b')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20),legend.position = 'none')

tmeanplot

mydf <- ggpredict(bglmm2sp, terms = c("tmin_jja"))
str(mydf)

mydf$tmin_jja<-mydf$x
mydf$melanin.prop<-mydf$predicted
mydfl<-mydf
mydfl$melanin.prop<-mydf$conf.low
mydfh<-mydf
mydfh$melanin.prop<-mydfh$conf.high

tminplot<-ggplot(morphs2,aes(x=tmin_jja,y=melanin.prop))+
  geom_density_2d_filled(alpha=0.7)+geom_point()+geom_line(data=mydf, size=1.5)+geom_line(data=mydfl, lty=2)+geom_line(data=mydfh, lty=2) + xlab("Minimum Temperature") + ylab("Proportion melanized")+ theme_classic()+ggtitle('b')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20),legend.position = 'none')

tminplot

tmean<- raster('~/Documents/Research/dissertation/coloration/wingcolor/maps/PRISM_tmean_30yr_normal_4kmM2_06_bil/PRISM_tmean_30yr_normal_4kmM2_06_bil.bil')


mydf <- ggpredict(bglmm7sp, terms = c("vpdmax_jja"))
str(mydf)

mydf$vpdmax_jja<-mydf$x
mydf$melanin.prop<-mydf$predicted
mydfl<-mydf
mydfl$melanin.prop<-mydf$conf.low
mydfh<-mydf
mydfh$melanin.prop<-mydfh$conf.high

vpdmaxplot<-ggplot(morphs2,aes(x=vpdmax_jja,y=melanin.prop))+geom_point()+geom_line(data=mydf, size=1.5)+geom_line(data=mydfl, lty=2)+geom_line(data=mydfh, lty=2) + xlab("Max VPD") + ylab("Proportion melanized")+ theme_classic()+ggtitle('c')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20),legend.position = 'none')

vpdmaxplot

tmean<- raster('/Documents/Research/dissertation/coloration/wingcolor/maps/PRISM_tmean_30yr_normal_4kmM2_06_bil/PRISM_tmean_30yr_normal_4kmM2_06_bil.bil')



#Mean temp map
tmeandf <- as.data.frame(tmean, xy = TRUE)
str(tmeandf)
tmeanmap<-ggplot()+
  geom_raster(data = tmeandf , aes(x = x, y = y, fill =PRISM_tmean_30yr_normal_4kmM2_06_bil))+
  coord_quickmap()+ 
  scale_fill_gradient2(low="blue",mid="white",high="red",name = "Mean Temp",
                       midpoint = 15, na.value = "white",)+
  new_scale_fill()+
  geom_scatterpie(data=morphs,aes(y=long, x=lat, r = .3),cols = c("light", "dark") )+
  scale_fill_manual(name="Color Morphs",
                    breaks = c("light", "dark"),
                    labels = c("Light", "Dark"),
                    values = c("light" = "orange",
                               "dark" = "black")
  ) +
  theme_classic()+
  xlim(-125,-116)+
  ylim(36,49)+
  xlab("Longitude")+
  ylab("Latitude")+ggtitle('a')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

tmeanmap

library(Cairo)
options(device = Cairo::CairoWin)


morphs3<-morphs

str(morphs3)
#Data additions/corrections for mapping
morphs3[24,1]<-'FBC'
morphs3[25,1]<-'JFW'
morphs3[26,1]<-'NET'
morphs3[27,1]<-'LWC'
##Lat-lon reversed
morphs3$lat[24]<--123.06355884753197
morphs3$long[24]<- 48.50009015210723
morphs3$lat[25]<--123.240030
morphs3$long[25]<-44.604668
morphs3$lat[26]<--123.87647089928129
morphs3$long[26]<-46.12049592638012
morphs3$lat[27]<--123.8710149839507
morphs3$long[27]<-46.01642716844461


tmeanmaplabels<-ggplot()+
  geom_raster(data = tmeandf , aes(x = x, y = y, fill =PRISM_tmean_30yr_normal_4kmM2_06_bil))+
  coord_quickmap()+ 
  scale_fill_gradient2(low="blue",mid="white",high="red",name = "Mean Temp",
                       midpoint = 15, na.value = "white",) +
  geom_text_repel(data=morphs3,aes(y=long, x=lat, label = site),nudge_x = c(0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0,0,0),nudge_y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0),min.segment.length = 0.1)+geom_point(data=morphs3,aes(y=long, x=lat))+
  xlim(-125,-116)+
  ylim(36,49)+
  xlab("Longitude")+
  theme_classic()+
  ylab("Latitude")+ggtitle('a')

morphs3$site
tmeanmaplabels 

jpeg(file="labelledmap.jpeg",width=1000, height=600)
tmeanmaplabels 
dev.off()


rep(0,20)

lay <- rbind(c(1,2),
             c(1,3))
grid.arrange(tmeanmap,tmeanplot,vpdmaxplot,layout_matrix=lay, widths=c(1.5,1))

