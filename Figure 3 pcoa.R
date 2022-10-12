library(vegan)
library(ape)
library(ggplot2)
library(plyr)
library(dplyr)
library(corrplot)
library(naniar)

#pcoa with euclidean distance
#read in data file
srm.dat<-read.csv('file for PCoA.csv', row.names=1, na.strings='#N/A', header=T)

#subset the variables indicating time point
exp.srm<-subset(srm.dat, select=Night:TS9)

#replace NA with 0
srm.dat[is.na(srm.dat)]<-0

#get rid of 14 columns containing time point data
srm.dat3<-srm.dat[c(15:457)]

#PCoA in situ
#data in rows 16-33
srm.TS<-data.frame(srm.dat3[16:33,])
exp.TS<-data.frame(exp.srm[16:33,])
exp.TS<-subset(exp.TS, select=c(Night, Day, TS5, TS6, TS7, TS8, TS9))

TS.dist<-vegdist(srm.TS, method='euclidean')

TS.pcoa<-pcoa(TS.dist)

biplot.pcoa(TS.pcoa)
biplot.pcoa(TS.pcoa, exp.TS)

pcoa_TS<-cmdscale(TS.dist, k=2)
pcoa_TS<-as.data.frame(pcoa_TS)
names(pcoa_TS)[1:2]<-c('PC1', 'PC2')
pcoa_TS$timeofday<-c(rep('Night', 6), rep('Day', 8), rep('Night', 4))
pcoa_TS$timepoint<-c(5,5,5,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9)
gg.pcoaTS<-ggplot(pcoa_TS, aes(x=PC1, y=PC2, colour=as.factor(timepoint), shape=as.factor(timepoint)))

find_hull<- function(pcoa_TS) pcoa_TS[chull(pcoa_TS$PC1, pcoa_TS$PC2), ]
hulls<-ddply(pcoa_TS, 'timepoint', find_hull)

gg.pcoa2TS<-gg.pcoaTS + geom_point(aes(x=PC1, y=PC2, colour=as.factor(timepoint), shape=as.factor(timepoint), fill=as.factor(timepoint)),size=6) + 
scale_shape_manual(name='timepoint', labels=c('D1 22:00', 'D2 02:00', 'D2 09:00', 'D2 15:00', 'D2 22:00'), values=c(rep(21,4),22)) +
	scale_colour_manual(name='timepoint', labels=c('D1 22:00', 'D2 02:00', 'D2 09:00', 'D2 15:00', 'D2 22:00'), values=c('blue4', 'blue', 'gold2', 'gold', 'blue4')) +
	scale_fill_manual(name='timepoint', labels=c('D1 22:00', 'D2 02:00', 'D2 09:00', 'D2 15:00', 'D2 22:00'), values=c('blue4', 'blue', 'gold2', 'gold', 'blue4')) +
	theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.background = element_rect(fill = NA, colour = NA),
      text=element_text(size=20))+ 
      geom_polygon(data=hulls, aes(x=PC1, y=PC2, fill=as.factor(timepoint)), alpha=.2)


#PCoA incubation
#rows 1-15
srm.C<-data.frame(srm.dat3[1:15,])
exp.C<-data.frame(exp.srm[1:15,])
exp.C<-subset(exp.C, select=c(Night, Day, C1, C2, C3, C4, C6))

C.dist<-vegdist(srm.C, method='euclidean')

C.pcoa<-pcoa(C.dist)

biplot.pcoa(C.pcoa)
biplot.pcoa(C.pcoa, exp.C)

pcoa_C<-cmdscale(C.dist, k=2)
pcoa_C<-as.data.frame(pcoa_C)
names(pcoa_C)[1:2]<-c('PC1', 'PC2')
pcoa_C$timeofday<-c(rep('Day', 6), rep('Night', 6), rep('Day', 3))
pcoa_C$timepoint<-c(1,1,1,2,2,2,3,3,3,4,4,4,6,6,6)
gg.pcoaC<-ggplot(pcoa_C, aes(x=PC1, y=PC2, colour=as.factor(timepoint), shape=as.factor(timepoint)))

find_hull<- function(pcoa_C) pcoa_C[chull(pcoa_C$PC1, pcoa_C$PC2), ]
hulls<-ddply(pcoa_C, 'timepoint', find_hull)

gg.pcoa2C<-gg.pcoaC + geom_point(aes(x=PC1, y=PC2, colour=as.factor(timepoint), shape=as.factor(timepoint), fill=as.factor(timepoint)),size=6) +
scale_shape_manual(name='timepoint', labels=c('06:00', '12:00', '18:00', '24:00', '36:00'), values=c(rep(21,3), rep(22,2))) +
	scale_colour_manual(name='timepoint', labels=c('06:00', '12:00', '18:00', '24:00', '36:00'), values=c('gold2', 'gold', 'blue', 'blue4', 'gold')) +
	scale_fill_manual(name='timepoint', labels=c('06:00', '12:00', '18:00', '24:00', '36:00'), values=c('gold2', 'gold', 'blue', 'blue4', 'gold')) +
	theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.background = element_rect(fill = NA, colour = NA),
      text=element_text(size=20)) +
      geom_polygon(data=hulls, aes(x=PC1, y=PC2, fill=as.factor(timepoint)), alpha=.2)