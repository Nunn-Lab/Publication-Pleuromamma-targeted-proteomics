source('biostats.R')
library(dplyr)
library(vegan)
library(ggplot2)
library(reshape)
library(tidyverse)

#read in the data file of protein abundance
circ.dda<-read.csv('circadian DDA proteins.csv', row.names=1)

#read in file with differentially abundant proteins (output of qspec) + annotations
diff.ab<-read.csv('DDA DAPs RK with annot copy.csv')

#average day and night abundances
dda.day<-subset(circ.dda, select=c(X62, X63, X64, X66, X68, X69, X72, X75, X77))
dda.night<-subset(circ.dda, select=c(X65, X67, X70, X71, X74, X76))

day.avg<-rowMeans(dda.day)
night.avg<-rowMeans(dda.night)

avg.dda<-data.frame(day.avg, night.avg)
avg.dda$Protein<-rownames(avg.dda)

merge.avg<-merge(x=diff.ab, y=avg.dda, by='Protein', all.x=T)

merge.sub<-subset(merge.avg, select=c(Protein, Name, day.avg, night.avg, LogFoldChange))

lfc.bar<-
  ggplot(merge.sub, aes(x=reorder(Protein, LogFoldChange), y=LogFoldChange, fill=LogFoldChange)) +
  geom_bar(stat='identity', position=position_dodge()) +
  coord_flip() +
  theme_bw() +
  xlab('Protein Name') +
  ylab('Log Fold Change (day/night)') +
  scale_x_discrete(labels=prot.names) +
  theme(axis.text=element_text(size=6)) +
  scale_fill_stepsn(colours = c( 'blue', 'gold1'))


