source('~/Documents/Multivariate Stats/biostats.R')
library(dplyr)
library(vegan)
library(ggplot2)
library(reshape)
library(ggthemes)

srm.dat<-read.csv('file for PCoA.csv', row.names=1, na.strings='#N/A', header=F)
srm.dat2<-read.csv('file for PCoA.csv', row.names=1, na.strings='#N/A', header=T)

colnames(srm.dat2)<-srm.dat[1,]

srm.dat2[is.na(srm.dat2)]<-0

#get rid of 14 columns with metadata
srm.dat3<-srm.dat2[c(15:457)]
srm.t<-data.frame(t(srm.dat3))

#average peptide transition abundances for each time point
C2T1<-cbind(srm.t$C2T1.1, srm.t$C2T1.2, srm.t$C2T1.3)
C2T2<-cbind(srm.t$C2T2.1, srm.t$C2T2.3, srm.t$C2T2.4)
C2T3<-cbind(srm.t$C2T3.2, srm.t$C2T3.3, srm.t$C2T3.4)
C2T4<-cbind(srm.t$C2T4.1, srm.t$C2T4.2, srm.t$C2T4.3)
C2T6<-cbind(srm.t$C2T6.1, srm.t$C2T6.3, srm.t$C2T6.4)
TS5<-cbind(srm.t$TS5.1, srm.t$TS5.2, srm.t$TS5.4)
TS6<-cbind(srm.t$TS6.1, srm.t$TS6.2, srm.t$TS6.4)
TS7<-cbind(srm.t$TS7.1, srm.t$TS7.2, srm.t$TS7.3)
TS8<-cbind(srm.t$TS8.1, srm.t$TS8.2, srm.t$TS8.3)
TS9<-cbind(srm.t$TS9.1, srm.t$TS9.2, srm.t$TS9.3, srm.t$TS9.4)

C2T1.avg<-rowMeans(C2T1)
C2T2.avg<-rowMeans(C2T2)
C2T3.avg<-rowMeans(C2T3)
C2T4.avg<-rowMeans(C2T4)
C2T6.avg<-rowMeans(C2T6)
TS5.avg<-rowMeans(TS5)
TS6.avg<-rowMeans(TS6)
TS7.avg<-rowMeans(TS7)
TS8.avg<-rowMeans(TS8)
TS9.avg<-rowMeans(TS9)

C2.avg.df<-data.frame(C2T1.avg, C2T2.avg, C2T3.avg, C2T4.avg, C2T6.avg)
rownames(C2.avg.df)<-rownames(srm.t)
TS.avg.df<-data.frame(TS5.avg, TS6.avg, TS7.avg, TS8.avg, TS9.avg)
rownames(TS.avg.df)<-rownames(srm.t)

C2.avg.df$sum<-rowSums(C2.avg.df)
C2.avg.df2<-subset(C2.avg.df, sum>0, select=-sum)
TS.avg.df$sum<-rowSums(TS.avg.df)
TS.avg.df2<-subset(TS.avg.df, sum>0, select=-sum)

#cluster analysis: incubation
C2.bray<-vegdist(C2.avg.df2, method='bray')
C2.avg.clust<-hclust(C2.bray, method='average')
plot(C2.avg.clust, labels=F)
rect.hclust(C2.avg.clust, h=0.3)
C2.cut<-cutree(C2.avg.clust, h=0.3)
write.csv(C2.cut, 'circadian SRM cluster assignments h 0.3.csv', quote=F)

prot.clust.C2<-read.csv('circadian SRM cluster assignments h 0.3.csv', header=T, row.names=1)
clust.nsaf.C2<-merge(x=prot.clust.C2, y=C2.avg.df2, by='row.names', all.x=T)

melt.C2<-melt(clust.nsaf.C2, id.vars=c('Row.names', 'cluster'))
gg.C2<-
  ggplot(melt.C2, aes(x=variable, y=value, group=Row.names)) +
  geom_line(alpha=0.3) +
  theme_bw() +
  facet_wrap(~cluster, scales='free_y') +
  labs(x='Time Point', y='Averaged Normalized Spectral Abundance Factor') +
  theme(axis.text.x = element_text(angle=90,size=6)) +
  geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6) +
  annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0.3, fill='yellow') +
  annotate("rect", xmin=2.5, xmax=4.5, ymin=0, ymax=Inf, alpha=0.3, fill='blue') +
  annotate("rect", xmin=4.5, xmax=5.5, ymin=0, ymax=Inf, alpha=0.3, fill='yellow') +
  scale_x_discrete('Time point',labels=c('06:00', '12:00', '18:00', '24:00', '36:00'))

#cluster analysis: in situ
TS.bray<-vegdist(TS.avg.df2, method='bray')
TS.avg.clust<-hclust(TS.bray, method='average')
plot(TS.avg.clust, labels=F)
rect.hclust(TS.avg.clust, h=0.3)
TS.cut<-cutree(TS.avg.clust, h=0.3)
write.csv(TS.cut, 'time series SRM cluster assignments h 0.3.csv', quote=F)

prot.clust.TS<-read.csv('time series SRM cluster assignments h 0.3.csv', header=T, row.names=1)
clust.nsaf.TS<-merge(x=prot.clust.TS, y=TS.avg.df2, by='row.names', all.x=T)

melt.TS<-melt(clust.nsaf.TS, id.vars=c('Row.names', 'cluster'))
gg.TS<-
  ggplot(melt.TS, aes(x=variable, y=value, group=Row.names)) +
  geom_line(alpha=0.3) +
  theme_bw() +
  facet_wrap(~cluster, scales='free_y') +
  labs(x='Time Point', y='Averaged Normalized Spectral Abundance Factor') +
  theme(axis.text.x = element_text(angle=90,size=6)) +
  geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6) +
  annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0.3, fill='blue') +
  annotate("rect", xmin=2.5, xmax=4.5, ymin=0, ymax=Inf, alpha=0.3, fill='yellow') +
  annotate("rect", xmin=4.5, xmax=5.5, ymin=0, ymax=Inf, alpha=0.3, fill='blue') +
  scale_x_discrete('Time point',labels=c('22:00', '02:00', '09:00', '15:00', '22:00'))
