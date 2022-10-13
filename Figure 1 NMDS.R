#load biostats
source('biostats.R')
library(dplyr)
library(vegan)

#read in data
qspec=read.csv("ABACUS_copepods_metazoan_output.csv")

# keep only proteins that are not contaminants or symbionts
copepod.prot<-subset(qspec, grepl(paste('TRINITY', collapse="|"), qspec$PROTID))

# keeping rows where ALL_NUMPEPSUNIQ > 1
subsetted=subset(copepod.prot, ALL_NUMPEPSUNIQ > 1)

# keeping columns that are PROTID, PROTLEN, and the ADJNSAF columns for each sample
filtered <- select(subsetted, contains(c("PROTID", "PROTLEN", 'ADJNSAF'))) 
filtered1 <- filtered[c(1:2, 64:78)]

cop.dat <- filtered1[3:17] 
rownames(cop.dat)<-filtered1$PROTID

#transpose data
cop.t<-as.data.frame(t(cop.dat))

#log-transform data
cop.tra<-data.trans((cop.t+1), method='log', plot=F)

#run NMDS on a bray-curtis dissimilarity matrix
nmds.cop<-metaMDS(cop.tra, distance='bray', trymax=10, autotransform=F)

#make object grps that is a list of group assignments for each of the samples
#these need to follow the order of your columns in the original data frame
grps<-c("06:00", "12:00 D3", "12:00 D3", "18:00","12:00 D3", "18:00","06:00", "06:00","18:00", "00:00", "12:00 D4", "00:00", "12:00 D4", "00:00", "12:00 D4")

#plot NMDS
fig.cop<-ordiplot(nmds.cop, choices=c(1,2), type='none', display='sites')
points(fig.cop, 'sites', pch=c(rep(21,9), rep(22,6)) , bg=c("gold", "gold2", "gold2", "blue", "gold2", "blue", "gold", "gold", "blue", "blue4", "gold2", "blue4", "gold2", "blue4", "gold2") , col='grey50', cex=1.5)
ordihull(fig.cop, groups=grps, draw='lines',col='grey75', label=T)

#ANOSIM
#this tests statistics of group assignments

#row standardize
cop.row<-data.stand(cop.t, method='total', margin='row', plot=F)

#bray-curtis dissimilarity
cop.d<-vegdist(cop.row, 'bray')

#run anosim
#save the output of summary
all.anosim<-anosim(cop.d, grouping=grps)
summary(all.anosim)

#export the eigenvalues of proteins
#this is essentially the weight, or value, of each protein on each axis
#if a protein has a strong weight, it is a strong contributor to the group separation on that axis
eigen<-envfit(nmds.cop$points, cop.tra, perm=1000)
#this may take a long time
eigen