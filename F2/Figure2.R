library(ape)
library(ade4)
library(adephylo)
library(treeio)
library(ggtree)
library(phylogram)
library(phytools)
library(data.tree)
library(tidytree)
library(dplyr)
library(paco)
library(ggplot2)
library(seqinr)
library(castor)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load subtrees generated in 'Figure1.R'
pv = read.tree('../shared_data/subtree_full.tree')
plot (pv)

# load host trees created using TimeTree.org
host = read.tree('th2.tree')

# load virus-host association .csv files 

HP = read.csv('HP-list.csv', header=FALSE)
HP
cophy = cophylo(pv, host, assoc=HP, rotate = T)
tangle = plot(cophy)

plot(cophy)


pdf('tanglegram.pdf')
plot(cophy)
dev.off()

HP1 = read.csv('HP-pruned.csv', row.names = 1)
HP2 = read.csv('HP-full.csv', row.names = 1)

host1 = read.tree('th1.tree')
host2 = read.tree('th2.tree')
pv1 = read.tree('../shared_data/subtree_pruned.tree')
pv2 = read.tree('../shared_data/subtree_full.tree')



dist1 = tree_distance(pv1, host1, metric='WassersteinLaplacianSpectrum', normalize=TRUE, NLeigenvalues = 0)
dist2 = tree_distance(pv2, host2, metric='WassersteinLaplacianSpectrum', normalize=TRUE, NLeigenvalues = 0)

dist1
dist2


htree1<-cophenetic(host1)
vtree1<-cophenetic(pv1)

htree2<-cophenetic(host2)
vtree2<-cophenetic(pv2)


# PACo - pruned subtree

D1 <- prepare_paco_data(H=htree1, P=vtree1,HP=HP1)
D1 = add_pcoord(D1,correction = 'lingoes')
D1 = PACo (D1,nperm=1000,seed=12,method='r2', symmetric = FALSE, shuffled = TRUE)
D1 = paco_links(D1)
#D
#p-value
D1$gof$p
D1$gof$ss
S1 = as.data.frame(D1$shuffled)


q5.1 = quantile(S1[,1],c(0.05))
q95.1 = quantile(S1[,1],c(0.95))


p1 = ggplot(S1, aes(x=S1[,1])) +
  geom_density(fill='grey70') +
  geom_vline(xintercept = D1$gof$ss, color= 'red') +
  geom_vline(xintercept = q5.1, color= 'black') +
  geom_vline(xintercept = q95.1, color= 'black') +
  xlab('Procrustes sum of squared residuals')

p1
pdf('subset_paco.pdf') 
p1
dev.off()


# PACo - full subtree

D2 <- prepare_paco_data(H=htree2, P=vtree2,HP=HP2)
D2 = add_pcoord(D2,correction = 'lingoes')
D2 = PACo (D2,nperm=1000,seed=12,method='r2', symmetric = FALSE, shuffled = TRUE)
D2 = paco_links(D2)
#D
#p-value
D2$gof$p
D2$gof$ss
S2 = as.data.frame(D2$shuffled)


q5.2 = quantile(S2[,1],c(0.05))
q95.2 = quantile(S2[,1],c(0.95))


p2 = ggplot(S2, aes(x=S2[,1])) +
  geom_density(fill='grey70') +
  geom_vline(xintercept = D2$gof$ss, color= 'red') +
  geom_vline(xintercept = q5.2, color= 'black') +
  geom_vline(xintercept = q95.2, color= 'black') +
  xlab('Procrustes sum of squared residuals')

p2

pdf('full_paco.pdf') 
p2
dev.off()


