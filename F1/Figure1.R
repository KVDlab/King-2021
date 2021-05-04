# library
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
library(viridis)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load Maximum Likelihood tree 
tree = read.beast('phylotree.phylotree')

# drop clade containing SaPV1, MaegPV1, MaegPV2, LcamPV1  
tree2 = drop.tip(tree, c('SaPV1','MaegPV1','MaegPV2','LcamPV1'))

# add bootstrap node support
fac <- c(as.numeric(tree2@data[['bootstrap']]))
col <- cut(fac, breaks=c(0,50,75,90,99,100))
col <- factor(col, levels=rev(levels(col)))
col <- factor(col, labels=c("100%","91-99%","76-90%","51-75%","0-50%"))
cols = c("#481567FF","#404788FF","#1F968BFF","#95D840FF","#B8DE29FF")

# build full tree with added support values
t = ggtree(tree2, ladderize = TRUE, lwd = 0.5,) +
  geom_nodepoint(aes(color = col), size = 1.5, show.legend = TRUE) + 
  scale_colour_manual(na.translate = F, name="bootstrap support", values=cols) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme_tree2(legend.position = c(0.11, 0.6))

# collapse individual nodes by genus
# Alphapapillomaviridae
t = t %>% collapse(node=441) + 
  geom_point2(aes(subset=(node==441)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==441), label = "Alphapapillomaviridae (83)"),  size = 2, hjust = -.1)

# Betapapillomaviridae
t = t %>% collapse(node=774) + 
  geom_point2(aes(subset=(node==774)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==774), label = "Betapapillomaviridae (55)"),  size = 2, hjust = -.1)

# Gammapapillomaviridae
t = t %>% collapse(node=647) + 
  geom_point2(aes(subset=(node==647)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==647), label = "Gammapapillomaviridae (87)"),  size = 2, hjust = -.1)
t = t %>% collapse(node=729) + 
  geom_point2(aes(subset=(node==729)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==729), label = "Gammapapillomaviridae (13)"),  size = 2, hjust = -.1)

# Mupapillomaviridae
t = t %>% collapse(node=638) + 
  geom_point2(aes(subset=(node==638)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==638), label = "Mupapillomaviridae (3)"),  size = 2, hjust = -.1)

# Upsilonpapillomaviridae
t = t %>% collapse(node=529) + 
  geom_point2(aes(subset=(node==529)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==529), label = "Upsilonpapillomaviridae (7)"),  size = 2, hjust = -.1)

# Omikronpapillomaviridae
t = t %>% collapse(node=542) + 
  geom_point2(aes(subset=(node==542)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==542), label = "Omikronpapillomaviridae (5)"),  size = 2, hjust = -.1)

# Dyopipapillomaviridae
t = t %>% collapse(node=540) + 
  geom_point2(aes(subset=(node==540)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==540), label = "Dyopipapillomaviridae (2)"),  size = 2, hjust = -.1)

# Omegapapillomaviridae
t = t %>% collapse(node=535) + 
  geom_point2(aes(subset=(node==535)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==535), label = "Omegapapillomaviridae (5)"),  size = 2, hjust = -.1)

# Deltapapillomaviridae
t = t %>% collapse(node=558) + 
  geom_point2(aes(subset=(node==558)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==558), label = "Deltapapillomaviridae (15)"),  size = 2, hjust = -.1)
t = t %>% collapse(node=575) + 
  geom_point2(aes(subset=(node==575)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==575), label = "Deltapapillomaviridae (2)"),  size = 2, hjust = -.1)

# Epsilonpapillomaviridae
t = t %>% collapse(node=572) + 
  geom_point2(aes(subset=(node==572)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==572), label = "Epsilonpapillomaviridae (4)"),  size = 2, hjust = -.1)

# Pipapillomaviridae
t = t %>% collapse(node=748) + 
  geom_point2(aes(subset=(node==748)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==748), label = "Pipapillomaviridae (8)"), size = 2, hjust = -.1)

# Taupapillomaviridae
t = t %>% collapse(node=757) + 
  geom_point2(aes(subset=(node==757)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==757), label = "Taupapillomaviridae (13)"), size = 2, hjust = -.1)

# Xipapillomaviridae
t = t %>% collapse(node=829) + 
  geom_point2(aes(subset=(node==829)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==829), label = "Xipapillomaviridae (15)"), size = 2, hjust = -.1)

# Chipapillomaviridae
t = t %>% collapse(node=587) + 
  geom_point2(aes(subset=(node==587)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==587), label = "Chipapillomaviridae (13)"), size = 2, hjust = -.1)

# Dyozetapapillomaviridae
t = t %>% collapse(node=860) + 
  geom_point2(aes(subset=(node==860)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==860), label = "Dyozetapapillomaviridae (2)"), size = 2, hjust = -.1)

# Lambdapapillomaviridae
t = t %>% collapse(node=616) + 
  geom_point2(aes(subset=(node==616)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==616), label = "Lambdapapillomaviridae (12)"), size = 2, hjust = -.1)

# Etapapillomaviridae
t = t %>% collapse(node=859) + 
  geom_point2(aes(subset=(node==859)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==859), label = "Etapapillomaviridae (2)"), size = 2, hjust = -.1)

# Treisepsilonpapillomaviridae
t = t %>% collapse(node=857) + 
  geom_point2(aes(subset=(node==857)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==857), label = "Treisepsilonpapillomaviridae (2)"), size = 2, hjust = -.1)

# Zetapapillomaviridae
t = t %>% collapse(node=581) + 
  geom_point2(aes(subset=(node==581)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==581), label = "Zetapapillomaviridae (2)"), size = 2, hjust = -.1)

# Dyorhopapillomaviridae
t = t %>% collapse(node=579) + 
  geom_point2(aes(subset=(node==579)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==579), label = "Dyorhopapillomaviridae (3)"), size = 2, hjust = -.1)

# Dyoiotapapillomaviridae
t = t %>% collapse(node=578) + 
  geom_point2(aes(subset=(node==578)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==578), label = "Dyoiotapapillomaviridae (3)"), size = 2, hjust = -.1) +
  geom_tiplab(aes(subset=(node==133), label = "Dyoiotapapillomaviridae (3)"), size = 2, hjust = -.1)

# Iotapapillomaviridae
t = t %>% collapse(node=583) + 
  geom_point2(aes(subset=(node==583)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==583), label = "Iotapapillomaviridae (4)"), size = 2, hjust = -.1)

# Dyonupapillomaviridae
t = t %>% collapse(node=599) + 
  geom_point2(aes(subset=(node==599)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==599), label = "Dyonupapillomaviridae (3)"), size = 2, hjust = -.1)

# Dyokappapapillomaviridae
t = t %>% collapse(node=604) + 
  geom_point2(aes(subset=(node==604)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==604), label = "Dyokappapapillomaviridae (4)"), size = 2, hjust = -.1) +
  geom_tiplab(aes(subset=(node==168), label = "Dyokappapapillomaviridae (1)"), size = 2, hjust = -.1)

# Rhopapillomaviridae
t = t %>% collapse(node=609) + 
  geom_point2(aes(subset=(node==609)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==609), label = "Rhopapillomaviridae (4)"), size = 2, hjust = -.1)

# Dyoomikronpapillomaviridae
t = t %>% collapse(node=523) + 
  geom_point2(aes(subset=(node==523)), shape=23, size=1.5) +
  geom_tiplab(aes(subset=(node==523), label = "Dyoomikronpapillomaviridae (4)"), size = 2, hjust = -.1)

t2 = t + 
  theme_tree2(legend.position = c(0.1, 0.8)) +
  geom_tiplab(size = 2, align = TRUE) +
  geom_treescale(x = 2.3, y = 0, fontsize = 3, width = 0.2) +
  xlim(0,4)
t2

pdf('MLtreeR.pdf')
t2
dev.off()

# extract subtrees for tanglegram analysis (Figure2.R)

phylo = as.phylo(t2)

plot(phylo)

pv.pruned = extract.clade(phylo, 614)
plot(pv.pruned)

tip=c('PlPV1','LwPV1','ElPV1','AmPV4','ElPV1','AmPV4','LwPV1','LwiePV1','CcrPV1','PcPV1','UuPV1','PlpPV1','FcaPV1','CPV1','CPV6')
pv.pruned= drop.tip(pv.pruned, tip)
plot(pv.pruned)



pv.full = extract.clade(phylo, 612)
plot(pv.full)

write.tree(pv.pruned, '../shared_data/subtree_pruned.tree')
write.tree(pv.full, '../shared_data/subtree_full.tree')



# build heatmap for clade containing Yang and Yin groups 

library(pheatmap)
library(seqinr)

# distance matrix generated from pairwise L1 sequence alignments
L1 = read.csv('pairwise_matrix.csv', row.names = 1)
mat = as.matrix(L1)

# create heatmap
h = pheatmap(mat, color = viridis(10, begin = 0.15, end = 0.88, direction = 1), cellheight = 12, cellwidth = 12, border_color = 'black',
              fontsize_row = 7, fontsize_col = 7)

pdf('L1heatmap.pdf') 
h 
dev.off()



