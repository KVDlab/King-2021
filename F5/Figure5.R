# library
library(ggtree)
library(car)
library(FSA)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(seqinr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#analysis groups
Yang = c("TbraPV3","EsPV1","TbraPV1","TbraPV2","EsPV3","MscPV2")
Yin = c("RfPV1","EhPV1","EdPV1","HPV41")
other  = c("TbelPV1","HPV204","HPV1","HPV63","CcanPV1","OcPV1","SfPV1","TePV1","CPV1","LwiePV1","CcrPV1","FcaPV1","PlpPV1","UuPV1","LrPV1","PcPV1","CPV6","AmPV4","LwPV1","ElPV1","PlPV1")
all = c(Yin,Yang,other)

# RSCU (Relative Synonymous Codon Usage)
seqs = read.fasta('concatenated.fas') 
rscu = NULL;
for (i in seqs) {
  tmp = uco(i, frame = 0, index = c('rscu'), as.data.frame = FALSE, NA.rscu = NA)
  rscu = rbind(rscu, tmp)
}
rscu = data.frame(rscu)
row.names(rscu) = names(seqs)
rscu$virus = row.names(rscu)
write.csv(rscu, 'rscu.csv') 


# convert rscu.cg data to long format and select codons for amino acids in which one or more contains a 'CG'
rscu = read.csv('rscu.csv', row.names = 1) # read RSCU data .csv file
AA = c('gct','gcc','gca','gcg','cct','ccc','cca','ccg','cgt','cgc','cga','cgg','aga','agg','tct','tcc','tca','tcg','agt','agc','act','acc','aca','acg')
rscu.cg = rscu[AA] # create subset of amino acids encoded by CG-containing codons + synonymous codons
rscu.cg$virus = row.names(rscu.cg)
rscu.cg$virus <- factor(rscu.cg$virus)
long = gather(rscu.cg, codon, rscu, gct:acg, factor_key=FALSE) # convert RSCU table to long format


# create data table groups for plotting 
letters = c(rep('A', 124), rep('P', 124), rep('R', 186), rep('S',186), rep('T', 124))
long$aa = letters # amino acid letter codes
group = c(rep('other', 21),'yang','yin','yang','yang','yin','yang','yang','yang','yin','yin')
group = rep(group, 24) # analysis groups (yang, yin, other)
long$group = group
long[[2]] <- toupper(long[[2]]) # convert codons to upper case

write.csv(long, 'rscu_plots.csv')


# create boxplot of synonymous-CG codon RSCU values faceted by amino acid
rscu.plots = read.csv('rscu_plots.csv', row.names = 1) # read RSCU .csv file for plotting
head(rscu.plots) # quick check of data

rscu.box = ggplot(rscu.plots, aes(x=codon, y=rscu, fill=group)) + 
  geom_boxplot(outlier.shape = 21) +
  geom_hline(aes(yintercept = 0.6, colour = 'darkred')) +
  facet_wrap(~aa, scale="free_x") +
  theme_bw() 

rscu.box

pdf('rscu.pdf')
rscu.box
dev.off() # write figure to .pdf (further editing done in Adobe Illustrator)


# stats


d1 = read.csv("rscu.csv")  # read csv file

#add group
d1$color <- ifelse(d1$X %in% Yang, 'Yang', 
                   ifelse (d1$X %in% Yin,'Yin',
                           ifelse (d1$X %in% other, 'other', "NA")))
d1$color <- factor(d1$color,levels = c("Yang", "Yin", "other","NA"))
d1 =(d1 %>%
       filter(color != "NA") %>%
       select(everything()))
d1 <- d1[c("gcg","ccg","cga","cgc","cgg","cgt", "tcg","acg","color")]
d1$color <- factor(d1$color)

formulae <- lapply(colnames(d1)[1:8], function(x) as.formula(paste0(x, " ~ color")))

res <- lapply(formulae, function(x) summary(aov(x, data = d1)))
tuk <- lapply(formulae, function(x) TukeyHSD(aov(x, data = d1)))
names(res) <- format(formulae)
names(tuk) <- format(formulae)

unlist(lapply(res, function(x) x[[1]]$"Pr(>F)"[1]))


tuk

dat = read.csv("input4R.csv", header = F)
d = transpose(dat[,-1])
colnames(d) <- dat[,1]
d
#convert to long format
d1 <- reshape(d, 
              direction = "long",
              varying = list(names(d)[1:6]),
              v.names = c("cusp"),
              timevar = "groups",
              times = c("Yang-Yang","Yin-Yin","other-other","Yang-Yin","Yang-other","other-Yin"))

d1$groups <- factor(d1$groups, levels = c("Yang-Yang","Yin-Yin","other-other","Yang-Yin","Yang-other","other-Yin"))



#build boxplot
box <- ggplot(d1, aes(groups,cusp)) +
  geom_jitter(height = 0, width = 0.05, size = 0.5, color='grey') +
  geom_boxplot() +
  theme_bw()

box

pdf(file="box.pdf")
box
dev.off()

#test anova assumptions
d1.lm <- lm(cusp ~ groups, data = d1)
d1.av <- aov(d1.lm)
plot(d1.av, 1)
leveneTest(cusp ~ groups, data = d1)
plot(d1.av, 2)
d1.av_residuals <- residuals(object = d1.av)
shapiro.test(x = d1.av_residuals)

#anova rejected, so use Kruskal Wallis and Dunn's post-test
kruskal.test(cusp ~ groups, data = d1) #p-value = 0.0001164
PT = dunnTest(cusp ~ groups,
              data=d1,
              method="bh") #Benjamini-Hochberg method
PT

#             Comparison           Z      P.unadj        P.adj
#1    other-other - other-Yin_  1.12611170 2.601182e-01 4.877217e-01
#2    other-other - Yang-other  0.66610882 5.053415e-01 6.891021e-01
#3    other-Yin_ - Yang-other -0.33405495 7.383381e-01 8.519286e-01
#4    other-other - Yang-Yang  4.83113859 1.357545e-06 2.036317e-05
#5    other-Yin_ - Yang-Yang  3.75591617 1.727085e-04 8.635425e-04
#6    Yang-other - Yang-Yang  3.94910641 7.844348e-05 5.883261e-04
#7    other-other - Yang-Yin_  1.66240023 9.643255e-02 2.892977e-01
#8    other-Yin_ - Yang-Yin_  1.21738665 2.234571e-01 4.788367e-01
#9    Yang-other - Yang-Yin_  1.35993285 1.738512e-01 4.346279e-01
#10   Yang-Yang - Yang-Yin_ -1.04363270 2.966554e-01 4.944256e-01
#11   other-other - Yin_-Yin_  0.46723624 6.403309e-01 8.004136e-01
#12   other-Yin_ - Yin_-Yin_  0.03069903 9.755096e-01 9.755096e-01
#13   Yang-other - Yin_-Yin_  0.19075227 8.487197e-01 9.093425e-01
#14   Yang-Yang - Yin_-Yin_ -2.27250026 2.305631e-02 8.646117e-02
#15   Yang-Yin_ - Yin_-Yin_ -0.92281513 3.561036e-01 5.341553e-01


aa = read.csv("aa_composition_values.csv")
aa

aa$color <- factor(aa$color, levels = c("Yang","Yin","Other"))

means.barplot = ggplot(aa, aes(x=as.factor(AA), y=mean, fill=color, order=as.factor(color))) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))

means.barplot

pdf(file="means.barplot.pdf")
means.barplot
dev.off()



