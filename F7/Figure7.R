# library
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(ggtree)
library(phytools)
library(plyr)
library(ggpubr)
library(forcats)
library(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# GENERATE .CSVs FOR PLOTS

# fourmer OE ratio barplot for Yang group
OE = read.csv('../shared_data/4_PaVE.complete-genomes.fas_ratio.csv', row.names = 1)
# filter tetramers N-CG-N
tOE = t(OE)
tOE = data.frame(tOE)
tOE['motif'] = row.names(tOE)
tOEf = dplyr::filter(tOE, grepl('.CG.', motif))
tOEf = tOEf[-436]
OEf = t(tOEf)
OEf = data.frame(OEf)

# filter Yang bat PVs
Yang = c("TbraPV3","EsPV1","TbraPV1","TbraPV2","EsPV3","MscPV2")
Yin = c("RfPV1","EhPV1","EdPV1","HPV41")
other  = c("TbelPV1","HPV204","HPV1","HPV63","CcanPV1","OcPV1","SfPV1","TePV1","CPV1","LwiePV1","CcrPV1","FcaPV1","PlpPV1","UuPV1","LrPV1","PcPV1","CPV6","AmPV4","LwPV1","ElPV1","PlPV1")
all = c(Yin,Yang,other)

OEf['pv'] = row.names(OEf)
OEf.yang = subset(OEf, pv %in% Yang)
OEf.yang = OEf.yang[-17]

# get averages for each Yang tetramer OE value and create table for figure 5A
OEf.yang['YANG',] = colMeans(x=OEf.yang, na.rm = TRUE)
yang.avg = data.frame(OEf.yang['YANG',])
yang.avg = data.frame(yang.avg)
long1 = gather(yang.avg, motif, OE.ratio, ACGA:TCGT, factor_key = TRUE) # convert OE averages table to long format

write.csv(long1, 'yang_4OE.csv') # .csv containing the averaged Yang PV OE (observed/expected) values for each N-CG-N tetramer



# get averages for each Yin, Yang, and other tetramer OE values and calculate YANG/YIN, YANG/OTHER, and YIN/OTHER OE ratios 
OEf.yin = subset(OEf, pv %in% Yin)
OEf.yin = OEf.yin[-17]
OEf.yin['YIN',] = colMeans(x=OEf.yin, na.rm = TRUE)
yin.avg = data.frame(t(OEf.yin['YIN',]))

OEf.other = subset(OEf, pv %in% other)
OEf.other = OEf.other[-17]
OEf.other['OTHER',] = colMeans(x=OEf.other, na.rm = TRUE)
other.avg = data.frame(t(OEf.other['OTHER',]))

OE.average = data.frame(t(yang.avg), yin.avg, other.avg)
OE.average['YANG.YIN'] = OE.average$YANG/OE.average$YIN
OE.average['YANG.OTHER'] = OE.average$YANG/OE.average$OTHER
OE.average['YIN.OTHER'] = OE.average$YIN/OE.average$OTHER



# create table for figure 5B
comp.plots = data.frame(OE.average[4], OE.average[5], OE.average[6])
comp.plots$motif = row.names(comp.plots)
comp.plots$motif = factor(comp.plots$motif)
long2 = gather(comp.plots, group, OE.ratio, YANG.YIN:YIN.OTHER, factor_key = TRUE) # convert OE averages table to long format

write.csv(long2, '4mer-comp_plots.csv') # .csv containing the averaged Yang, Yin, and other OE (observed/expected) ratios for each N-CG-N tetramer



# PLOTS

# FIGURE 5A: calculate standard deviation (sd) and standard error (se) 
tyang.avg = data.frame(t(yang.avg))
sd = sd(tyang.avg$YANG) # standard deviation
se = sd(tyang.avg$YANG)/sqrt(length(tyang.avg$YANG)) # standard error

# create barplot of Yang bat N-CG-N tetramer OE ratios
yang.plot = read.csv('yang_4OE.csv', row.names = 1)
y1 = ggplot(data = yang.plot, aes(x = fct_reorder(motif, OE.ratio), y = OE.ratio, fill = motif)) +
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.7), show.legend = FALSE) +
  geom_errorbar(aes(x=motif, ymin=OE.ratio-se, ymax=OE.ratio+se), width=0.3, colour="black", alpha=0.9, size=1) +
  theme_bw() 
y1

pdf('barplot_F5A.pdf')
y1
dev.off() # write figure to .pdf 




# FIGURE 5B: create barplot of comparative OE ratios for each N-CG-N tetramer
comp.plots = read.csv('4mer-comp_plots.csv', row.names = 1) # read .csv containing ratios (yang/other, yang/yin, yin/other) of the OE (observed/expected) ratio for each N-CG-N tetramer
y2 = ggplot(data = comp.plots, aes(x = motif, y = OE.ratio, fill = group)) +
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.7), show.legend = TRUE) +
  theme_bw() 
y2

pdf('barplot_F5B.pdf')
y2
dev.off() # write figure to .pdf 




#data generated from python
d1 = read.csv("ratio_test-4.csv")  # read csv file
#convert to long format
d2 <- reshape(d1, 
        direction = "long",
        varying = list(names(d1)),
        v.names = c("proportion"),
        timevar = "dinucleotide",
        times = colnames(d1))
#select .CG.
d3 = (d2 %>%
        filter(grepl(".CG.", dinucleotide)) %>%
        select(dinucleotide,proportion))


#data generated from python
d4 = read.csv("ratio_sample-4.csv")  # read csv file
#convert to long format
d5 <- reshape(d4, 
              direction = "long",
              varying = list(names(d4)),
              v.names = c("proportion"),
              timevar = "dinucleotide",
              times = colnames(d4))
#select .CG.
d6 = (d5 %>%
        filter(grepl(".CG.", dinucleotide)) %>%
        select(dinucleotide,proportion))



#build violin plot of dinucs
y3 <- ggplot(d3, aes(x = fct_reorder(dinucleotide, proportion),proportion)) +
  geom_violin() +
  geom_point(data = d6, aes(x = fct_reorder(dinucleotide, proportion),proportion), color='red') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylim(0,1.5)

y3

pdf("violin.pdf")
y3
dev.off() # write figure to .pdf 

ACGT =  (d6 %>%
        filter(grepl("ACGT", dinucleotide)) %>%
        select(dinucleotide,proportion))
ACGT

#dinucleotide proportion
#         ACGT  0.2771289

ACGT_dis = d1$ACGT
min(ACGT_dis)
#0.235841
quantile(ACGT_dis, probs = c(0.01, 0.05)) 
# 1%        5% 
# 0.2921766 0.3378846 


y4 = grid.arrange(y1, y2, y3, nrow = 3)

pdf('figure5.pdf')
y4
dev.off() # write figure panel to .pdf (further editing done in Adobe Illustrator)

