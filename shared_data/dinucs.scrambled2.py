import itertools, os, csv, random, re
from altschulEriksonDinuclShuffle import dinuclShuffle
from Bio import SeqIO
from Bio.SeqUtils import GC

from itertools import zip_longest
from statistics import mean

import numpy as np


def mean_of_lists(lst):
    return [mean(y for y in t if y != -1) for t in zip_longest(*lst, fillvalue=-1)]



viruses= ["HPV41", "EdPV1", "RfPV1", "EhPV1","TbraPV1","TbraPV2","TbraPV3","EsPV3","MscPV2","EsPV1"]
yang = ["TbraPV1","TbraPV2","TbraPV3","EsPV3","MscPV2","EsPV1"]
yin = ["HPV41", "EdPV1", "RfPV1", "EhPV1"]

word_size = 4
motif = ".CG."

in_files = ["PaVE.complete-genomes.fas"]

with open("new_stats.csv","w") as out:
    dinuc=[]
    for i in itertools.product("ACGT", repeat=word_size):
        dinuc.append("".join(i))
    print (",".join(dinuc), file=out)
    
    for r in range (100):
        print (r)
        ratios = {}
        for v in viruses:
            for in_file in in_files:
                for record in SeqIO.parse(in_file,"fasta"):
                    ID=record.id
                    if ID == v:
                        print (v)
                        with open("temp.fasta.fas","w") as temp:
                                seq =  str(record.seq).upper()
                                new_seq =  (dinuclShuffle(seq))
                                print (f">{ID}\n{''.join(new_seq)}", file = temp)
                        os.system( "compseq -sequence temp.fasta.fas -word "+str(word_size)+" -outfile temp.comp.txt -calcfreq Y" )
                        temp = []
                        with open("temp.comp.txt","r") as f:
                            for line in f:
                                line = line.strip().split("\t")
                                if len(line) > 4:
                                    if "#" not in line[0]:
                                        temp.append (float(line[-1]))
                        print (",".join(map(str,[v]+temp[:-1])), file = out)

        #ratios[v] = (temp[:-1])
    #print (ratios)
