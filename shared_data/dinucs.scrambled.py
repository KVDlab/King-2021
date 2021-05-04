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
running_count = []
for r in range (1000):
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
                    #print (temp)

        ratios[v] = (temp[:-1])

    filterByKey = lambda keys: {x: ratios[x] for x in keys}
    yang = filterByKey(yang)

    yang_mean = (mean_of_lists(list(yang.values())))
    
    
    yin = filterByKey(yin)
    yin_mean = (mean_of_lists(list(yin.values())))

    for n, i in enumerate(yin_mean):
        if i == 0:
            yin_mean[n] = 0.0000000000000001
    
    
    
    results = [i / j for i, j in zip(yang_mean, yin_mean)]
    running_count.append (results)

with open("ratio_test-"+str(word_size)+".csv","w") as out:
    dinuc=[]
    for i in itertools.product("ACGT", repeat=word_size):
        dinuc.append("".join(i))
    print (",".join(dinuc), file=out)
    for r in (running_count):
        print (",".join(map(str,r)),file=out)


with open ("ratio_test-"+str(word_size)+"-percentiles.csv", "w") as out:
    print (",".join(map(str,["hexamer","n","mean","1%","5%","min","max"])),file=out)
    a = zip(*csv.reader(open("ratio_test-"+str(word_size)+".csv", "r")))
    for l in a:
        m = re.search(motif,l[0])
        if m:
            l_float = [float(i) for i in l[1:]]
            l_float = [i for i in l_float if i < 5]
            l_float = [i for i in l_float if i > 0]

            #print (l[0], l_float)
            #print (",".join(map(str,[l[0], len(l_float), np.mean(l_float), np.percentile(l_float, 1), np.percentile(l_float, 5), min(l_float)])))
            print (",".join(map(str,[l[0], len(l_float), np.mean(l_float), np.percentile(l_float, 1), np.percentile(l_float, 5), min(l_float),max(l_float)])),file=out)







