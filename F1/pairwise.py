from Bio import SeqIO
from Bio.Seq import Seq
import os, itertools, csv, random, numpy


with open("pairwise.csv","w") as pw:
    viruses=[]
    for r in SeqIO.parse("L1nt.fas","fasta"):
        viruses.append(r.id)
    
    #viruses=['HPV41','EdPV1']
    
    for c in (list(itertools.combinations(viruses, 2))):
        with open("temp.fas","w") as out:
            for r in SeqIO.parse("L1nt.fas","fasta"):
                if r.id in c:
                    print (f'>{r.id}\n{r.seq.translate()}',file=out)
    
        os.system("linsi temp.fas > temp.linsi.fas")
        pairwise = []
        for aa in SeqIO.parse("temp.linsi.fas","fasta"):
            x=0
            temp=[]
            for n in SeqIO.parse("L1nt.fas","fasta"):
                if n.id == aa.id:
                    for a in aa.seq:
                        if a != "-":
                            y=x*3
                            #print (n.id, n.seq[y:y+3], n.seq[y:y+3].translate(), a)
                            temp.append(str(n.seq[y:y+3]))
                            x=x+1
                        else:
                            #print (n.id, "fail"+a)
                            temp.append("---")
            pairwise.append ("".join(temp))
        count=0
        for p in range(len(pairwise[0])):
            if (pairwise[0][p] == pairwise[1][p]):
                count=count+1
        print (",".join([c[0], c[1], str(100*(count/len(pairwise[0])))]),file=pw)
            
names=[]
similarities={}
for r in csv.reader(open('pairwise.csv',"r")):
    if r[0] not in names:
       names.append(r[0])
    if r[1] not in names:
        names.append(r[1])
    name = r[0]+"|"+r[1]
    similarities[name] = r[2]




matrix=[]
for x in itertools.product(names,repeat=2): 
    if x[0] == x[1]:
        matrix.append (100)
    else:
        name=x[0]+"|"+x[1]
        if name in similarities:
            matrix.append (similarities[name])
        else:
            name=x[1]+"|"+x[0]
            matrix.append (similarities[name])


with open("pairwise_matrix.csv","w") as out:
    print (",".join(['']+names), file=out)
    for n in range(0,len(matrix),10):
        print (",".join(map(str,[names[int(n/10)]]+matrix[n:n+10])), file=out)



