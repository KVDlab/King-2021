import csv, re, os, itertools
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

viruses = []

for r in SeqIO.parse("concatenated.fas","fasta"):
    if r.id not in viruses:
        viruses.append(r.id)

red_bats = ["TbraPV2","EsPV1","TbraPV1","TbraPV3","EsPV3","MscPV2"]
related = ["RfPV1","EhPV1","EdPV1","HPV41"]
other = []





for v in viruses:
    if v not in red_bats and v not in related:
        other.append(v)



red_2 = (list(itertools.permutations((red_bats),2)))
related_2 = (list(itertools.permutations((related),2)))
other_2 = (list(itertools.permutations((other),2)))
red_related = list(itertools.product(red_bats, related))
red_other = list(itertools.product(red_bats, other))
related_other = list(itertools.product(related, other))








for r in SeqIO.parse("concatenated.fas","fasta"):
        with open("temp.fas","w") as fas:
            print (f">{r.id}\n{r.seq}", file = fas)
        os.system("cusp temp.fas "+r.id+".cusp")
        

cusp = []

for file in os.listdir("."):
    if file.endswith(".cusp"):
        cusp.append(file)
for l in (list(itertools.combinations(cusp, 2))):
   os.system ("codcmp "+l[0]+" "+l[1]+" "+l[0].split(".")[0]+"-"+l[1].split(".")[0]+".codcmp")
 
codcmp=[]
 
for file in os.listdir("."):
    if file.endswith(".codcmp"):
        codcmp.append(file)
        
        
with open("codcmp.out.csv","w") as out:
    print (",".join(["virus 1","virus 2","Sum Squared Difference","Mean Squared Difference","Root Mean Squared Difference","Sum Difference ","Mean Difference ","Codons not appearing"]),file=out)
    for c in codcmp:
        with open(c,"r") as f:
            v1 = (c.split(".")[0].split("-")[0])
            v2 = (c.split(".")[0].split("-")[1])
            temp = [v1,v2]
            for l in f:
                l = l.strip()
                if "=" in l:
                   temp.append (l.split (" = ")[1])
            print (",".join(temp),file=out)
        
os.system ("rm *.cusp")
os.system ("rm *.codcmp")


RMSD = []
viruses = []
z=30
x=29
y=0
for r in csv.reader(open("codcmp.out.csv","r")):
    if "virus" not in r[0]:
        RMSD.append(r[4])
        if r[0] not in viruses:
            viruses.append(r[0])
        if r[1] not in viruses:
            viruses.append(r[1])



with open("codcmp.matrix.csv","w") as out:
    top = [""]
    top = top+viruses
    print (",".join(top), file = out)
    for r in range(x):
        t = z-x
        temp = [viruses[r]]
        temp =temp+ (t * [""])
        temp = temp + (RMSD[y:y+x])
        print (temp)
        print (",".join(temp), file = out)
        y=y+x
        x=x-1

red_red = []
related_related = []
other_other = []
related_red = []
other_red = []
other_related = []


for r in csv.reader(open("codcmp.out.csv","r")):
    x = ((r[0],r[1]))
    if x in red_2:
        red_red.append( (r[4]) )
    elif x in related_2:
        related_related.append( (r[4]) )
    elif x in other_2:
        other_other.append( (r[4]) )
    elif x in red_related:
        related_red.append((r[4]))
    elif x in red_other:
        other_red.append((r[4]))
    elif x in related_other:
        other_related.append((r[4]))

with open("input4R.csv","w") as out:               
    print (f",".join(["red-red"]+red_red), file = out)
    print (f",".join(["related-related"]+related_related), file = out)
    print (f",".join(["other-other"]+other_other), file = out)
    print (f",".join(["red-related"]+related_red), file = out)
    print (f",".join(["red-other"]+other_red), file = out)
    print (f",".join(["other-related"]+other_related), file = out)




#calculate aa composition as described https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2590925/

aa=["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T"]




with open("aa_composition.csv","w") as out:
    print (f",".join(["virus"]+aa), file =out)
    for r in SeqIO.parse("concatenated.fas","fasta"):
        protein = str(r.seq.translate())
        aa_count=[]
        for a in aa:
            aa_count.append( 100*(protein.count(a)/len(protein)) )
        
        print (f",".join(map(str,[r.id]+aa_count )), file =out)
        
        
        
P={}
S={}
T={}
A={}
R={}
        
for r in csv.reader(open("aa_composition.csv")):
    if r[0] != 'virus':
       P[r[0]] = r[2]
       A[r[0]] = r[3]
       R[r[0]] = r[14]
       S[r[0]] = r[19]
       T[r[0]] = r[20]
 
P_RB = []
P_Rel = []
P_O = [] 
       
for p in P:
    if p in red_bats:
        P_RB.append ( float((P[p])) )
    elif p in related:
        P_Rel.append ( float(P[p]))
    else:
        P_O.append (float(P[p]))
        
A_RB = []
A_Rel = []
A_O = [] 
       
for a in A:
    if a in red_bats:
        A_RB.append ( float(A[a]) )
    elif a in related:
        A_Rel.append ( float(A[a]) )
    else:
        A_O.append (float(A[a]))


R_RB = []
R_Rel = []
R_O = [] 
       
for r in R:
    if r in red_bats:
        R_RB.append (float(R[r])) 
    elif r in related:
        R_Rel.append (float(R[r]))
    else:
        R_O.append (float(R[r]))

T_RB = []
T_Rel = []
T_O = [] 
       
for t in T:
    if t in red_bats:
        T_RB.append (float(T[t]))
    elif t in related:
        T_Rel.append (float(T[t]))
    else:
        T_O.append (float(T[t]))
        
S_RB = []
S_Rel = []
S_O = [] 
       
for s in S:
    if s in red_bats:
        S_RB.append (float(S[s]))
    elif s in related:
        S_Rel.append (float(S[s]))
    else:
        S_O.append (float(S[s]))


        
with open("aa_composition_values.csv","w") as out:
    print (f",".join(map(str,["AA","color","mean","sd"]+[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])), file=out)
    
    print (f",".join(map(str,["P","Yang",np.mean(P_RB), np.std(P_RB)]+P_RB)), file=out)
    print (f",".join(map(str,["P","Yin",np.mean(P_Rel), np.std(P_Rel)]+P_Rel)), file=out)
    print (f",".join(map(str,["P","Other",np.mean(P_O), np.std(P_O)]+P_O)), file=out)
    
    print (f",".join(map(str,["S","Yang",np.mean(S_RB), np.std(S_RB)]+S_RB)), file=out)
    print (f",".join(map(str,["S","Yin",np.mean(S_Rel), np.std(S_Rel)]+S_Rel)), file=out)
    print (f",".join(map(str,["S","Other",np.mean(S_O), np.std(S_O)]+S_O)), file=out)
    
    print (f",".join(map(str,["T","Yang",np.mean(T_RB), np.std(T_RB)]+T_RB)), file=out)
    print (f",".join(map(str,["T","Yin",np.mean(T_Rel), np.std(T_Rel)]+T_Rel)), file=out)
    print (f",".join(map(str,["T","Other",np.mean(T_O), np.std(T_O)]+T_O)), file=out)

    print (f",".join(map(str,["R","Yang",np.mean(R_RB), np.std(R_RB)]+R_RB)), file=out)
    print (f",".join(map(str,["R","Yon",np.mean(R_Rel), np.std(R_Rel)]+R_Rel)), file=out)
    print (f",".join(map(str,["R","Other",np.mean(R_O), np.std(R_O)]+R_O)), file=out)
    
    print (f",".join(map(str,["A","Yang",np.mean(A_RB), np.std(A_RB)]+A_RB)), file=out)
    print (f",".join(map(str,["A","Yon",np.mean(A_Rel), np.std(A_Rel)]+A_Rel)), file=out)
    print (f",".join(map(str,["A","Other",np.mean(A_O), np.std(A_O)]+A_O)), file=out)

    


    
