import itertools, os, csv
from Bio import SeqIO
from Bio.SeqUtils import GC


in_files = ["PaVE.complete-genomes.fas"]
for in_file in in_files:
    wordsize=[2,4]
    for w in wordsize:
        dinuc=[]
        for i in itertools.product("ACGT", repeat=w):
            dinuc.append("".join(i))
        topline=[""]+dinuc+["%GC"]
    
        print (dinuc)
        print ("calculating the frequency of sequences with wordsize "+str(w))
        with open(str(w)+"_"+in_file+"_count.csv","w") as c, open(str(w)+"_"+in_file+"_observed.csv","w") as o, open(str(w)+"_"+in_file+"_expected.csv","w") as e, open(str(w)+"_"+in_file+"_ratio.csv","w") as r:
            print (",".join(topline),file=c)
            print (",".join(topline),file=o)
            print (",".join(topline),file=e)
            print (",".join(topline+["%GC"]),file=r)
            for record in SeqIO.parse(in_file,"fasta"):
                with open("temp.fasta.fas","w") as temp:
                    ID=record.id
                    print (ID)
                    print (f">{ID}\n{record.seq}",file=temp)
                    #print (ID, GC(record.seq))
                    
                os.system( "compseq -sequence temp.fasta.fas -word "+str(w)+" -outfile temp.comp.txt -calcfreq Y" )
            
                with open("temp.comp.txt","rU") as f:
                    c1=[]
                    o1=[]
                    e1=[]
                    r1=[]
                    
                    c1.append(ID)
                    o1.append(ID)
                    e1.append(ID)
                    r1.append(ID)
                    
                    for line in f:
                        line=line.strip().split("\t")
                        if line[0] in dinuc:
                            c1.append(line[1])
                            o1.append(line[3])
                            e1.append(line[4])
                            r1.append(line[5])
                    r1.append(str(GC(record.seq)))
                    #print (results)
                    print (",".join(c1),file=c)
                    print (",".join(o1),file=o)
                    print (",".join(e1),file=e)
                    print (",".join(r1),file=r)

