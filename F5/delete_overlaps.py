import csv, re

def trim_E6 (genes):
    A = genes['E6'].split("..")
    Ar = range(int(A[0]),int(A[1]))
    if 'E7' in genes:
        B = genes['E7'].split("..")
        Br = range(int(B[0]),int(B[1]))
    else:
        Br = [0]
    C=None
    if 'E10' in genes:
        C=genes['E10'].split("..")
    CG = genes['CG']

    
   
    
    overlap = (set(Ar).intersection(set(Br)))
    if len(overlap) % 3 == 0:
        delete = len(overlap)+3
        if C == None:
            return (CG[int(A[0])-1:int(A[1])-delete])
        else:
            return (CG[int(A[0])-1:int(C[0])-2]+CG[int(C[1])+2:int(A[1])-delete])
    elif len(overlap) % 3 == 1:
        delete = len(overlap)+2
        if C == None:
            return (CG[int(A[0])-1:int(A[1])-delete])
        else:
            return (CG[int(A[0])-1:int(C[0])-2]+CG[int(C[1])+2:int(A[1])-delete])

def trim_E7 (genes):
    if E6 in genes:
        A = genes['E6'].split("..")
        Ar = range(int(A[0]),int(A[1]))
    else:
        Ar=[0]
        
        
    B = genes['E7'].split("..")
    C = genes['E1'].split('..')
    CG = genes['CG']
    
    

    Br = range(int(B[0]),int(B[1]))
    Cr = range(int(C[0]),int(C[1]))
    
    overlap5 = set(Ar).intersection(set(Br))
    overlap3 = set(Br).intersection(set(Cr))
    if len(overlap5) % 3 == 0:
        delete = len(overlap5)+2
        if len(overlap3) == 0:
            return (CG[int(B[0])+delete:int(B[1])-len(overlap3)-3])         
        else:
            return (CG[int(B[0])+delete:int(B[1])-len(overlap3)-2])
    elif len(overlap5) % 3 == 1:
        delete = len(overlap5)+1
        if len(overlap3) == 0:
            return (CG[int(B[0])+delete:int(B[1])-len(overlap3)-3])         
        else:
            return (CG[int(B[0])+delete:int(B[1])-len(overlap3)-2])

def trim_E1(genes):

    if E7 in genes:
        A = genes['E7'].split("..")
        Ar = range(int(A[0]),int(A[1]))
    else:
        Ar=[0]
    B = genes['E1'].split("..")
    C = genes['E2'].split('..')
    if "E8^E2" in genes:
        m=re.search("(\d*)\.\.(\d*)\+(\d*)\.\.(\d*)", genes['E8^E2'])
        if m:
            D = (m.groups()[:2])
    else:
        D = None



    CG = genes['CG']
    

    Br = range(int(B[0]),int(B[1]))
    Cr = range(int(C[0]),int(C[1]))
    
    

    
    overlap5 = set(Ar).intersection(set(Br))
    overlap3 = set(Br).intersection(set(Cr))
    
    if D == None:
        if len(overlap5) == 0:
            return (CG[int(B[0])+len(overlap5)-1:int(B[1])-len(overlap3)-2])
        else: 
            return (CG[int(B[0])+len(overlap5)+1:int(B[1])-len(overlap3)-2])
    else:
        if len(overlap5) == 0:
            return (CG[int(B[0])+len(overlap5)-1:int(D[0])-2]+CG[int(D[1]):int(B[1])-len(overlap3)-2])
        else: 
            return (CG[int(B[0])+len(overlap5)+1:int(D[0])-2]+CG[int(D[1]):int(B[1])-len(overlap3)-2])

def trim_E2(genes):
    A = genes['E1'].split("..")
    B = genes['E2'].split("..")
    m=re.search("(\d*)\.\.(\d*)\+(\d*)\.\.(\d*)", genes['E1^E4'])
    if m:
        D = (m.groups()[2:])
    CG = genes['CG']
    
    Ar = range(int(A[0]),int(A[1]))
    Br = range(int(B[0]),int(B[1]))

    
    overlap5 = set(Ar).intersection(set(Br))

    gene = CG[int(B[0])+len(overlap5)+1:int(D[0])-3]+CG[int(D[1])+2:int(B[1])]
    return (gene[:-3])

def trim_L2 (genes):
    A = genes['L2'].split("..")
    B = genes['L1'].split("..")
    CG = genes['CG']
    
    Ar = range(int(A[0]),int(A[1]))
    Br = range(int(B[0]),int(B[1]))
    
    overlap = (set(Ar).intersection(set(Br)))
    if len(overlap) % 3 == 0:
        delete = len(overlap)+3

        return (CG[int(A[0])-1:int(A[1])-delete])
    elif len(overlap) % 3 == 1:
        delete = len(overlap)+2
        return (CG[int(A[0])-1:int(A[1])-delete])

def trim_L1 (genes):
    A = genes['L2'].split("..")
    B = genes['L1'].split("..")

    CG = genes['CG']
    
    Ar = range(int(A[0]),int(A[1]))
    Br = range(int(B[0]),int(B[1]))

    
    overlap5 = set(Ar).intersection(set(Br))
    if (len(overlap5)) == 0:
        gene = CG[int(B[0])+len(overlap5)-1:int(B[1])]

    else:
        gene = CG[int(B[0])+len(overlap5)+1:int(B[1])]
    return (gene[:-3])




viruses = ["TbelPV1","HPV204","HPV1","HPV63","CcanPV1","OcPV1","SfPV1","TePV1","CPV1","LwiePV1","CcrPV1","FcaPV1","PlpPV1","UuPV1","LrPV1","PcPV1","CPV6","AmPV4","LwPV1","ElPV1","PlPV1","MscPV2","RfPV1","EsPV1","EsPV3","EhPV1","TbraPV1","TbraPV2","TbraPV3","EdPV1","HPV41"]
accessions = {}
for v in viruses:
    for r in csv.reader(open('../shared_data/reference.csv',"r")):
        if r[0] == v:
            accessions[r[1]] = v



with open("E6.fas","w") as E6, open("E7.fas","w") as E7, open("E1.fas","w") as E1, open("E2.fas","w") as E2, open("L2.fas","w") as L2, open("L1.fas","w") as L1, open("concatenated.fas","w") as conc:

    for a in accessions:
        concatenated = []
        print (accessions[a])
        genes = {}
        for r in csv.reader(open("../shared_data/PaVE.csv","r")):
            if a == r[0]:
                if r[1] == "CG":
                    genes[r[1]] = r[3]
                else:
                    genes[r[1]] = r[2]
    
        if "E6" in genes:
            print (f">{accessions[a]}\n{trim_E6(genes)}", file= E6)
            concatenated.append(trim_E6(genes))
        if "E7" in genes:
            print (f">{accessions[a]}\n{trim_E7(genes)}", file= E7)
            concatenated.append(trim_E7(genes))
        if "E1" in genes:
            print (f">{accessions[a]}\n{trim_E1(genes)}", file= E1)
            concatenated.append(trim_E1(genes))
        if "E2" in genes:      
            print (f">{accessions[a]}\n{trim_E2(genes)}", file= E2)
            concatenated.append(trim_E2(genes))
        if "L2" in genes:       
            print (f">{accessions[a]}\n{trim_L2(genes)}", file= L2)
            concatenated.append(trim_L2(genes))
        if "L1" in genes:       
            print (f">{accessions[a]}\n{trim_L1(genes)}", file= L1)
            concatenated.append(trim_L1(genes))
        print (f">{accessions[a]}\n{''.join(concatenated)}", file=conc)
    
