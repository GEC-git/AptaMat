def file_reader():
    dot="(.-"
    file_from="" #PLEASE PUT HERE PATH TO THE EXPORTED .TXT FROM BEAGLE2.
    lines=open(file_from).readlines()
    results_dictionnary={}
    for i,line in enumerate(lines):
        if line.startswith(">"):
            second=False
            first=True
            line=line.replace(">","")
            id1=""
            j=0
            while line[j]!="|":
                id1+=line[j]
                j+=1
            id2=""
            j+=1
            while line[j]!="|":
                id2+=line[j]
                j+=1
            results_dictionnary[id1+id2]=[None,None]
            
        elif line[0] in dot and "A" not in line and "C" not in line and "T" not in line and "G" not in line and "U" not in line and ":" not in line:
            if first:
                results_dictionnary[id1+id2][0]=line.replace("\n","")
                first=False
                second=True
            elif second:
                results_dictionnary[id1+id2][1]=line.replace("\n","")
                second=False
                
    return results_dictionnary

results_dictionnary=file_reader()

def get_db(struct1,struct2):
    if struct1.id==struct2.id:
        return struct1.dotbracket, struct2.dotbracket
    try:
        db1=results_dictionnary[struct1.id+struct2.id][0]
        db2=results_dictionnary[struct1.id+struct2.id][1]
    except KeyError:
        db2=results_dictionnary[struct2.id+struct1.id][0]
        db1=results_dictionnary[struct2.id+struct1.id][1]

    return db1, db2