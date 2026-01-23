import os
from tqdm import tqdm

path = "/home/bcuvillier/Documents/Datasets/RNAStrAlign_bpseq"

def create_file_list(path):
    """Parcourt récursivement une arborescence et affiche le chemin de chaque fichier et répertoire."""
    all_files=[]
    for repertoire_courant, sous_repertoires, fichiers in os.walk(path):
        for nom in fichiers:
            if nom != ".DS_Store":
                all_files.append(os.path.join(repertoire_courant, nom))
    
    return all_files
            
all_files = create_file_list(path)

def dotbracket_verif(input_str):
    """
        This function verifies if all the opened brackets are actually closed.
        It appears that sometimes, the output of random_gen is invalid, we want to ignore those outputs.
    """
    nb=0
    for elt in input_str:
        if elt=="(":
            nb+=1
        if elt==")":
            nb-=1
    
    add=0
    for elt in input_str:
        if elt==".":
            add+=1
    if add == len(input_str):
        nb=1
    
    return not bool(nb)

def file_converter(path):
    with open(path) as bp:
        lines=bp.readlines()
    dbnot=["." for i in range(len(lines))]
    ACTG=""
    done=[]
    ignore=False
    for line in lines:
        line = line.replace("\n","")
        letter=""
        bp1=""
        bp2=""
        i=0
        f_bp1=True
        f_letter=False
        while i != len(line):
            if f_bp1:
                if line[i]==" ":
                    f_bp1=False
                    f_letter=True
                else:
                    bp1+=line[i]
            elif f_letter:
                if line[i]==" ":
                    f_letter=False
                else:
                    letter+=line[i]
            else:
                bp2+=line[i]
            i+=1
        
        ACTG+=letter
        if (int(bp1)!=0 and int(bp2)!=0) and (int(bp1) not in done) and (int(bp2) not in done) and not ignore:
            done.append(int(bp1))
            done.append(int(bp2))
            if int(bp1)>len(dbnot) or int(bp2)>len(dbnot):
                print("invalid bpseq file, ignoring:",path)
                ignore=True
            if not ignore:
                dbnot[int(bp1)-1]="("
                dbnot[int(bp2)-1]=")"
           
    if ignore :
        return False, path
    
    db=""
    for elt in dbnot:
        db+=elt
    
    if not dotbracket_verif(db):
        print("Invalid dotbracket, ignoring:",path)
        return False, path
    
    doss=path.split("/")
    name=doss[-1].replace(".bpseq","")
    family=doss[6]
    
    return True, name, family, ACTG, db
    

    

def parser(all_files):
    all_structs=[]
    ignored=[]
    for paths in tqdm(all_files):
        struct = file_converter(paths)
        if struct[0]:
            all_structs.append([struct[1],struct[2],struct[3],struct[4]])
        else:
            ignored.append(struct[1])
        
    return all_structs, ignored


all_structs, ignored=parser(all_files)

def create_CLUSTER(all_structs,nb_max_per_fam):
    f = open("/home/bcuvillier/Documents/AptaMat/clustering/RNAstralign_clust_reduced.dat",'a')
    tbw=""
    dict_fam={}
    for struct in all_structs:
        family=struct[1]
        if family not in dict_fam.keys():
            dict_fam[family]=1
        else:
            dict_fam[family]+=1
        
        if dict_fam[family] <= nb_max_per_fam:
            name=struct[0]
            sequence=struct[2]
            dotbracket=struct[3]
            tbw+=family+"    "+name+"    "+sequence+"    "+dotbracket+"\n"

    f.write("FAMILY    NAME    SEQUENCE    DOTBRACKET\n")   
    f.write(tbw)
    f.close()
    
create_CLUSTER(all_structs,1000)