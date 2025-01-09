import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir,"..","..",'aptafast'))
sys.path.append(root_path)

import AptaFast as AF

def fam_compteur_fa(file_from):
    lines=open(file_from).readlines()
    dico={}
    for i, line in enumerate(lines):
        if line.startswith(">"):
            line1=line.strip(">\n")
            line1=line1.replace(" ","_")
            line1=line1.strip("0123456789")
            if line1 not in dico:
                dico[line1]=1
            else:
                dico[line1]+=1
    return dico


def file_converter_FASTA_to_CLUSTER():
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_RNALIGNED_150x8_CLEANED.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED/data_clustering_BPRNA_RFAM_scraped_10000_FASTA_NONALIGNED_8x150_CLEANED.fa"
    lines=open(file_from).readlines()
    tbw=""
    f.write("FAMILY    dbn    SEQUENCE    DOTBRACKET\n")
    for i, line in enumerate(lines):
        if line.startswith(">"):
            line1=line.strip(">\n")
            line1=line1.replace(" ","_")
            tbw+=line1.strip("0123456789")+"    "+line1+".dbn    "
        elif not AF.Dotbracket.is_dotbracket(line):
            tbw+=line.strip("\n")+"    "
        else:
            tbw+=line
    f.write(tbw)
    f.close()


def find_spaces(string):
    index=[]
    for i,elt in enumerate(list(string)):
        if elt==" ":
            index.append(i)
    return index

def file_converter_FASTA_to_CLUSTER_choices(nb_fam_max,nb_max_per_fam):
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/datasets/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_NONALIGNED_CHOICES.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/datasets/data_clustering_BPRNA_RFAM_scraped_10000.fa"
    lines=open(file_from).readlines()
    tbw=""
    f.write("FAMILY    dbn    SEQUENCE    DOTBRACKET\n")
    
    fam=fam_compteur_fa(file_from)
    dico_inter={}
    for i, line in enumerate(lines):
        if nb_fam_max == 0:
            f.write(tbw)
            f.close()
        else:
            if line.startswith(">"):
                line1=line.strip(">\n")
                line1=line1.replace(" ","_")
                line1_fam=line1.strip("0123456789")
                if line1_fam not in dico_inter:
                    if fam[line1_fam]>nb_max_per_fam:
                        dico_inter[line1_fam]=1
                        nb_fam_max-=1
                        allowed=True
                        tbw+=line1_fam+"    "+line1+".dbn    "
                    else:
                        allowed = False
                else:
                    if fam[line1_fam]>nb_max_per_fam:
                        if dico_inter[line1_fam]<nb_max_per_fam:
                            dico_inter[line1_fam]+=1
                            allowed=True
                            tbw+=line1_fam+"    "+line1+".dbn    "
                        else:
                            allowed = False
                    else:
                        allowed = False
                
            elif not AF.Dotbracket.is_dotbracket(line) and allowed:
                tbw+=line.strip("\n")+"    "
            elif allowed:
                tbw+=line

           
def file_converter_CLUSTER_to_FASTA():
    f_created=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/data_clustering_BPRNA_RFAM_scraped_10000_FASTA_NONALIGNED_8x150.fa",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_NONALIGNED_8x150.dat"
    lines=open(file_from).readlines()
    tbw=""
    for i, line in enumerate(lines):
        if line.startswith("FAMILY"):
            pass
        else:
            spaces=find_spaces(line)
            family=""
            name=""
            sequence=""
            dotbracket=""
            nb_space=0
            for j,elt in enumerate(list(line)):
                if j in spaces:
                    if list(line)[j-1] != " ":
                        nb_space+=1
                elif nb_space==0:
                    family+=elt
                elif nb_space==1:
                    name+=elt
                elif nb_space==2:
                    sequence+=elt
                elif nb_space==3:
                    dotbracket+=elt
            name=name.strip(".dbn")
            tbw+=">"+family+str(i)+"\n"+sequence+"\n"+dotbracket
    f_created.write(tbw)
    f_created.close()
            
    
    
    
    
    
    
    
    
    