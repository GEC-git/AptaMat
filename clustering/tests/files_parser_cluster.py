import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir,"..",'aptafast'))
sys.path.append(root_path)

import AptaFast as AF

def file_parser_clustering():
    f=open("data_clustering_test_150x10_TBA_RNAlign.txt",'a')
    f.write("FAMILY    PDB_chain    SEQUENCE    DOTBRACKET\n")
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_tmRNA_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=str(line1)
                tbw+="    "
                tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line.strip("\n")
        tbw+="\n"
        f.write(tbw)        
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_CRW_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=str(line1)
                tbw+="    "
                tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line.strip("\n")
        tbw+="\n"
        f.write(tbw)
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_RFAM_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=str(line1)
                tbw+="    "
                tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line.strip("\n")
        tbw+="\n"
        f.write(tbw)
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_SRP_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=str(line1)
                tbw+="    "
                tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line.strip("\n")
        tbw+="\n"
        f.write(tbw)
    f.close()

def file_parser_RNAlign():
    f=open("data_clustering_test_150x10_TBA_RNAlign.txt",'a')
    #f.write("FAMILY    PDB_chain    SEQUENCE    DOTBRACKET\n")
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_tmRNA_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=">"+str(line1)+str(i)+"\n"
                #tbw+="    "
                #tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    #tbw+=line.strip("\n")+"    "
                    tbw+=line
                else:
                    #tbw+=line.strip("\n")
                    tbw+=line
        #tbw+="\n"
        f.write(tbw)        
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_CRW_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=">"+str(line1)+str(i)+"\n"
                #tbw+="    "
                #tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    #tbw+=line.strip("\n")+"    "
                    tbw+=line
                else:
                    #tbw+=line.strip("\n")
                    tbw+=line
        #tbw+="\n"
        f.write(tbw)
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_RFAM_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=">"+str(line1)+str(i)+"\n"
                #tbw+="    "
                #tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line
                    #tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line
                    #tbw+=line.strip("\n")
        #tbw+="\n"
        f.write(tbw)
    
    for i in range(1,151):
        file = "/home/bcuvillier/Téléchargements/dbnFiles/bpRNA_SRP_"+str(i)+".dbn"
        lines = open(file).readlines()
        tbw=""
        for j, line in enumerate(lines):
            if j == 0:
                line1=line.strip("#Name: _"+str(i)+"\n")
                line2=line.strip("#Name: \n")
                tbw+=">"+str(line1)+str(i)+"\n"
                #tbw+="    "
                #tbw+=str(line2)+".dbn    "
            elif line.startswith("#"):
                pass
            else:
                if not AF.Dotbracket.is_dotbracket(line):
                    tbw+=line
                    #tbw+=line.strip("\n")+"    "
                else:
                    tbw+=line
                    #tbw+=line.strip("\n")
        #tbw+="\n"
        f.write(tbw)
    f.close()    

def file_converter_FASTA_to_CLUSTER():
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/first_tests/clustering_dataset_CLUSTER_ALIGNED_pseudomode.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/first_tests/clustering_dataset_FASTA_ALIGNED_pseudomode.dat"
    lines=open(file_from).readlines()
    tbw=""
    f.write("FAMILY    dbn    SEQUENCE    DOTBRACKET\n")
    for i, line in enumerate(lines):
        if line.startswith(">"):
            line1=line.strip(">\n")
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

           
def file_converter_CLUSTER_to_FASTA():
    f_created=open("data_clustering_test_29x4_FASTA_NON_ALIGNED.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/first tests/clustering_dataset_NON_ALIGNED.dat"
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
            tbw+=">"+name+"\n"+sequence+"\n"+dotbracket
    f_created.write(tbw)
            
    
    
    
    
    
    
    
    
    