import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir,"..",'aptafast'))
sys.path.append(root_path)

import AptaFast as AF

f=open("data_clustering_test.dat",'a')
f.write("FAMILY    PDB_chain    SEQUENCE    DOTBRACKET\n")
for i in range(1,30):
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


for i in range(1,30):
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


for i in range(1,30):
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

for i in range(1,30):
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

