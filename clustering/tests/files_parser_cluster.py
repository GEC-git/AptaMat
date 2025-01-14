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

def pseudoknots_annihilator():
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED + NOPSEUDOKNOTS/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_NONALIGNED_150x8_ULTRACLEANED.dat"
    f_created=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/data_clustering_BPRNA_RFAM_scraped_10000_FASTA_NONALIGNED_150x8_ULTRACLEANED_NOPSEUDOKNOTS.fa",'a')
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
            ignore=False
            for elt in dotbracket:
                if elt in "[]{}<>123456789":
                    ignore=True
            if not(ignore):
                tbw+=">"+family+str(i)+"\n"+sequence+"\n"+dotbracket
    f_created.write(tbw)
    f_created.close()

def file_converter_BEAGLEexport_to_cluster():
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_150x8_ULTRACLEANED_NOPSEUDOKNOTS_REDUCED_beagle2aligned.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/export_beagle.txt"
    lines=open(file_from).readlines()
    score=0
    score_after=0
    results_dictionnary={}
    for i,line in enumerate(lines):
        if line.startswith(">"):
            score=score_after
            new_line=True
            line=line.strip(">")
            j=0
            id=""
            while line[j]!="|":
                id+=str(line[j])
                j+=1
            if id not in results_dictionnary.keys():
                results_dictionnary[id]=[None,None]
                score=0
                score_after=0
            score_after=float(line[line.rfind("Z-score:")+8:])
            
        elif score_after>=score and new_line:
            new_line=False
            results_dictionnary[id][0]=line.strip("\n")
            line_after=True
        elif score_after>=score and not new_line and line_after:
            results_dictionnary[id][1]=line.strip("\n")
            line_after=False
    
    tbw=""
    for items in results_dictionnary.items():
        family=items[0][0:9]
        dbn=items[0]
        sequence=items[1][0]
        dotbracket=items[1][1]
        tbw+=family+"    "+dbn+"    "+sequence+"    "+dotbracket+"\n"
    
    f.write("FAMILY    dbn    SEQUENCE    DOTBRACKET\n")   
    f.write(tbw)
    f.close()

def file_converter_PP_to_CLUSTER():
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_150x8_ULTRACLEANED_NOPSEUDOKNOTS_REDUCED_locarnaligned.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/data_clustering_BPRNA_RFAM_scraped_10000_FASTA_NONALIGNED_150x8_ULTRACLEANED_NOPSEUDOKNOTS_REDUCED.out/results/result.pp"
    lines=open(file_from).readlines()
    tbw=""
    f.write("FAMILY    dbn    SEQUENCE    DOTBRACKET\n")
    for i, line in enumerate(lines):
        line=line.replace("\n","")
        if not(line.startswith("#")):
            family=line[0:9]
            dbn=line[0:13].replace(" ","")
            sequence=False
            dotbracket=""
            ACTG=""
            first=True
            for k in range(13,len(line)):
                j=line[k]
                if j!=" ":
                    if j in '.()[]' and first:
                        sequence=True
                        first=False
                    
                    if sequence:
                        dotbracket+=j
                    else:
                        ACTG+=j
            tbw+=family+"    "+dbn+"    "+ACTG+"    "+dotbracket+"\n"
    f.write(tbw)
    f.close()
    
def file_converter_FASTA_to_CLUSTER():
    f=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_150x8_ULTRACLEANED_NOPSEUDOKNOTS_REDUCED_foresteraligned.dat",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED+NOPSEUDOKNOTS+REDUCED/data_clustering_BPRNA_RFAM_scraped_10000_FASTA_150x8_ULTRACLEANED_NOPSEUDOKNOTS_REDUCED_foresteraligned.fa"
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
    f_created=open("/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/fasta_test.fa",'a')
    file_from="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_LOCARNALIGNED_150x8_ULTRACLEANED.dat"
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
            
