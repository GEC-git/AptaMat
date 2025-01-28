import os
import subprocess

def rnalign2d_pairwise(struct1,struct2):

    tbw=">"+struct1.id+"\n"
    tbw+=struct1.sequence+"\n"
    tbw+=struct1.dotbracket+"\n"
    tbw+=">"+struct2.id+"\n"
    tbw+=struct2.sequence+"\n"
    tbw+=struct2.dotbracket
    filename="temp_struct"+struct2.id+struct1.id+str(os.getpid())+".fa"
    f_created=open(filename,'w+')
    f_created.write(tbw)
    f_created.close()
    
    dotbracket1al, dotbracket2al = align_from_file(filename)
    
    
    os.remove(filename)
    
    return dotbracket1al, dotbracket2al

def align_from_file(file):
    result=subprocess.run(["rnalign2d","-i",file,"-o",file.replace(".fa","")+"_output.fa"],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    lines = open(file.replace(".fa","")+"_output.fa").readlines()
    os.remove(file.replace(".fa","")+"_output.fa")

    return lines[2].replace("\n",""), lines[5].replace("\n","")

