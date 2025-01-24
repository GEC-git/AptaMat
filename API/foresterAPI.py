import os
import subprocess

def forester_pairwise(struct1,struct2):

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
    
    ident1, dotbracket1al, ident2, dotbracket2al = forester_from_file(filename)
    
    
    os.remove(filename)
    
    return dotbracket1al, dotbracket2al

def forester_from_file(file):
    lines_from=open(file,'r').read()
    result=subprocess.run(["RNAforester","--fasta"],input=lines_from+"\n&",capture_output=True,text=True)
    lines=result.stdout.splitlines()
    return lines[14].replace("\n",""),lines[16].replace("\n",""),lines[17].replace("\n",""),lines[19].replace("\n","")

