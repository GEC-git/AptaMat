import os
import subprocess



def forester_pairwise(struct1,struct2):

    tbw=">"+struct1.id+"\n"
    tbw+=struct1.sequence+"\n"
    tbw+=struct1.dotbracket+"\n"
    tbw+=">"+struct2.id+"\n"
    tbw+=struct2.dotbracket+"\n"
    tbw+=struct2.dotbracket
    f_created=open("temp_struct"+struct2.id+str(os.getpid())+".fa",'a')
    f_created.write(tbw)
    f_created.close()
    
    ident1, dotbracket1al, ident2, dotbracket2al = forester_from_file("temp_struct"+struct2.id+str(os.getpid())+".fa")
    
    os.remove("temp_struct1"+struct1.id+str(os.getpid())+".fa")
    os.remove("temp_struct2"+struct2.id+str(os.getpid())+".fa")
    
    return dotbracket1al, dotbracket2al

def forester_from_file(file):
    result=subprocess.run(["RNAforester","-f="+file,"fasta"],capture_output=True,text=True)
    lines=result.stdout.splitlines()
    print(result.stdout)
    #print(lines)